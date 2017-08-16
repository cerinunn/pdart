#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Make 'chains' by checking the framecount and inserting gaps.
Where the framecounts are consecutive and the timestamps are within
a tolerance, the method will make a chain.
Initially, the method will check for a perfect match - for 4 consecutive
framecounts. If these are not found with the first 4 framecounts, it will
loop through until it finds 4 suitable records.
Gaps of less than 3 samples will be included in the chain.
The starttime is the first timestamp in the chain.
More than 4 consecutive mismatches, mean that the chain will be broken,
and a new chain started if possible.

:copyright:
    The pdart Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from datetime import timedelta
import os
import io
import gzip
import glob
import numpy as np
import numpy.ma as ma
import logging, logging.handlers
import csv
import fnmatch
import shutil

from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, Trace, Stats, read
from pdart.util import stream_select, trace_eq

global DELTA
DELTA = 0.15094
SECONDS_PER_DAY = 3600.0 * 24.0
INVALID = -99999
ABSOLUTE_ADJUST_TIME = 2.
MIN_SAMPLE_LENGTH = 361

def call_build_chains(
    stations=['S11','S12','S14','S15','S16'],
    starttime0=None,
    framecount0=None,
    adjust0=None,
    obs_delta0=DELTA,
    timestamp0=None,
    framecount_adjust=None,
    raw_dir='.',
    chain_dir='.',
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T:21:00.000000Z'),
    read_gzip=True,
    write_gzip=True):
    '''
    Calls build_chains()
    '''

    log_filename = 'logs/build_chains.log'
    logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
    # logging.basicConfig(filename=log_filename, filemode='w', level=logging.DEBUG)

    for station in stations:

        # check that the overall directory exists
        if not os.path.exists(chain_dir):
            msg = ("The directory {} doesn't exist".format(chain_dir))
            raise IOError(msg)
        else:
            # make the subdirectory with the station name
            chain_dir_station =  os.path.join(chain_dir, station)
            if not os.path.exists(chain_dir_station):
                os.makedirs(chain_dir_station)

    # build chains for each station
    for station in stations:

        raw_dir_station =  os.path.join(raw_dir, station)
        chain_dir_station =  os.path.join(chain_dir, station)

        time_interval = timedelta(hours=3)
        start = start_time
        while start < end_time:

            # work out the filenames
            stream_filename = '%s_%s.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
            if read_gzip:
                stream_filename = '%s.gz' % (stream_filename)
            chain_filename = '%s_%s.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
            gzip_chain_filename = '%s.gz' % (chain_filename)
            raw_dir_filename = os.path.join(raw_dir_station, stream_filename)
            chain_dir_filename = os.path.join(chain_dir_station, chain_filename)
            gzip_chain_dir_filename = os.path.join(chain_dir_station,
              gzip_chain_filename)

            # read in the raw SEED file
            try:
                stream = read(raw_dir_filename)
            except FileNotFoundError:
                msg = 'chains.py cannot find file: {}'.format(raw_dir_filename)
                print(msg)
                logging.info(msg)
                # increment the time interval
                start += time_interval
                continue

            # select for just this station, (not necessary, but just in case)
            stream = stream.select(station=station)


            if len(stream) > 0:

                # build the chains (for this station)
                stream1 = build_chains(stream=stream)
                # get rid of any short traces that overlap with other
                # streams
                stream2 = discard_short_traces(stream1)

                if len(stream2) > 0:
                    # split if the streams have masked data (because there
                    # are gaps
                    stream2 = stream2.split()
                    stream2.write(chain_dir_filename, 'MSEED')

                    if write_gzip:
                        with open(chain_dir_filename, 'rb') as f_in, gzip.open(gzip_chain_dir_filename, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                        os.unlink(chain_dir_filename)

            # increment the time interval
            start += time_interval

def build_chains(stream):
    '''
    Make 'chains' by checking the framecount and inserting gaps.
    '''
    log_filename = 'logs/build_chains.log'
    logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
    # logging.basicConfig(filename=log_filename, filemode='w', level=logging.DEBUG)

    # quick check to make sure only one station
    station = stream[0].stats.station

    if len(stream) != len(stream.select(station=station)):
        raise ValueError("More than one station in the stream")

    # TODO only need 3 if running SPZ
    # original_data = np.full((n,3), INVALID, 'int32')

    # begin by selecting the raw trace, and sorting by starttime
    FR_stream = stream.select(channel='_FR')
    FR_stream = FR_stream.sort(keys=['starttime'])

    return_stream = Stream()

    # for each of the raw streams
    for fs in FR_stream:

        # find the matching streams with the same start and end time
        original = stream_select(stream,network=fs.stats.network, station=fs.stats.station,
          location=fs.stats.location,starttime=fs.stats.starttime,endtime=fs.stats.endtime)

        # get the stream for the timing trace
        TT_stream = original.select(channel='_TT')

        start_pointer = 0
        pointer = 0
        consecutive_invalid = 0
        valid_chain = False
        len_data = len(fs.data)

        # TODO remove invalid nasty data

        while start_pointer < len_data:
            # loop through data from start pointer to the end
            msg = ('Valid chain: {} Start pointer{} {}'.format(valid_chain, start_pointer, fs.data[start_pointer:]))
            logging.debug(msg)
            for i, framecount1 in enumerate(fs.data[start_pointer:]):
                # pointer is the CURRENT index
                pointer = start_pointer + i

                match_timestamp1 = TT_stream[0].data[pointer]

                # first step - look for a short chain assuming the first one in the trace is ok
                if i == 0:
                    msg = ('i = 0, {} {} {} {}'.format(i, start_pointer, framecount1, UTCDateTime(match_timestamp1)))
                    logging.debug(msg)
                    chain_framecount0 = fs.data[start_pointer]
                    chain_timestamp0 = TT_stream[0].data[start_pointer]
                    chain_pointer0 = start_pointer
                    # make a pointer_array
                    n = fs.stats.npts
                    pointer_array = np.full(n, INVALID, 'int32')
                    pointer_array[start_pointer] = start_pointer
                    valid_chain = False
                    consecutive_invalid = 0
                    # if the framecount is out of range, continue
                    if framecount1 < 0 or framecount1 > 89.75:
                        msg = 'invalid framecount'
                        logging.debug(msg)
                        break

                else: # records where i < 0

                    # use for debugging
                    # if i < 6:
                    #     msg = ('i = {} {} {}'.format(i, start_pointer, framecount1))
                    #     logging.debug(msg)

                    # if the framecount is out of range, start again
                    if framecount1 < 0 or framecount1 > 89.75:
                        if i < 4:
                            # unable to make a chain of 4, so break out
                            valid_chain = False
                            msg = 'invalid framecount, less than 4 {}'.format(start_pointer)
                            logging.debug(msg)
                            break
                        else:
                            # we just ignore it
                            continue
                    # if the frame range is valid
                    else:
                        # check for the correct sample index from the framecount
                        # and timestamp
                        sample_diff = _calc_match(match_timestamp1, chain_timestamp0, chain_framecount0, framecount1, obs_delta0=DELTA)

                        if sample_diff is not None:
                            pointer_array[pointer] = sample_diff + chain_pointer0
                            # last_idx = pointer_array[pointer-1]
                            msg = ('Sample, i, framecount and pointer array, ', sample_diff, i, framecount1, str(pointer_array[0:7]))
                            logging.debug(msg)

                            # check for consecutive values within the frame
                            if pointer_array[pointer-1]+1 != pointer_array[pointer]:
                                consecutive_invalid += 1
                            else:
                                consecutive_invalid = 0

                            # the current one is valid, so update the
                            # timestamp and framecount
                            chain_timestamp0 = match_timestamp1
                            chain_framecount0 = framecount1
                            chain_pointer0 += sample_diff
                        else:
                            consecutive_invalid += 1


                        # if i is 3 and the framecounts have been consecutive
                        # then mark the chain as invalid
                        if i == 3 and consecutive_invalid == 0:
                            valid_chain = True

                        if i < 4 and consecutive_invalid != 0:
                            # unable to make a chain of 4, so break out
                            msg = ('Unable to make a chain of 4, so break out {} {} {}'.format(sample_diff, pointer, framecount1))
                            logging.debug(msg)
                            break

                        # break the chain if more than 3 framecounts have been
                        # invalid
                        if consecutive_invalid > 3:
                            break


            # make a chain from previous execution the loop if a
            # valid one exists
            if valid_chain:
                re_stream, last_pointer = _reconstruct_streams(original, pointer_array)
                return_stream += re_stream
                start_pointer = last_pointer + 1
                valid_chain = False
                msg = ('Start pointer after valid {} {} {}'.format(start_pointer, last_pointer, len_data))
                logging.debug(msg)
            else:
                start_pointer = start_pointer + 1
                msg = ('Start pointer after invalid {} {}'.format(start_pointer, len_data))
                logging.debug(msg)

    return return_stream


def _samples_to_shift(framecount1, est_framecount1):
    """
    Find number of samples to shift - used to make chains and in splicing.
    We pass the framecount of the trace (framecount1). An estimate of the
    correct framecount is passed in with est_framecount1.
    The function returns the number of samples to shift.

    Example:
    3.00 3.25 3.50 3.75 4.00 4.25  the estimated framecount is 3.75
    xxxx xxxx xxxx 4.00 4.25 4.50  framecount1 is 4
    so we have to shift the second trace one sample (positive idx_diff)
    xxxx xxxx xxxx xxxx 4.00 4.25 4.50
    """

    if framecount1 < 0 or framecount1 > 89.75:
        raise ValueError('Framecount is invalid: %s' % framecount1)

    frm_diff1 = framecount1 - est_framecount1

    # if they are the same, no shift is required
    if frm_diff1 == 0:
        idx_diff = 0
        return idx_diff

    # we need to count forwards and backwards, and adjust for the end of the
    # frame
    if frm_diff1 > 0:
        frm_diff2 = frm_diff1 - 90
    elif frm_diff1 < 0:
        frm_diff2 = frm_diff1 + 90

    # find the smallest frame difference (postive or negative )
    if abs(frm_diff1) < abs(frm_diff2):
        frm_diff = frm_diff1
    else:
        frm_diff = frm_diff2

    # change to sample numbers (4 samples per frame)
    idx_diff = frm_diff * 4

    return idx_diff

def _calc_match_samp_frame(starttime1, starttime0, framecount0, framecount1, obs_delta0):
    """
    Calculate the estimated index number and frame number.
    The estimate is based on matching the timestamp of the trace we
    are editing (starttime1), with the starttime and framecount of the first one in
    the group. The method calculates the sample index and the framecount.
    """
    tim_diff = starttime1 - starttime0
    if tim_diff < 0 :
        raise ValueError('Time difference is negative.')
    sample_idx = int(round(tim_diff/obs_delta0))

    # there are 90*4 = 360 samples per block
    # we can ignore whole blocks so just use the modulo operator
    part_block_samples = sample_idx % 360
    # round to nearest 0.25
    part_block_samples = (round(part_block_samples*4)/4)
    # change to framecounts (4 samples per frame number)
    part_block = part_block_samples/4
    # calculate the correct framecount for the time of the
    # new trace
    est_framecount1 = framecount0 + part_block
    # if the number is greater than 89.75, take 90 from it to get the correct number
    if est_framecount1 > 89.75:
        est_framecount1 = est_framecount1 - 90

    if framecount1 != est_framecount1:
        # if the trace needs shifting sample numbers, do it here:
        idx_diff = _samples_to_shift(framecount1, est_framecount1)
        if abs(idx_diff) < 12:
            sample_idx = sample_idx + idx_diff
            # make some more checks
            if sample_idx != 0:
                observed_delta = tim_diff / sample_idx
                # obs_delta_percent = round(100.*(DELTA-observed_delta)/DELTA,3)
                obs_delta_percent = 100.*(DELTA-observed_delta)/DELTA
                msg = 'Average sampling interval:{0} Percent error from nominal interval:{1}%'.format(observed_delta, obs_delta_percent)
                logging.info(msg)
                if abs(obs_delta_percent) > 1:
                    msg = 'Reject automatically'
                    logging.debug(msg)
                    sample_idx = None
            else:
                sample_idx = None
        else:
            sample_idx = None

    return sample_idx



def discard_short_traces(stream):
    """
    Discard short traces which overlap with longer traces.
    Short traces are often quite poor quality - but sometimes occur when
    a longer trace exists. If so, they can be safely discarded.
    """

    # copy the original stream
    return_stream = stream.copy()

    # sort a stream (from the timing channel) with the number of samples
    sorted_stream = return_stream.select(channel='ATT').sort(keys=['npts'])

    # if there is more than one trace in sorted_stream, see if there are any
    # traces to discard.
    if len(sorted_stream) > 1:

        # outer loop of traces, sorted number of samples
        for tr in sorted_stream:
            # if the trace is short
            if tr.stats.npts < MIN_SAMPLE_LENGTH:
                start_timestamp = tr.data[0]
                end_timestamp = tr.data[-1]
                # inner loop of traces, to check against
                for tr1 in sorted_stream:
                    remove_flag = False
                    # if the inner and outer trace are the same, do nothing
                    if trace_eq(tr,tr1):
                        continue
                    start_timestamp_check = tr1.data[0]
                    end_timestamp_check = tr1.data[-1]
                    # check the short trace overlaps both ends of another trace
                    if ( start_timestamp > start_timestamp_check and
                      end_timestamp < end_timestamp_check ):
                        remove_flag = True
                        msg = ('Removing short trace: ', tr)
                        logging.debug(msg)
                        stream_short = stream_select(sorted_stream,network=tr.stats.network, station=tr.stats.station,
                          location=tr.stats.location,starttime=tr.stats.starttime,endtime=tr.stats.endtime)
                        for tr2 in stream_short:
                            # remove from the return_stream
                            return_stream.remove(tr2)

                        if remove_flag:
                            # break the inner loop (and continue the outer one)
                            break

                if remove_flag:
                    # if we removed the trace, we can move to the next short sample
                    continue

            # the stream is ordered by trace length, so we can stop execution
            # when the traces are too long
            else:
                break

    return return_stream

def _reconstruct_streams(stream, pointer_array):
    """
    Reconstruct streams based on the correct positioning.
    The positioning comes from the pointer array. For example,

    original:
    012367

    pointer_array:
    [0,1,2,3,6,7]

    reconstructed:
    0123xx67

    """

    re_stream = stream.copy()

    # find the first valid pointer
    for first_pointer in pointer_array:
        if first_pointer != INVALID:
            break

    # find the last valid pointer
    for last_pointer in pointer_array[::-1]:
        if last_pointer != INVALID:
            break

    # number of samples
    n = last_pointer - first_pointer + 1

    # find the starttime
    for tr in re_stream.select(channel='_TT'):
        starttime0 = UTCDateTime(tr.data[first_pointer])

    # look through the whole stream
    for tr in re_stream:
        # change to the new channel name
        if tr.stats.channel == '_FR':
            tr.stats.channel = 'AFR'
        if tr.stats.channel == '_TT':
            tr.stats.channel = 'ATT'
        if tr.stats.channel == '_MZ':
            tr.stats.channel = 'MHZ'
        if tr.stats.channel == '_M1':
            tr.stats.channel = 'MH1'
        if tr.stats.channel == '_M2':
            tr.stats.channel = 'MH2'
        if tr.stats.channel == '_SZ':
            tr.stats.channel = 'SHZ'

        # set the data type
        if tr.stats.channel in ('AFR', 'ATT'):
            dtype = 'float64'
        else:
            dtype = 'int32'

        # make an array (with all invalid values)
        data_array = np.full(n, INVALID, dtype)

        # loop through the data from the trace
        for i, data in enumerate(tr.data[first_pointer:last_pointer+1]):
            pointer = first_pointer + i
            idx = pointer_array[pointer]
            if idx != INVALID:
                idx = idx - first_pointer
                # in a new array, set the data at the point given
                # by the correct index
                data_array[idx] = data

        # if any of the values are invalid, make a masked trace
        if INVALID in data_array:
            data_array = ma.masked_equal(data_array, INVALID)

        # update the trace data with the new data
        tr.data = data_array

        # update the starttime of the trace
        tr.stats.starttime = starttime0

    return re_stream, last_pointer

def _calc_match(starttime1, starttime0, framecount0, framecount1, obs_delta0):
    """
    Calculate the number of indices to move the sample on.

    framecount0 and starttime0 should be the last successful matches.

    """
    tim_diff = starttime1 - starttime0
    if tim_diff < 0 :
        raise ValueError('Time difference is negative.')
    sample_idx = tim_diff/obs_delta0

    # there are 90*4 = 360 samples per block
    # we can ignore whole blocks so just use the modulo operator
    part_block_samples = sample_idx % 360
    # round to nearest 0.25
    part_block_samples = (round(part_block_samples*4)/4)
    # change to framecounts (4 samples per frame number)
    part_block = part_block_samples/4
    # calculate the correct framecount for the time of the
    # new trace
    est_framecount1 = framecount0 + part_block
    # if the number is greater than 89.75, take 90 from it to get the correct number
    if est_framecount1 > 89.75:
        est_framecount1 = est_framecount1 - 90

    sample_idx = round(sample_idx)
    if sample_idx != 0 and framecount1 == est_framecount1:
    # make some more checks
        observed_delta = tim_diff / sample_idx
        # obs_delta_percent = round(100.*(DELTA-observed_delta)/DELTA,3)
        obs_delta_percent = 100.*(DELTA-observed_delta)/DELTA
        msg = 'Sampling interval:{0} Percent error from nominal interval:{1}%'.format(observed_delta, obs_delta_percent)
        logging.debug(msg)
        if abs(obs_delta_percent) > 1:
            msg = 'Reject automatically:{0} Percent error from nominal interval:{1}%'.format(observed_delta, obs_delta_percent)
            logging.debug(msg)
            sample_idx = None
    else:
        sample_idx = None

    return sample_idx
