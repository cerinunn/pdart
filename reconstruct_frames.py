#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reconstruct the frames with the correct gaps.
It is sometimes necessary to manually exclude files.

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
from pdart.view import plot_from_stream
from pdart.splice import _samples_to_shift, splice_chains

import time


global DELTA
DELTA = 0.15094
SECONDS_PER_DAY = 3600.0 * 24.0
INVALID = -99999
ABSOLUTE_ADJUST_TIME = 2.
MIN_SAMPLE_LENGTH = 361




def call_reconstruct_frames(
    stations=['S11','S12','S14','S15','S16'],
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T:21:00.000000Z'),
    starttime0=None,
    framecount0=None,
    adjust0=None,
    obs_delta0=None,
    raw_dir='.',
    frame_dir='.',
    splice_dir='.',
    manual_remove_file=None,
    read_gzip=True,
    write_gzip=True,
    view_traces_to_remove=False,
    show_frames=True,
    save_frames_pdf=False,
    save_frames=True,
    show_splice=True,
    save_splice_pdf=False,
    save_splice=True,
    reset24=True,
    reset=False):
    '''
    Calls reconstruct_frames()
    '''

    log_filename = 'logs/reconstruct_frames.log'
    logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
    # logging.basicConfig(filename=log_filename, filemode='w', level=logging.DEBUG)


    for station in stations:

        # check that the overall directory exists
        if not os.path.exists(frame_dir):
            msg = ("The directory {} doesn't exist".format(frame_dir))
            raise IOError(msg)
        else:
            # make the subdirectory with the station name
            frame_dir_station =  os.path.join(frame_dir, station)
            if not os.path.exists(frame_dir_station):
                os.makedirs(frame_dir_station)

        # check that the overall directory exists
        if not os.path.exists(splice_dir):
            msg = ("The directory {} doesn't exist".format(splice_dir))
            raise IOError(msg)
        else:
            # make the subdirectory with the station name
            splice_dir_station =  os.path.join(splice_dir, station)
            if not os.path.exists(splice_dir_station):
                os.makedirs(splice_dir_station)

    if obs_delta0 is None:
        obs_delta0 = DELTA

    # build chains for each station
    for station in stations:

        raw_dir_station =  os.path.join(raw_dir, station)
        frame_dir_station =  os.path.join(frame_dir, station)
        splice_dir_station =  os.path.join(splice_dir, station)

        manual_remove_list = []
        if manual_remove_file is not None:
            with open(manual_remove_file, 'rt') as fh:
                for line in fh:
                    if line[0:19] != 'INFO:root:########:':
                        line = line.replace('INFO:root:######:','')
                        manual_remove_list.append(line)

        time_interval = timedelta(hours=3)
        start = start_time
        while start < end_time:

            if starttime0 is not None:
                msg = "starttime0=UTCDateTime('{}'); framecount0={}; adjust0={}; obs_delta0={}".format(str(starttime0), framecount0, adjust0, obs_delta0)
            else:
                msg = "starttime0=None; framecount0=None; adjust0=None"
            logging.info(msg)
            print(msg)
            msg = "start_time=UTCDateTime('{}')".format(str(start))
            logging.info(msg)
            print(msg)

            raw_filename = '%s_%s' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
            raw_filename = os.path.join(raw_dir_station, raw_filename)
            if read_gzip:
                raw_filename = '%s.MINISEED.gz' % (raw_filename)
            else:
                raw_filename = '%s.MINISEED' % (raw_filename)


            frame_filename = '%s_%s.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
            frame_filename = os.path.join(frame_dir_station, frame_filename)

            if write_gzip:
                frame_filename_gzip = '%s.gz' % (frame_filename)

            splice_filename = '%s_%s.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
            splice_filename = os.path.join(splice_dir_station, frame_filename)

            if write_gzip:
                splice_filename_gzip = '%s.gz' % (frame_filename)

            if save_frames_pdf:
                frame_outfile = '%s_%s.pdf' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
                frame_outfile = os.path.join(splice_dir_station, frame_outfile)
            else:
                frame_outfile = None

            if save_splice_pdf:
                splice_outfile = '%s_%s.pdf' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
                splice_outfile = os.path.join(splice_dir_station, splice_outfile)
            else:
                splice_outfile = None


            # read in the raw SEED file
            try:
                stream = read(raw_filename)
            except FileNotFoundError:
                msg = 'chains.py cannot find file: {}'.format(raw_filename)
                print(msg)
                logging.info(msg)
                # increment the time interval
                start += time_interval
                continue

            # select for just this station, (not necessary, but just in case)
            stream = stream.select(station=station)

            if len(stream) > 0:

                # build the chains (for this station)
                stream2 = reconstruct_frames(stream=stream, view_traces_to_remove=view_traces_to_remove, manual_remove_list=manual_remove_list)

                view_stream = stream + stream2

                plot_from_stream(
                    stream=view_stream,
                    stations=[station],
                    channels=['_FR','_TT','AFR','ATT'],
                    merge_locations=False,
                    time_interval=time_interval,
                    outfile=frame_outfile,
                    show=True)

                if len(stream2) > 0:
                    if save_frames:
                        stream2 = stream2.split()
                        stream2.write(frame_filename, 'MSEED')
                        # merge the stream again in case we need it later
                        stream2 = stream2.merge()

                        # this is slow
                        if write_gzip:
                            with open(frame_filename, 'rb') as f_in, gzip.open(frame_filename_gzip, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                            os.unlink(frame_filename)

                    if reset:
                        # use this when you want to check each file separately
                        starttime0=None
                        framecount0=None
                        adjust0=None
                        obs_delta0=DELTA

                    if reset24:
                        if start.hour == 0:
                            # use this when you want to check 24 hours
                            starttime0=None
                            framecount0=None
                            adjust0=None
                            obs_delta0=DELTA

                    # splice the chains for the period
                    return_stream, starttime0, framecount0, adjust0  = splice_chains(stream=stream2,
                      manual_remove_list=None,
                      starttime0=starttime0, framecount0=framecount0, adjust0=adjust0, obs_delta0=obs_delta0)

                    if len(return_stream) > 0:

                        plot_from_stream(
                            stream=return_stream,
                            stations=[station],
                            channels=['AFR','ATT','MHZ'],
                            merge_locations=True,
                            time_interval=time_interval,
                            outfile=splice_outfile,
                            show=show_splice)

                        if save_splice:
                            return_stream = return_stream.split()
                            return_stream.write(splice_filename, 'MSEED')

                            # this is slow
                            if write_gzip:
                                with open(splice_filename, 'rb') as f_in, gzip.open(splice_filename_gzip, 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                                os.unlink(splice_filename)

            # increment the time interval
            start += time_interval

def reconstruct_frames(stream, view_traces_to_remove=False, manual_remove_list=None):
    '''
    Reconstruct frames by checking the framecount and inserting gaps.
    '''
    log_filename = 'logs/reconstruct_frames.log'
    logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
    # logging.basicConfig(filename=log_filename, filemode='w', level=logging.DEBUG)

    # quick check to make sure only one station
    station = stream[0].stats.station

    if len(stream) != len(stream.select(station=station)):
        raise ValueError("More than one station in the stream")

    # TODO only need 3 if running SPZ
    # original_data = np.full((n,3), INVALID, 'int32')

    # get rid of any short traces that overlap with other
    # streams

    stream = discard_short_traces(stream)

    # if we need anything to remove, remove it here
    if manual_remove_list is not None:
        stream = manual_trim(stream,manual_remove_list)

    # begin by selecting the raw trace, and sorting by starttime
    FR_stream = stream.select(channel='_FR')
    FR_stream = FR_stream.sort(keys=['starttime'])

    return_stream = Stream()

    # for each of the raw streams
    for fs in FR_stream:

        if view_traces_to_remove:
            msg = '######:{} {} {} {} {} {}'.format(fs.stats.network,
              fs.stats.station, fs.stats.location, fs.stats.starttime,
              fs.stats.endtime, fs.stats.npts)
            logging.info(msg)

        # find the matching streams with the same start and end time
        original = stream_select(stream,network=fs.stats.network, station=fs.stats.station,
          location=fs.stats.location,starttime=fs.stats.starttime,endtime=fs.stats.endtime)

        # get the stream for the timing trace
        TT_stream = original.select(channel='_TT')

        # TODO something fancy with the invalid frames - especially at the start.

        # make a pointer_array
        n = fs.stats.npts
        pointer_array = np.full(n, INVALID, 'int32')

        # loop through data
        for pointer, framecount1 in enumerate(fs.data):

            if pointer == 0:
                timestamp0 = TT_stream[0].data[0]
                framecount0 = framecount1
                index_pointer = pointer
                # update the pointers in the pointer array
                pointer_array[pointer] = index_pointer
            else:
                timestamp1 = TT_stream[0].data[pointer]

                # use for debugging
                # if i < 6:
                #     msg = ('i = {} {} {}'.format(i, start_pointer, framecount1))
                #     logging.debug(msg)

            # if the framecount is out of range, start again
                if framecount1 < 0 or framecount1 > 89.75:
                    # we just ignore it
                    continue
                # if the frame range is valid
                else:
                    # check for the correct sample index from the framecount
                    # and timestamp
                    sample_diff = _calc_match(timestamp1, timestamp0, framecount0, framecount1, obs_delta0=DELTA)

                    if sample_diff is not None:
                        index_pointer += sample_diff
                        # if ixx > 1580:
                        #     print(index_pointer,sample_diff, ixx )

                        # the current one is valid, so update the
                        # timestamp and framecount
                        timestamp0 = timestamp1
                        framecount0 = framecount1

                        # update the pointers in the pointer array
                        pointer_array[pointer] = index_pointer

        # reconstruct the stream
        re_stream, last_pointer = _reconstruct_streams(original, pointer_array)

        return_stream += re_stream

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
    # don't forget that the pointer values can be larger then the original
    # data length!
    for last_pointer in pointer_array[::-1]:
        if last_pointer != INVALID:
            break

    # number of samples
    n = last_pointer - first_pointer + 1

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
        for pointer, data in enumerate(tr.data):
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

        # calculate statistics for the NEW trace
        # do this after reforming the traces
        if tr.stats.channel == 'ATT':
            starttime0 = UTCDateTime(tr.data[0])
            end = UTCDateTime(tr.data[-1])
            obs_delta0 = (end - starttime0)/(tr.stats.npts - 1)
            msg = 'starttime0 {} end {} {} obs_delta0 {}'.format(starttime0, end, tr.data[-1], obs_delta0)
            logging.debug(msg)

    # loop again, now we know the starttime
    for tr in re_stream:
        # update the starttime
        tr.stats.starttime = starttime0

    return re_stream, last_pointer

def _calc_match(timestamp1, timestamp0, framecount0, framecount1, obs_delta0):
    """
    Calculate the number of indices to move the sample on.

    The match is based on the FRAMECOUNT, with checks based on the
    time difference.

    framecount0 and timestamp0 should be the last successful matches.

    """
    tim_diff = timestamp1 - timestamp0
    if tim_diff < 0 :
        raise ValueError('Time difference is negative.')
    sample_idx = tim_diff/obs_delta0

    # there are 90*4 = 360 samples per block
    # we can ignore whole blocks so just use the modulo operator
    part_block_samples = sample_idx % 360
    # print(part_block_samples)
    # round to nearest 0.25
    # part_block_samples = (round(part_block_samples*4)/4)
    # round the sampling
    part_block_samples = round(part_block_samples)
    # change to framecounts (4 samples per frame number)
    part_block = part_block_samples/4
    # print('part block samples ', sample_idx, part_block_samples,part_block)
    # calculate the correct framecount for the time of the
    # new trace
    est_framecount1 = framecount0 + part_block
    # if the number is greater than 89.75, take 90 from it to get the correct number
    if est_framecount1 > 89.75:
        est_framecount1 = est_framecount1 - 90

    sample_idx = round(sample_idx)
    # if sample_idx != 0 and framecount1 == est_framecount1:
    # # make some more checks
    #     observed_delta = tim_diff / sample_idx
    #     # obs_delta_percent = round(100.*(DELTA-observed_delta)/DELTA,3)
    #     obs_delta_percent = 100.*(DELTA-observed_delta)/DELTA
    #     msg = 'Sampling interval:{0} Percent error from nominal interval:{1}%'.format(observed_delta, obs_delta_percent)
    #     logging.debug(msg)
    #     # up to a 10% error (we make a new chain if data is not consecutive -
    #     # so this only applies to chains that are OK.
    #     if abs(obs_delta_percent) > 10:
    #         msg = 'Reject automatically:{0} Percent error from nominal interval:{1}%'.format(observed_delta, obs_delta_percent)
    #         logging.debug(msg)
    #         sample_idx = None
    # else:
    #     sample_idx = None
    #

    if sample_idx < 1 or framecount1 != est_framecount1:
        msg = 'Framecount does not match {} {} {} {}'.format(framecount1,est_framecount1, sample_idx, UTCDateTime(timestamp1))
        sample_idx = None
        logging.debug(msg)

    return sample_idx

def discard_short_traces(stream):
    """
    Discard short traces which overlap with longer traces.
    Short traces are often quite poor quality - but sometimes occur when
    a longer trace exists. If so, they can be safely discarded.
    """

    # copy the original stream
    return_stream = stream.copy()

    # delete any traces from location 00, but give a warning
    loc = return_stream.select(location='00')
    for tr in loc:
        if tr.stats.channel == '_TT':
            msg = 'Removing trace because it is from location 00: {}'.format(str(tr))
            logging.warning(msg)
            print(msg)
        return_stream.remove(tr)


    # sort a stream (from the timing channel) with the number of samples
    sorted_stream = return_stream.select(channel='_TT').sort(keys=['npts'])

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
                        # print(msg)
                        logging.debug(msg)
                        stream_short = stream_select(stream,network=tr.stats.network, station=tr.stats.station,
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

# def manual_remove(stream, manual_remove_list):
#     for line in manual_remove_list:
#         (network, station, location, starttime, endtime, npts) = line.split()
#         starttime = UTCDateTime(starttime)
#         endtime = UTCDateTime(endtime)
#
#         # print (network, station, location, starttime, endtime, npts)
#
#         match = stream_select(stream,network=network,station=station,location=location,starttime=starttime,endtime=endtime)
#         for mat in match:
#             stream.remove(mat)
#
#     return stream

def manual_trim(stream, manual_remove_list):
    for line in manual_remove_list:
        (network, station, location, starttime, endtime, npts) = line.split()
        begin_trim = UTCDateTime(starttime)
        end_trim = UTCDateTime(endtime)

        for tr in stream:
            starttime = tr.stats.starttime
            endtime = tr.stats.endtime
            tr1 = None
            tr2 = None
            if begin_trim >= starttime and begin_trim <= endtime:
                tr1 = tr.slice(endtime=begin_trim-DELTA)
                if tr1.stats.npts > 1:
                    stream.append(tr1)
            if end_trim >= starttime and end_trim <= endtime:
                tr2 = tr.slice(starttime=end_trim+DELTA)
                if tr2.stats.npts > 1:
                    stream.append(tr2)
            if (tr1 is not None) or (tr2 is not None):
                if tr1.stats.channel == '_TT':
                    print('Trace in manual remove file - removing. ', tr)
                stream.remove(tr)
    return stream
