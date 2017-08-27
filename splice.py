#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Splices the records (chains) together.

:copyright:
    The PDART Development Team & Ceri Nunn
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
from urllib.request import Request, build_opener, HTTPError

from pdart.chains import _samples_to_shift
from pdart.util import stream_select, trace_eq

global DELTA
DELTA = 0.15094
SECONDS_PER_DAY = 3600.0 * 24.0
INVALID = -99999
ABSOLUTE_ADJUST_TIME = 2.

def call_splice_chains(
    stations=['S11','S12','S14','S15','S16'],
    chain_dir='.',
    splice_dir='.',
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T:21:00.000000Z'),
    read_gzip=True,
    write_gzip=True,
    manual_remove_file=None,
    reset=True):

    '''
    Calls splice_chains()
    Remember that the absolute timing depends on the previous stream,
    so run call_splice_chains in date order.
    '''

    log_filename = 'logs/splice.log'
    logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
    # logging.basicConfig(filename=log_filename, filemode='w', level=logging.DEBUG)

    for station in stations:

        # check that the overall directory exists
        if not os.path.exists(splice_dir):
            msg = ("The directory {} doesn't exist".format(splice_dir))
            raise IOError(msg)
        else:
            # make the subdirectory with the station name
            splice_dir_station =  os.path.join(splice_dir, station)
            if not os.path.exists(splice_dir_station):
                os.makedirs(splice_dir_station)

        manual_remove_list = []
        if manual_remove_file is not None:
            with open(manual_remove_file, 'rt') as fh:
                for line in fh:
                    if line[0:19] != 'INFO:root:########:':
                        line = line.replace('INFO:root:######:','')
                        manual_remove_list.append(line)

    # run for each station
    for station in stations:

        starttime0=None
        framecount0=None
        adjust0=None

        chain_dir_station =  os.path.join(chain_dir, station)
        splice_dir_station =  os.path.join(splice_dir, station)

        time_interval = timedelta(hours=3)
        start = start_time
        while start < end_time:

            chain_filename = '%s_%s*.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
            if read_gzip:
                chain_filename = '%s*.gz' % (chain_filename)

            splice_filename = '%s_%s.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
            splice_dir_filename = os.path.join(splice_dir_station,
              splice_filename)
            if write_gzip:
                gzip_splice_filename = '%s.gz' % (splice_filename)
                gzip_splice_dir_filename = os.path.join(splice_dir_station,
                  gzip_splice_filename)


            chain_dir_filename = os.path.join(chain_dir_station, chain_filename)

            msg = '########: {}'.format(chain_dir_filename)
            logging.info(msg)

            # check that the file pattern exists
            if not glob.glob(chain_dir_filename):
                msg = 'splice.py cannot find file: {}'.format(chain_dir_filename)
                print(msg)
                logging.info(msg)
                # increment the time interval
                start += time_interval
                continue

            # read in the file
            stream = read(chain_dir_filename)

            # select for just this station, (not necessary, but just in case)
            stream = stream.select(station=station)

            if reset:
                # use this when you want to check each file separately
                starttime0=None
                framecount0=None
                adjust0=None

            if len(stream) > 0:

                # splice the chains for the period
                return_stream, starttime0, framecount0, adjust0  = splice_chains(stream=stream,
                  manual_remove_list=manual_remove_list,
                  starttime0=starttime0, framecount0=framecount0, adjust0=adjust0)

                msg = "{} starttime0=UTCDateTime('{}'); framecount0={}; adjust0={}".format(chain_dir_filename, str(starttime0), framecount0, adjust0)
                logging.info(msg)

                # it returns starttime0, framecount0, adjust0, so that we
                # can run correctly for the next period

                if len(return_stream) > 0:
                    # split the streams in order to save them
                    return_stream = return_stream.split()
                    # save the streams
                    return_stream.write(splice_dir_filename, 'MSEED')
                    if write_gzip:
                        with open(splice_dir_filename, 'rb') as f_in, gzip.open(gzip_splice_dir_filename, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                        os.unlink(splice_dir_filename)
                else:
                    msg = 'No spliced records for {}.'.format(chain_dir_filename)
                    logging.info(msg)

            else:
                msg = 'No records found for {}.'.format(chain_dir_filename)
                logging.info(msg)

            # increment the time interval
            start += time_interval

def splice_chains(stream, manual_remove_list=None, starttime0=None, framecount0=None, adjust0=None):
    '''
    Splice the records (chains) together.
    Remember that the absolute timing depends on the previous stream,
    so run call_splice_chains in date order.
    '''
    log_filename = 'logs/splice.log'
    logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
    # logging.basicConfig(filename=log_filename, filemode='w', level=logging.DEBUG)

    if adjust0 is None:
        if starttime0 is not None or framecount0 is not None:
            msg = 'If adjust0 is not set, starttime0 and framecount0 must also be set.'
            logging.warning(msg)
            exit()

    # quick check to make sure only one station
    station = stream[0].stats.station
    if len(stream) != len(stream.select(station=station)):
        raise ValueError("More than one station in the stream")

    # begin by selecting the frame traces, and sorting by starttime
    frm_stream = stream.select(channel='AFR')
    frm_stream = frm_stream.sort(keys=['starttime'])

    # if we need anything to remove, remove it here
    if manual_remove_list is not None:
        frm_stream = manual_remove(frm_stream,manual_remove_list)

    return_stream = Stream()

    # for i, fs in enumerate(frm_stream):
    #     fs.plot()
    #     if i > 4:
    #         exit()
    # exit()

    # for ix, fs in enumerate(frm_stream):
    #     if ix > 4:
    #         exit()
    for fs in frm_stream:

        if type(fs.data) == ma.MaskedArray:
            percent_invalid = ma.count_masked(fs.data)/len(fs.data)*100.
        else:
            percent_invalid = 0

        msg = ('######: {} {} {} {} {} {} {}%'.format(fs.stats.network, fs.stats.station, fs.stats.location, fs.stats.starttime, fs.stats.endtime, fs.stats.npts, round(percent_invalid,2)))
        logging.info(msg)

        # if starttime0, framecount0 and adjust0
        # have not already been set, we set them from the earliest
        # trace in the stream
        if adjust0 is None:
            starttime0 = fs.stats.starttime
            framecount0 = fs.data[0]
            adjust0 = 0
            # find the matching traces
            match = stream_select(stream,network=fs.stats.network, station=fs.stats.station,
              location=fs.stats.location,starttime=fs.stats.starttime,endtime=fs.stats.endtime)
            return_stream += match
            continue

        # find the starttime and framecount of the current trace
        starttime1 = fs.stats.starttime
        endtime1 = fs.stats.endtime
        framecount1 = fs.data[0]

        # adjust the starttime
        adjust_starttime0 = starttime0 - adjust0

        # estimate the sample number of the current trace, assuming
        # it continues from the successful trace
        sample_idx = _calc_match_samp_frame(starttime1, adjust_starttime0, framecount0, framecount1, obs_delta0=DELTA)

        valid_chain = False

        # sample_idx is None when it is invalid
        if sample_idx is not None:
            # estimate the new starttime for the current trace
            est_starttime1 = starttime0 + (sample_idx * DELTA)
            # check that the adjust time is not getting too large
            adjust1 = est_starttime1 - starttime1
            msg = '{} {} {} {} {}'.format(est_starttime1,starttime1, adjust1, adjust0, abs(adjust1 - adjust0))
            logging.debug(msg)
            if abs(adjust1 - adjust0) < ABSOLUTE_ADJUST_TIME:
                valid_chain = True
        # else:
        #     logging.debug('None')

        # update the startimes for the traces which match the other details
        st_update = stream_select(stream,network=fs.stats.network, station=fs.stats.station,
                  location=fs.stats.location,starttime=fs.stats.starttime,endtime=fs.stats.endtime)

        # loop through the matching traces
        for i, tr in enumerate(st_update):
            # traces are either updated or removed from the original stream
            if valid_chain:
                if i==0:
                    # record the change in the log for the first trace only
                    msg = 'adjust_time:{}, for station: {}, location: {}, starttime: {}, endtime: {}'.format(adjust1,
                      tr.stats.station, tr.stats.location, starttime1, endtime1 )
                    logging.debug(msg)
                    # update starttime0, framecount0, adjust0 with details from this trace
                    starttime0 = est_starttime1
                    framecount0 = framecount1
                    adjust0 = adjust1
                # adjust trace starttime
                tr.stats.starttime = est_starttime1
            else:
                # throw the trace away
                if i==0:
                    msg = 'Remvoing this trace: for station: {}, location: {}, , starttime: {}, endtime: {}'.format(
                      tr.stats.station, tr.stats.location,  tr.stats.starttime, tr.stats.endtime )
                    logging.debug(msg)
                # remove the stream from the trace we passed in to the method
                st_update.remove(tr)
        return_stream += st_update


    if len(return_stream) > 0:
        # calculate the length and write to the log file
        length_stream = return_stream.select(channel='ATT')
        length_stream = length_stream.sort(keys=['starttime'])
        length_stream2 = length_stream.copy()
        length_stream2 = length_stream2.sort(keys=['endtime'], reverse=True)
        elapsed_time = length_stream2[0].stats.endtime - length_stream[0].stats.starttime
        elapsed_timestamps = length_stream2[0].data[-1] - length_stream[0].data[0]
        obs_delta = elapsed_timestamps / ((elapsed_time)/DELTA)
        msg = ('elapsed_time: {} elapsed_timestamps: {} obs_delta: {}'.format(round(elapsed_time,3), round(elapsed_timestamps,3),obs_delta))
        if elapsed_timestamps > 10801. or elapsed_timestamps < 10799:
            logging.warning(msg)
        else:
            logging.info(msg)


    return return_stream, starttime0, framecount0, adjust0

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

    idx_diff = 0

    if framecount1 != est_framecount1:
        # if the trace needs shifting sample numbers, do it here:
        idx_diff = _samples_to_shift(framecount1, est_framecount1)

        if abs(idx_diff) < 12:
            sample_idx = sample_idx + idx_diff
        else:
            sample_idx = None
            msg = 'idx_diff == {} - reject automatically'.format(idx_diff)
            logging.info(msg)

    msg = 'sample_idx {} idx_diff {} framecount1 {} est_framecount1 {}'.format(sample_idx, idx_diff, framecount1, est_framecount1)
    logging.debug(msg)
    return sample_idx


def manual_remove(stream, manual_remove_list):
    for line in manual_remove_list:
        (network, station, location, starttime, endtime, npts, percent_invalid) = line.split()
        starttime = UTCDateTime(starttime)
        endtime = UTCDateTime(endtime)

        match = stream_select(stream,network=network,station=station,location=location,starttime=starttime,endtime=endtime)
        for mat in match:
            stream.remove(mat)

    return stream
