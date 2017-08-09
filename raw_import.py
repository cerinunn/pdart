#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run the initial import from the csv files to the raw SEED files.
The raw files have not had the frames reconstructed yet.

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

global DELTA
DELTA = 0.15094
INVALID = -99999


def call_raw_import(
    stations = ['S11','S12','S14','S15','S16'],
    csv_dir='.',
    raw_dir='.',
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T:21:00.000000Z'),
    write_gzip=True):
    '''
    Makes 'raw' stream files from gzipped csv files.
    Calls raw_import()
    '''

    for station in stations:
        # check that the overall directory exists
        if not os.path.exists(raw_dir):
            msg = ("The directory {} doesn't exist".format(raw_dir))
            raise IOError(msg)
        else:
            # make the subdirectory with the station name
            raw_dir_station =  os.path.join(raw_dir, station)
            if not os.path.exists(raw_dir_station):
                os.makedirs(raw_dir_station)

    time_interval = timedelta(hours=3)
    start = start_time

    while start < end_time:

        # find csv gzip importfile  with the start time
        csv_file_gzip = '%s.csv.gz' % start.strftime("%Y-%m-%dT%H:%M:%S")
        csv_file_gzip = os.path.join(csv_dir, csv_file_gzip)

        # import the stream from the gzipped file
        stream = raw_import(csv_file_gzip)

        if len(stream) > 0:

            stations = []
            for tr in stream:
                stations.append(tr.stats.station)
            stations = list(set(stations))

            for station in stations:
                station_stream = stream.select(station=station)
                # station_stream = station_stream.split()
                stream_file = '%s_%s.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
                gzip_file = '%s.gz' % (stream_file)
                raw_dir_station = os.path.join(raw_dir, station)
                stream_file = os.path.join(raw_dir_station, stream_file)
                gzip_file = os.path.join(raw_dir_station, gzip_file)
                stream.write(stream_file, 'MSEED')
                if write_gzip:
                    with open(stream_file, 'rb') as f_in, gzip.open(gzip_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    os.unlink(stream_file)

        # increment the time interval
        start += time_interval

    print('Done')

def raw_import(gzip_filename):
    """
    Makes a 'raw' stream file from the gzipped csv file.
    The csv file has been downloaded from the JAXA website.
    The method makes a raw stream which does not yet have the frames
    reconstructed.

    :type gzip_filename: str
    :param gzip_filename: gzipped filename of the CSV file to be read.
    :rtype: :class:`~obspy.core.stream.Stream`
    :return: A ObsPy Stream object.

    """
    # read the gzipped csv file
    with gzip.open(gzip_filename, 'rt') as fh:
        # read file
        buf = []
        header = next(fh).split(',')

        # read the header
        # it should contain either 1 channel or 3
        if len(header) == 8:
            # the RESP files use either 'MH1', 'MH2', 'MHZ'
            # the JAXA files use 'LPX', 'LPY', 'LPZ'
            # X should point north, Y east, but this is not always the case
            # so we rename LPX to MH1, and LPY to MH2
            channels = ['MH1', 'MH2', 'MHZ']
            raw_channels = ['_M1', '_M2', '_MZ']
            for line in fh:
                temp = line.split(',')

                try:
                    temp[4] = UTCDateTime(temp[4])
                except ValueError as e:
                    # this is a specific error which is found in the csv file
                    if temp[4] == '1975-49-11 19:13:04.232000':
                        temp[4] = UTCDateTime('1975-09-11 19:13:04.232000')
                    else:
                        raise

                try:
                    temp[0] = int(temp[0])
                except ValueError as e:
                    # this is a specific error which is found in the csv file
                    if temp[4] == UTCDateTime('1975-09-15 12:53:36.849000') and temp[0] == '<3':
                        temp[0] = 83
                    else:
                        raise

                buf.append((temp[1], temp[2], temp[4],
                  int(temp[0]), int(temp[3]), int(temp[5]), int(temp[6]),
                  int(temp[7])))

        elif len(header) == 6:
            channels = ['SPZ']
            raw_channels = ['_SZ']
            for line in fh:
                # check the manual list of points which have been removed
                if line in remove_manually:
                    continue

                temp = line.split(',')
                # the original order:
                # frame_count, ap_station, ground_station, nc, time, spz
                # make a tuple (in a new order so that it can be sorted):
                # ap_station, ground_station, time, frame_count, nc, spz
                buf.append((temp[1], temp[2],
                  UTCDateTime(temp[4]),
                  int(temp[0]), int(temp[3]), int(temp[5])))

    # sort by ap_station, ground_station and time (and also everything else,
    # but that won't matter)
    buf.sort()

    stream = Stream()
    data_x = []
    data_y = []
    data_z = []
    data_sz = []
    abs_times = []
    frame_count_ncs = []
    corr_frame_count_ncs = []

    stats = Stats()
    stats.delta = DELTA
    network = 'XA'
    last_id = None

    for data in buf:

        # read in the data from the buffer
        station = data[0].rjust(3,'S')
        ground_station = data[1].rjust(2,'0')
        time = data[2]

        frame_count = data[3]
        nc = data[4]
        # create a combination of frame count and nc - from 0.0 to 89.75
        frame_count_nc = float(frame_count) + (float(nc) - 1.) * 0.25

        id = "{0:s}.{1:s}.{2:s}.{3:s}".format(
          network, station, ground_station, channels[0])

        # check whether we are adding to an existing one, or creating a new one
        if (last_id is None or last_id != id):
                # before creating the new one, add previous trace(s) to the stream
                if len(abs_times) > 0:
                    _make_traces(stream=stream, stats=stats, header=header,channels=raw_channels, data_x=data_x,data_y=data_y, data_z=data_z, data_sz=data_sz,abs_times=abs_times,frame_count_ncs=frame_count_ncs)

                data_x = []
                data_y = []
                data_z = []
                data_sz = []
                abs_times = []
                frame_count_ncs = []

                stats = Stats()
                stats.delta = DELTA
                stats.starttime = time
                stats.network = network
                stats.station = station
                stats.location = ground_station

        # add the data) from any line
        if len(header) == 8:
            data_x.append(data[5])
            data_y.append(data[6])
            data_z.append(data[7])
        else:
            data_sz.append(data[5])
        abs_times.append(time.timestamp)
        frame_count_ncs.append(frame_count_nc)

        last_id = id

    # add the last one
    if len(abs_times) > 0:
        _make_traces(stream=stream, stats=stats, header=header,channels=raw_channels, data_x=data_x,data_y=data_y, data_z=data_z, data_sz=data_sz,abs_times=abs_times,frame_count_ncs=frame_count_ncs)

    return stream


def _make_traces(stream=None, stats=None, header=None, channels=None,
        data_x=None, data_y=None, data_z=None, data_sz=None,
        abs_times=None,frame_count_ncs=None):
    '''
    Make the traces from the lists imported from the csv file.
    '''

    if stream is None:
        stream = Stream()

    if len(abs_times) > 1:
        if len(header) == 8:
            stats.channel = channels[0]
            _append_stream(stream, data_x, stats)
            stats_y = stats.copy()
            stats_y.channel = channels[1]
            _append_stream(stream, data_y, stats_y)
            stats_z = stats.copy()
            stats_z.channel = channels[2]
            _append_stream(stream, data_z, stats_z)
        else:
            stats.channel = channels[0]
            _append_stream(stream, data_sz, stats)

        stats_times = stats.copy()
        stats_times.channel = '_TT'
        _append_stream(stream, abs_times, stats_times)

        stats_frames= stats.copy()
        stats_frames.channel = '_FR'
        _append_stream(stream, frame_count_ncs, stats_frames)


def _append_stream(stream, data_list, stats):
    '''
    Add the trace to the stream (whilst checking the data type)
    '''
    if type(data_list[0]) == float:
        parsed_data = np.array(data_list,'float64')
    else:
        parsed_data = np.array(data_list,'int32')

    if INVALID in parsed_data:
        parsed_data = ma.masked_equal(parsed_data, INVALID)

    stats.npts = len(parsed_data)

    # append the trace
    stream.append(Trace(data=parsed_data, header=stats))

    return stream
