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
from pdart.util import stream_select

global DELTA
DELTA = 0.15094
SECONDS_PER_DAY = 3600.0 * 24.0
INVALID = -99999
ABSOLUTE_ADJUST_TIME = 2.
MIN_SAMPLE_LENGTH = 361

def rename(
    stations=['S11','S12','S14','S15','S16'],
    starttime0=None,
    framecount0=None,
    adjust0=None,
    obs_delta0=DELTA,
    timestamp0=None,
    framecount_adjust=None,
    raw_dir='.',
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T:21:00.000000Z'),
    read_gzip=True,
    write_gzip=True):
    '''
    Calls build_chains()
    '''

    # log_filename = 'logs/rename.log'
    # logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)

    # build chains for each station
    for station in stations:

        raw_dir_station =  os.path.join(raw_dir, station)

        time_interval = timedelta(hours=3)
        start = start_time
        while start < end_time:

            # work out the filenames
            stream_filename = '%s_%s.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
            gzip_filename = '%s_%s.MINISEED.gz' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)

            stream_dir_filename = os.path.join(raw_dir_station, stream_filename)
            gzip_stream_dir_filename = os.path.join(raw_dir_station, gzip_filename)


            # read in the raw SEED file
            try:
                print(gzip_stream_dir_filename)
                stream = read(gzip_stream_dir_filename)
            except FileNotFoundError:
                msg = 'rename.py cannot find file: {}'.format(gzip_stream_dir_filename)
                print(msg)
                logging.info(msg)
                # increment the time interval
                start += time_interval
                continue

            if len(stream) > 0:
                station_stream = stream.select(station=station)
                for tr in station_stream:
                    if tr.stats.channel == '_fr':
                        tr.stats.channel = '_FR'
                    elif tr.stats.channel == '_ti':
                        tr.stats.channel = '_TT'
                    else:
                        continue

                station_stream.write(stream_dir_filename, 'MSEED')
                if write_gzip:
                    with open(stream_dir_filename, 'rb') as f_in, gzip.open(gzip_stream_dir_filename, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    os.unlink(stream_dir_filename)


            # increment the time interval
            start += time_interval

if __name__ == "__main__":
    start_time=UTCDateTime('1977-09-30T18:00:00.000000Z')
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z')

    start_time=UTCDateTime('1976-01-02T00:00:00.000000Z')



    raw_dir='/Users/nunn/lunar_data/JAXA_RAW'

    rename(
        stations=['S12'],
        start_time=start_time,
        raw_dir=raw_dir)
