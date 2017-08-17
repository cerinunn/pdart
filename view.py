#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
View the results.

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

from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, Trace, Stats, read
from urllib.request import Request, build_opener, HTTPError

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def plot_from_stream(
    stream=None,
    stations=['S11','S12','S14','S15','S16'],
    channels=['AFR'],
    merge_locations=False,
    plot_type='normal',
    time_interval=3):
    '''
    View the traces by passing in a stream.
    '''

    stream_all = stream

    for station in stations:

        # view by station
        stream = stream_all.select(station=station)
        stream = stream.sort(keys=['starttime'])

        for i, tr in enumerate(stream):
            # if i == 0:
            #     starttime0 = tr.stats.starttime
            if tr.stats.channel not in channels:
                stream.remove(tr)
            if merge_locations:
                tr.stats.location = ''

            if tr.stats.channel in ('ATT', '_TT'):
                # print(UTCDateTime()
                # this one shows a constant(ish) slope
                # tr.data = tr.data - starttime0.timestamp
                # this one should be close to zero
                tr.data = tr.data - tr.stats.starttime.timestamp - tr.times()

        # print("Sometimes it doesn't merge - I don't know why")
        # stream.merge()
        # stream.plot(size=(1200,600),method='full')

        if len(stream) > 0:
            if plot_type == 'dayplot':
                stream.plot(type ='dayplot',time_offset=0,show_y_UTC_label=False,
                  vertical_scaling_range=300)
            else:
                stream.plot(size=(1200,600),method='full')

def plot_from_file(
    stations=['S11','S12','S14','S15','S16'],
    channels=['AFR'],
    file_dir='.',
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T:21:00.000000Z'),
    merge_locations=False,
    plot_type='normal',
    time_interval=3,
    read_gzip=True):
    '''
    Read in files between the start and end times and plot.
    Calls plot_from_stream()
    '''


    if plot_type == 'dayplot':
        time_interval = timedelta(hours=24)
    else:
        time_interval = timedelta(hours=time_interval)

    for station in stations:

        file_dir_station =  os.path.join(file_dir, station)

        start = start_time
        while start < end_time:

            # use the same filename
            if plot_type == 'dayplot':
                stream_filename = '%s*_%s.MINISEED' % (start.strftime("%Y-%m-%dT"), station)
            elif plot_type == 'normal' and time_interval == timedelta(hours=24):
                stream_filename = '%s*_%s.MINISEED' % (start.strftime("%Y-%m-%dT"), station)
            else:
                stream_filename = '%s_%s.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
            if read_gzip:
                stream_filename = '%s.gz' % (stream_filename)
            stream_filename = os.path.join(file_dir_station, stream_filename)

            # read in the file
            try:
                stream = read(stream_filename)
                print(stream_filename)
            except FileNotFoundError:
                msg = 'view.py cannot find file: {}'.format(stream_filename)
                print(msg)
                logging.info(msg)
                # increment the time interval
                start += time_interval
                continue



            # select for this station
            stream = stream.select(station=station)
            stream = stream.sort(keys=['starttime'])

            plot_from_stream(
                stream=stream,
                stations=[station],
                channels=channels,
                merge_locations=merge_locations,
                plot_type=plot_type,
                time_interval=time_interval)

            # increment the time interval
            start += time_interval
