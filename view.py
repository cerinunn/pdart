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
from matplotlib.dates import AutoDateLocator, MonthLocator
from dateutil.rrule import MONTHLY, YEARLY


from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, Trace, Stats, read
from obspy.imaging.util import ObsPyAutoDateFormatter
from urllib.request import Request, build_opener, HTTPError


from pdart.util import stream_select
from pdart.chains import discard_short_traces

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

global DELTA
DELTA = 0.15094

def plot_from_stream(
    stream=None,
    stations=['S11','S12','S14','S15','S16'],
    channels=['AFR'],
    merge_locations=False,
    plot_type='normal',
    time_interval=3,
    outfile=None,
    show=True,
    view_traces_to_remove=False):
    '''
    View the traces by passing in a stream.
    '''

    stream_all = stream

    for station in stations:

        # view by station (make a copy because we modify the trace)
        stream = stream_all.select(station=station).copy()
        stream = stream.sort(keys=['starttime'])

        if view_traces_to_remove:
            stream = discard_short_traces(stream)

        for i, tr in enumerate(stream):
            # if i == 0:
            #     starttime0 = tr.stats.starttime
            if view_traces_to_remove:
                if tr.stats.channel == '_FR':
                    msg = '{} {} {} {} {} {}'.format(tr.stats.network,
                      tr.stats.station, tr.stats.location, tr.stats.starttime,
                      tr.stats.endtime, tr.stats.npts)
                    logging.info(msg)

            if tr.stats.channel not in channels:
                stream.remove(tr)

            if merge_locations:
                tr.stats.location = ''

            # if tr.stats.channel in ('ATT'):
            #     obs_delta = (tr.data[-1] - tr.data[0])/(tr.stats.npts - 1)
                # print('obs delta ', obs_delta, tr.id, tr.stats.starttime)

            if tr.stats.channel in ('ATT', '_TT'):
                # print(UTCDateTime()
                # this one shows a constant(ish) slope
                # tr.data = tr.data - starttime0.timestamp
                # this one shows the drift between times() and the timestamp
                # (should be only a few 10ths of second in 24 hours)
                tr.data = tr.data - tr.stats.starttime.timestamp - tr.times()


        # print("Sometimes it doesn't merge - I don't know why")
        # stream.merge()
        # stream.plot(size=(1200,600),method='full')

        if len(stream) > 0:
            if plot_type == 'dayplot':
                if outfile is not None:
                    stream.plot(type ='dayplot',time_offset=0,show_y_UTC_label=False,
                      vertical_scaling_range=300,outfile=outfile)
                if show:
                    stream.plot(type ='dayplot',time_offset=0,show_y_UTC_label=False,
                      vertical_scaling_range=300)
            else:
                # there is a problem which means choosing both options doesn't work
                if outfile is not None:
                    stream.plot(size=(1200,600),method='full', equal_scale=False, outfile=outfile,)
                if show:
                    stream.plot(size=(1200,600),method='full',equal_scale=False)



def plot_from_file(
    stations=['S11','S12','S14','S15','S16'],
    channels=['AFR'],
    file_dir='.',
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T:21:00.000000Z'),
    merge_locations=False,
    plot_type='normal',
    time_interval=3,
    read_gzip=True,
    save_pdf=False,
    view_traces_to_remove=False):
    '''
    Read in files between the start and end times and plot.
    Calls plot_from_stream()
    '''

    if view_traces_to_remove:
        log_filename = 'logs/view_traces_to_remove.log'
        # logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
        logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)

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
                stream_filename = '%s*_%s*.MINISEED' % (start.strftime("%Y-%m-%dT"), station)
            elif plot_type == 'normal' and time_interval == timedelta(hours=24):
                stream_filename = '%s*_%s*.MINISEED' % (start.strftime("%Y-%m-%dT"), station)
            else:
                stream_filename = '%s_%s*.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)

            if read_gzip:
                stream_filename = '%s.gz' % (stream_filename)
            stream_filename = os.path.join(file_dir_station, stream_filename)

            if save_pdf:
                outfile = '%s_%s.pdf' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
                outfile = os.path.join(file_dir_station, outfile)
            else:
                outfile = None

            # check that the file pattern exists
            if not glob.glob(stream_filename):
                msg = 'view.py cannot find file: {}'.format(stream_filename)
                print(msg)
                logging.info(msg)
                # increment the time interval
                start += time_interval
                continue

            stream = read(stream_filename)
            #
            # print('Temporary!!!')
            # stream = stream_select(stream,
            #            starttime=UTCDateTime('1971-08-01T16:10:38.923000Z'))

            # select for this station
            stream = stream.select(station=station)
            stream = stream.sort(keys=['starttime'])


            plot_from_stream(
                stream=stream,
                stations=[station],
                channels=channels,
                merge_locations=merge_locations,
                plot_type=plot_type,
                time_interval=time_interval,
                outfile=outfile,
                view_traces_to_remove=view_traces_to_remove)

            # increment the time interval
            start += time_interval






##########

def plot_from_stream_all_stations(
    stream=None,
    channels=['AFR'],
    merge_locations=False,
    time_diffs=True,
    save_pdf=False,
    outfile=None):
    '''
    View the traces by passing in a stream.
    '''

    for i, tr in enumerate(stream):

        if merge_locations:
            tr.stats.location = ''

        if tr.stats.channel in ('ATT', '_TT'):
            # delta_times = tr.data.copy()
            # for i in range(len(delta_times)):
            #     delta_times[i] = i * obs_delta

            # calculate the difference between consecutive values in the
            # time trace.
            # this is 1 sample shorter than the others
            tr.data = np.diff(tr.data)

            # now calculate the average for the frame
            n_frame = 360
            # if you want to view it without averages, set n_frame=1
            # n_frame = 1
            N = tr.stats.npts//n_frame
            if N > 0:
                means = np.zeros(N)
                for i in range(N):
                    means[i] = tr.data[i*n_frame:i*n_frame+n_frame].mean()

                tr.data = means
                tr.stats.delta = DELTA * n_frame
            else:
                stream.remove(tr)

    # print("Sometimes it doesn't merge - I don't know why")
    # stream.merge()
    # stream.plot(size=(1200,600),method='full')



    if len(stream) > 0:
        stream.write('/Users/nunn/lunar_data/temp2/all.MINISEED', 'MSEED')

        if time_diffs:
            fig = stream.plot(handle=True,size=(1200,600),show=False,
              starttime=UTCDateTime('1970-01-15T03:00:00.000000Z'),
              endtime=UTCDateTime('1977-09-30T:21:00.000000Z'))
            # clear the suptitle which has already been set
            fig.suptitle('')
            axes = fig.get_axes()

            # plt.ylabel('s')
            for i, ax in enumerate(axes):
                # view unaveraged time diffs
                # ax.set_ylim(0.14,0.16)
                # per frame time diffs
                if i == 0:
                    ax.set_title('Apparent Sampling Interval')
                ax.set_ylim(0.150,0.1515)
                ax.set_ylabel('seconds')

                # customise the tick locations
                locator = AutoDateLocator(minticks=3, maxticks=9)
                ax.xaxis.set_major_formatter(ObsPyAutoDateFormatter(locator))
                ax.xaxis.set_major_locator(locator)

                # customise the minor tick locations
                locator = MonthLocator()
                ax.xaxis.set_minor_locator(locator)

            # tight_layout to set the space around the plot properly
            # plt.tight_layout()
            # plt.subplots_adjust(left=0.08, right=0.96)
            # plt.subplots_adjust(left=0.08, right=0.96)
            plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.93)
            if save_pdf:
            # print(outfile)
                plt.savefig(outfile)
            else:
                plt.show()
    else:
        print('No streams found.')

def plot_from_wildcards(
    stream_filename='',
    channels=['ATT'],
    merge_locations=False,
    time_diffs=True,
    save_pdf=False,
    outfile=None):
    '''
    Read streams in from wildcards
    Calls plot_from_stream_all_stations()
    '''


    # nrow = ncol = 2
    # a = []
    # fig, axs = plt.subplots(nrows=nrow, ncols=ncol)
    # for i, row in enumerate(axs):
    #     for j, ax in enumerate(row):
    #         a.append(ax)
    #
    # for i, ax in enumerate(a):
    #     ax.set_ylabel(str(i))
    #
    # exit()

    # check that the file pattern exists
    if not glob.glob(stream_filename):
        msg = 'view.py cannot find file: {}'.format(stream_filename)
        print(msg)
        logging.info(msg)
        return()

    stream = read(stream_filename)
    if len(channels) == 1:
        stream = stream.select(channel=channels[0])
    else:
        for tr in stream:
            if tr.stats.channel not in channels:
                stream.remove(tr)

    stream = stream.sort(keys=['starttime'])


    plot_from_stream_all_stations(
        stream=stream,
        merge_locations=merge_locations,
        time_diffs=time_diffs,
        save_pdf=save_pdf,
        outfile=outfile)


def plot_from_2_dirs(
    stations=['S11','S12','S14','S15','S16'],
    channels=['AFR','_FR'],
    file_dir1='.',
    file_dir2='.',
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T:21:00.000000Z'),
    merge_locations=False,
    plot_type='normal',
    time_interval=3,
    read_gzip=True,
    show=True,
    save_pdf=False):
    '''
    Read in files between the start and end times and plot.
    Calls plot_from_stream()
    '''

    if plot_type == 'dayplot':
        time_interval = timedelta(hours=24)
    else:
        time_interval = timedelta(hours=time_interval)

    for station in stations:

        file_dir1_station =  os.path.join(file_dir1, station)
        file_dir2_station =  os.path.join(file_dir2, station)

        start = start_time
        while start < end_time:

            stream = Stream()

            # use the same filename
            if plot_type == 'dayplot':
                stream_filename = '%s*_%s*.MINISEED' % (start.strftime("%Y-%m-%dT"), station)
            elif plot_type == 'normal' and time_interval == timedelta(hours=24):
                stream_filename = '%s*_%s*.MINISEED' % (start.strftime("%Y-%m-%dT"), station)
            else:
                stream_filename = '%s_%s*.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)

            if read_gzip:
                stream_filename = '%s.gz' % (stream_filename)
            stream1_filename = os.path.join(file_dir1_station, stream_filename)
            stream2_filename = os.path.join(file_dir2_station, stream_filename)

            if save_pdf:
                outfile = '%s_%s.pdf' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
                outfile = os.path.join(file_dir2_station, outfile)
            else:
                outfile = None

            print(stream1_filename)
            # check that the file pattern exists
            if not glob.glob(stream1_filename):
                msg = 'view.py cannot find file: {}'.format(stream1_filename)
                print(msg)
                logging.info(msg)
            else:
                stream += read(stream1_filename)

            # check that the file pattern exists
            if not glob.glob(stream2_filename):
                msg = 'view.py cannot find file: {}'.format(stream2_filename)
                print(msg)
                logging.info(msg)
            else:
                stream += read(stream2_filename)

            if len(stream) > 0:
                # select for this station
                stream = stream.select(station=station)
                stream = stream.sort(keys=['starttime'])

                plot_from_stream(
                    stream=stream,
                    stations=[station],
                    channels=channels,
                    merge_locations=merge_locations,
                    plot_type=plot_type,
                    time_interval=time_interval,
                    outfile=outfile,
                    show=show)

            # increment the time interval
            start += time_interval
