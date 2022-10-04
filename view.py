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
from datetime import datetime, timedelta
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
from matplotlib.dates import AutoDateLocator, MonthLocator, date2num, num2date, YearLocator, epoch2num
from matplotlib.collections import LineCollection
from dateutil.rrule import MONTHLY, YEARLY
import calendar
from math import floor
import pandas as pd

from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, Trace, Stats, read
from obspy.imaging.util import ObsPyAutoDateFormatter
from urllib.request import Request, build_opener, HTTPError
from pdart.csv_join_work_tapes import make_dir
from pdart.util import remove_negative_ones


import matplotlib
import matplotlib.pyplot as plt


global DELTA
DELTA = 0.1509
FULL_DELTA = 16/106
SECONDS_PER_DAY = 3600.0 * 24.0

INVALID_99 = -99.0
# 
# def plot_from_stream(
#     stream=None,
#     stations=['S11','S12','S14','S15','S16'],
#     channels=['AFR'],
#     merge_locations=False,
#     plot_type='normal',
#     time_interval=3,# removed
#     outfile=None,
#     show=True):
#     '''
#     View the traces by passing in a stream.
#     '''
# 
# 
#     # make a copy because we modify the trace
#     new_stream = stream.copy()
# 
#     for station in stations:
# 
#         # view by station (make a copy because we modify the trace)
#         station_stream = new_stream.select(station=station)
# 
#         # if view_traces_to_remove:
#         #     station_stream = station_stream.sort(keys=['starttime'])
#         #     station_stream = discard_short_traces(station_stream)
# 
#         for i, tr in enumerate(station_stream):
#             # if i == 0:
#             #     starttime0 = tr.stats.starttime
#             # if view_traces_to_remove:
#             #     if tr.stats.channel == '_FR':
#             #         msg = '{} {} {} {} {} {}'.format(tr.stats.network,
#             #           tr.stats.station, tr.stats.location, tr.stats.starttime,
#             #           tr.stats.endtime, tr.stats.npts)
#             #         logging.info(msg)
#             #         print('INFO:root:######:{}'.format(msg))
# 
#             if tr.stats.channel not in channels:
#                 station_stream.remove(tr)
# 
#             if merge_locations:
#                 tr.stats.location = ''
# 
#             # if tr.stats.channel in ('ATT'):
#             #     obs_delta = (tr.data[-1] - tr.data[0])/(tr.stats.npts - 1)
#                 # print('obs delta ', obs_delta, tr.id, tr.stats.starttime)
# 
#             if tr.stats.channel == 'ATT' or tr.stats.channel == 'AT0':
#                 # print(UTCDateTime()
#                 # this one shows a constant(ish) slope
#                 # tr.data = tr.data - starttime0.timestamp
#                 # this one shows the drift between times() and the timestamp
#                 # (should be only a few 10ths of second in 24 hours)
#                 tr.data = tr.data - tr.stats.starttime.timestamp - tr.times()
# 
#         # print("Sometimes it doesn't merge - I don't know why")
#         # stream.merge()
#         # stream.plot(size=(1200,600),method='full')
# 
#         if len(station_stream) > 0:
#             if plot_type == 'dayplot':
#                 if show:
#                     station_stream.plot(type ='dayplot',time_offset=0,show_y_UTC_label=False,
#                       vertical_scaling_range=300, outfile=None)
#                 if outfile is not None:
#                     station_stream.plot(type ='dayplot',time_offset=0,show_y_UTC_label=False,
#                       vertical_scaling_range=300,outfile=outfile)
#             else:
#                 if show:
#                     station_stream.plot(size=(1200,600),method='full',equal_scale=False, outfile=None)
#                 if outfile is not None:
#                     station_stream.plot(size=(1200,600), equal_scale=False, outfile=outfile)
# 
#             plt.close('all')

# def plot_from_file(
#     stations=['S11','S12','S14','S15','S16'],
#     channels=['AFR'],
#     file_dir='.',
#     out_dir='.',
#     start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
#     end_time=UTCDateTime('1977-09-30T21:00:00.000000Z'),
#     merge_locations=False,
#     plot_type='normal',
#     time_interval=3,
#     read_gzip=True,
#     save_pdf=False,
#     view_traces_to_remove=False):
#     '''
#     Read in files between the start and end times and plot.
#     Calls plot_from_stream()
#     '''
#
#     if view_traces_to_remove:
#         log_filename = 'logs/view_traces_to_remove.log'
#         # logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
#         logging.basicConfig(filename=log_filename, filemode='w', level=logging.INFO)
#
#     if plot_type == 'dayplot':
#         time_interval = timedelta(hours=24)
#     else:
#         time_interval = timedelta(hours=time_interval)
#
#     start = start_time
#     while start < end_time:
#         stream = Stream()
#         for channel in channels:
#             for station in stations:
#                 year = start.year
#                 julday = start.julday
#
#                 print(plot_type, time_interval)
#
#                 # use the same filename
#                 # if plot_type == 'dayplot':
#                 #     stream_filename = find_filename_date(station, '*', channel, start)
#                 # elif plot_type == 'normal' and time_interval == timedelta(hours=24):
#                 #     stream_filename = find_filename_date(station, '*', channel, start)
#                 # else:
#                 #     stream_filename = find_filename_date(station, '*', channel, start)
#
#                 stream_filename = find_filename_date(station, '*', channel, start)
#
#                 stream_directory = find_dir(file_dir,year,station,channel)
#
#                 if read_gzip:
#                     stream_filename = '%s.gz' % (stream_filename)
#
#                 stream_filename = os.path.join(file_dir,str(year),'XA', station,channel,stream_filename)
#
#                 print(stream_filename)
#
#                 # check that the file pattern exists
#                 if not glob.glob(stream_filename):
#                     msg = 'view.py cannot find file: {}'.format(stream_filename)
#                     print(msg)
#                     logging.info(msg)
#                     # increment the time interval
#                     start += time_interval
#                     continue
#
#                 stream += read(stream_filename)
#                 # print(stream)
#
#
#         if save_pdf:
#             channel_string = "_".join(channels)
#             outfile = '%s_%s_%s.pdf' % (start.strftime("%Y-%m-%dT%H:%M:%S"),
#               station, channel_string)
#             outfile = os.path.join(out_dir, outfile)
#         else:
#             outfile = None
#
#         stream = stream.trim(start,start+time_interval)
#         stream = stream.sort(keys=['starttime'])
#
#         plot_from_stream(
#             stream=stream,
#             stations=[station],
#             channels=channels,
#             merge_locations=merge_locations,
#             plot_type=plot_type,
#             time_interval=time_interval,
#             outfile=outfile,
#             view_traces_to_remove=view_traces_to_remove)
#
#         # increment the time interval
#         start += time_interval

# def plot_from_wildcards(
#     stream_filename='',
#     channels=['ATT'],
#     merge_locations=False,
#     time_diffs=True,
#     save_pdf=False,
#     outfile=None):
#     '''
#     Read streams in from wildcards
#     Calls plot_from_stream_all_stations()
#     '''
#
#     # check that the file pattern exists
#     if not glob.glob(stream_filename):
#         msg = 'view.py cannot find file: {}'.format(stream_filename)
#         print(msg)
#         logging.info(msg)
#         return()
#
#     stream = read(stream_filename)
#     if len(channels) == 1:
#         stream = stream.select(channel=channels[0])
#     else:
#         for tr in stream:
#             if tr.stats.channel not in channels:
#                 stream.remove(tr)
#
#     stream = stream.sort(keys=['starttime'])
#
#
#     plot_from_stream_all_stations(
#         stream=stream,
#         merge_locations=merge_locations,
#         time_diffs=time_diffs,
#         save_pdf=save_pdf,
#         outfile=outfile)


# def plot_from_2_dirs(
#     stations=['S11','S12','S14','S15','S16'],
#     channels=['AFR','_FR'],
#     file_dir1='.',
#     file_dir2='.',
#     start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
#     end_time=UTCDateTime('1977-09-30T21:00:00.000000Z'),
#     merge_locations=False,
#     plot_type='normal',
#     time_interval=3,
#     read_gzip=True,
#     show=True,
#     save_pdf=False):
#     '''
#     Read in files between the start and end times and plot.
#     Calls plot_from_stream()
#     '''
#
#     if plot_type == 'dayplot':
#         time_interval = timedelta(hours=24)
#     else:
#         time_interval = timedelta(hours=time_interval)
#
#     for station in stations:
#
#         file_dir1_station =  os.path.join(file_dir1, station)
#         file_dir2_station =  os.path.join(file_dir2, station)
#
#         start = start_time
#         while start < end_time:
#
#             stream = Stream()
#
#             # use the same filename
#             if plot_type == 'dayplot':
#                 stream_filename = '%s*_%s*.MINISEED' % (start.strftime("%Y-%m-%dT"), station)
#             elif plot_type == 'normal' and time_interval == timedelta(hours=24):
#                 stream_filename = '%s*_%s*.MINISEED' % (start.strftime("%Y-%m-%dT"), station)
#             else:
#                 stream_filename = '%s_%s*.MINISEED' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
#
#             if read_gzip:
#                 stream_filename = '%s.gz' % (stream_filename)
#             stream1_filename = os.path.join(file_dir1_station, stream_filename)
#             stream2_filename = os.path.join(file_dir2_station, stream_filename)
#
#             if save_pdf:
#                 outfile = '%s_%s.pdf' % (start.strftime("%Y-%m-%dT%H:%M:%S"), station)
#                 outfile = os.path.join(file_dir2_station, outfile)
#             else:
#                 outfile = None
#
#             print(stream1_filename)
#             # check that the file pattern exists
#             if not glob.glob(stream1_filename):
#                 msg = 'view.py cannot find file: {}'.format(stream1_filename)
#                 print(msg)
#                 logging.info(msg)
#             else:
#                 stream += read(stream1_filename)
#
#             # check that the file pattern exists
#             if not glob.glob(stream2_filename):
#                 msg = 'view.py cannot find file: {}'.format(stream2_filename)
#                 print(msg)
#                 logging.info(msg)
#             else:
#                 stream += read(stream2_filename)
#
#             if len(stream) > 0:
#                 # select for this station
#                 stream = stream.select(station=station)
#                 stream = stream.sort(keys=['starttime'])
#
    #                 plot_from_stream(
#                     stream=stream,
#                     stations=[station],
#                     channels=channels,
#                     merge_locations=merge_locations,
#                     plot_type=plot_type,
#                     time_interval=time_interval,
#                     outfile=outfile,
#                     show=show)
#
#             # increment the time interval
#             start += time_interval


# def moving_average(a, n=3):
# # jerry from https://stackoverflow.com/users/1984065/jerry
#
#     if isinstance(a,np.ma.MaskedArray):
#         ret = np.cumsum(a.filled(0))
#         ret[n:] = ret[n:] - ret[:-n]
#         counts = np.cumsum(~a.mask)
#         counts[n:] = counts[n:] - counts[:-n]
#         ret[~a.mask] /= counts[~a.mask]
#         ret[a.mask] = np.nan
#     else:
#         ret = np.cumsum(a,dtype=float)
#         ret[n:] = ret[n:]-ret[:-n]
#         return ret[-1:]/n
#
#     return ret

def read_solar_presence(
    stations=['S11','S12','S14','S15','S16'],
    solar_presence_dir='../solar_presence/'
    ):

    solar_presence_times = []

    for station in stations:
        infile = '%s.txt' % (station)
        infile =  os.path.join(solar_presence_dir, infile)

        with open(infile, 'r') as f:
            parsing = False
            for line in f.readlines():
                if line.startswith('$$SOE'):
                    parsing = True
                elif line.startswith('$$EOE'):
                    break
                elif parsing:
                    # # 1977-Dec-05 05:44 2443482.738888889  s  271.3170  -0.3015
                    s_date, s_time, _, solar_presence, _, _ = line.split()

                    # parse the date
                    s_datetime = "%s %s" % (s_date, s_time)
                    s_datetime = datetime.strptime(s_datetime, "%Y-%b-%d %H:%M")
                    s_datetime = UTCDateTime(s_datetime)

                    if len(solar_presence) == 2:
                        solar_presence = solar_presence[1:2]
                    # only use rising and sunset
                    if solar_presence in ('r', 's'):
                        t = (station, s_datetime, solar_presence)
                        solar_presence_times.append(t)

            f.closed

    return solar_presence_times


# def read_time_delays(
#     stations=['S11','S12','S14','S15','S16'],
#     time_delays_dir='../time_delays/'
#     ):
# 
#     time_delays = []
# 
#     for station in stations:
#         infile = '%s.txt' % (station)
#         infile =  os.path.join(time_delays_dir, infile)
# 
#         with open(infile, 'r') as f:
#             parsing = False
#             for line in f.readlines():
#                 if line.startswith('$$SOE'):
#                     parsing = True
#                 elif line.startswith('$$EOE'):
#                     break
#                 elif parsing:
# 
#  # 1969-Jan-01 00:00:00.000 2440222.500000000 *x    0.022505
#                     s_date, s_time, _, _, time_delay = line.split()
# 
#                     # parse the date
#                     s_datetime = "%s %s" % (s_date, s_time)
#                     s_datetime = datetime.strptime(s_datetime, "%Y-%b-%d %H:%M:%S.%f")
# 
#                     time_dalay_date = UTCDateTime(s_datetime)
# 
#                     time_delay = round(float(time_delay),6)
#                     # put the time_delay into seconds
#                     time_delay = round(time_delay * 60.,5)
# 
#                     t = (station, time_dalay_date, time_delay)
# 
#                     time_delays.append(t)
# 
#             f.closed
# 
#     return time_delays

# def elapsed_time(stream):
#     # calculate elapsed times
# 
#     e_time = 3600
#     # approx samples per hour 23,850
#     samples_per_hour = floor(e_time / DELTA)
# 
#     for tr in stream:
#         tr.stats.location = ''
#     stream = stream.merge()
# 
#     if stream.count() > 1:
#         print('Too many streams')
#         exit()
#     else:
#         tr = stream[0]
# 
#     starttime = tr.stats.starttime
#     # check that this is approximately on the hour
#     if starttime.minute > 1:
#         if starttime.hour > 23:
#             return
#         else:
#             trim_starttime = UTCDateTime(starttime.year, starttime.month, starttime.day, starttime.hour+1, 0)
#             tr.trim(starttime=trim_starttime)
# 
#     # slice the array (start:stop:step notation )
#     tr.data = tr.data[0::samples_per_hour]
#     tr.stats.delta = e_time
#     tr.stats.channel = 'ASI'
# 
#     time_diff = np.diff(tr.data)
# 
#     tr.data = time_diff/samples_per_hour
#     tr.stats.starttime = tr.stats.starttime + timedelta(minutes=30)
# 
#     if tr.stats.npts < 1:
#         stream.remove(tr)
#     elif isinstance(tr.data,np.ma.MaskedArray) and tr.data.count() < 1:
#         stream.remove(tr)
# 
#     return stream

def plot_overview_from_stream(
    stream,
    channel='ASI',
    stations=['S11','S12','S14','S15','S16'],
    ylim = (480 - 100, 480 + 100),
    ylim_masking = False, # mask values outside the ylimits 
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z'),
    plot_solar_presence=True,
    solar_presence_dir='../Electronic_Supplement/files/solar_presence/',
    # plot_time_delays=False,
    # time_delays_dir='../time_delays/',
    merge_locations=False,
    save_pdf=False,
    outfile=None):
    '''
    Plot an overview
    '''


    stream = stream.sort(keys=['starttime'])

    if merge_locations:
        for tr in stream:
            tr.stats.location = ''

    # sometimes small (probably non existant) differences in the floats seem to throw errors
    for tr in stream:
        tr.stats.delta = round(tr.stats.delta,8)

    # remove the masks
    for tr in stream: 
        if tr.stats.channel in ('ASI'):
            tr.data = np.ma.masked_values(tr.data,-99.0)


    print(solar_presence_dir)
    if plot_solar_presence:
        solar_presence_times = read_solar_presence(solar_presence_dir=solar_presence_dir)

    # if plot_time_delays:
    #     time_delays = read_time_delays(time_delays_dir=time_delays_dir)

    # plot using the sampling rate of the mid-period channels (4 x smaller
    # than 'ATT')
    for tr in stream: 
        tr.data = tr.data/4.0

        if ylim_masking:
            tr.data = ma.masked_where(tr.data <= ylim[0], tr.data)
            tr.data = ma.masked_where(tr.data >= ylim[1], tr.data)


    stream.merge()

    if len(stream) > 0:

        # check which stations were found
        # stations = sorted(stations)
        # for station in stations:
        #     # make new list with just the ones with data
        #     stations[:] = [stat for stat in stations if (stream.select(station=stat).count() > 0)]

        fig = stream.plot(handle=True,size=(1200,600),show=False)
    # ,method='full')
        # clear the suptitle which has already been set
        fig.suptitle('')
        axes = fig.get_axes()
        if len(axes) == 3:
            height = 6 *0.75
        else:
            height = 6
        fig.set_size_inches(12, height)
        # for ax in axes:
        #     ax.axhline(y=DELTA, color='r')
        #
        #     if ylim is None:
        #         if channel == 'ASI':
        #             # ax.set_ylim(0.150925,0.150955)
        #             ax.set_ylim(0.116,0.189)
        #             # pass
        #     else:
        #         ax.set_ylim(ylim)
        #     # don't use scientific notation for the y axis
        #     ax.ticklabel_format(useOffset=False,axis='y')
        # plt.show()

        plt.ylabel('s')
        for i, ax in enumerate(axes):
            # view unaveraged time diffs
            # ax.set_ylim(0.14,0.16)
            # per frame time diffs

            # print(handles)
            # print(labels)

            if i == 0:
                if channel == 'ASI':
                    ax.set_title('Average Sampling Interval')
                else:
                    ax.set_title('Overview - {}'.format(channel))


            if channel == 'ASI':
                ax.axhline(y=FULL_DELTA, color='r')

            if ylim is None:
                if channel == 'ASI':
                    ax.set_ylim(0.150925,0.150955)
                    # ax.set_ylim(0.116,0.189)
            else:
                ax.set_ylim(ylim)
            # don't use scientific notation for the y axis
            ax.ticklabel_format(useOffset=False,axis='y')
            if channel == 'ASI':
                ax.set_ylabel('seconds')
            else:
                ax.set_ylabel('DU')
            # customise the tick locations
            locator = AutoDateLocator(minticks=3, maxticks=9)
            ax.xaxis.set_major_formatter(ObsPyAutoDateFormatter(locator))
            ax.xaxis.set_major_locator(locator)

            # customise the minor tick locations
            locator = MonthLocator()
            ax.xaxis.set_minor_locator(locator)

            # find the x axis (note the axis uses matplotlib dates)
            (xlim_min, xlim_max) = ax.get_xlim()

            # if plot_time_delays:
            #     # get the time_delays for this station
            #     station_time_delay_dates = []
            #     station_time_delays = []
            #     for (time_delay_station, time_delay_date, time_delay) in time_delays:
            #         if time_delay_station == stations[i]:
            #             # change UTCDateTime to the right format
            #             station_time_delay_dates.append(date2num(time_delay_date.datetime))
            #             station_time_delays.append(time_delay)

            if plot_solar_presence:
                # using the labels - we can find which station is plotted on
                # this axis
                for child in ax.get_children():
                    if isinstance(child, matplotlib.text.Text):
                        label = child.get_text()
                        if label != '' and label.split('.')[0] == 'XA':
                            ax_station = label.split('.')[1]

                            # print(solar_presence_times)

                            # first get the ones that apply
                            station_solar_presence_times = []
                            for (solar_presence_station, s_datetime, solar_presence) in solar_presence_times:
                                # print(solar_presence_station, ax_station)
                                if solar_presence_station == ax_station:
                                    # solar_mark = date2num(s_datetime.datetime)
                                    # print(solar_mark)
                                    # if solar_mark >= xlim_min and solar_mark <= xlim_max:
                                    station_solar_presence_times.append((solar_presence_station, s_datetime, solar_presence))
                                # else:
                                #     print(solar_presence_station)

                            # print(station_solar_presence_times)

                            for i, (solar_presence_station, s_datetime, solar_presence) in enumerate(station_solar_presence_times):
                                start_shade = None
                                end_shade = None
                                # sunrise
                                if solar_presence == 'r':
                                    end_shade = date2num(s_datetime.datetime)
                                    # find the last time the sun set (or the begining)
                                    if i == 0:
                                        start_shade = xlim_min
                                    else:
                                        start_shade = station_solar_presence_times[i-1][1]
                                        start_shade = date2num(start_shade.datetime)

                                    ax.axvspan(start_shade, end_shade, facecolor='b', alpha=0.3)

                            # deal with the last one
                            start_shade = None
                            end_shade = None
                            if solar_presence == 's':
                                start_shade = date2num(s_datetime.datetime)
                                end_shade = xlim_max
                                ax.axvspan(start_shade, end_shade, facecolor='b', alpha=0.3)

            # if plot_time_delays:
            #     ax2 = ax.twinx()
            #     ax2.plot(station_time_delay_dates, station_time_delays, 'g-')
            #     ax2.set_xlim((xlim_min, xlim_max))
            #     ax2.set_ylabel('seconds')
            #     ax2.set_ylim(0.95,1.75)

        # tight_layout to set the space around the plot properly
        # plt.tight_layout()
        xlim_min = date2num(start_time.datetime)
        xlim_max = date2num(end_time.datetime)
        plt.xlim(xlim_min, xlim_max)

        plt.subplots_adjust(left=0.09, bottom=0.08, right=0.95, top=0.93)
        if save_pdf:
            print(outfile)
            plt.savefig(outfile)
            plt.show()
        else:
            plt.show()
    else:
        print('No streams found.')

def plot_overview_from_file(
    stations=['S11','S12','S14','S15','S16'],
    channels=['ASI'],
    ylim=None,
    ylim_masking = False, # mask values outside the ylimits 
    overview_top_level_dir='.',
    pdf_dir='.',
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z'),
    read_gzip=True,
    plot_solar_presence=True,
    solar_presence_dir='../Electronic_Supplement/files/solar_presence/',
    # plot_time_delays=False,
    # time_delays_dir='../time_delays/',
    merge_locations=False,
    save_pdf=False):

    '''
    Calls the method which plots the overview - including the Sampling Interval
    '''

    # check that the overall directory exists
    if not os.path.exists(pdf_dir):
        msg = ("The directory {} doesn't exist".format(pdf_dir))
        raise IOError(msg)

    # original_channel = None
    # if channels == ['ASI']:
    #     channels = ['ATT']
    #     original_channel = 'ASI'

    time_interval = timedelta(hours=24)
    read_stream = Stream()
    start = start_time
    print(start_time, end_time) 
    while start < end_time:
        for channel in channels:
            for station in stations:

                year = start.year
                julday = start.julday

                file_dir = find_dir_upper(top_level_dir=overview_top_level_dir,station=station,starttime=start)
                av_filename = find_filename_date_upper(station, '*', 'ASI', start)
                filename = os.path.join(file_dir,av_filename)

                if read_gzip:
                    filename = '%s.gz' % (filename)

                # check that the file pattern exists
                if not glob.glob(filename):
                    msg = 'view.py cannot find file: {}'.format(filename)
                    print(msg)
                    logging.info(msg)
                    # # increment the time interval
                    # start += time_interval
                    continue

                stream = read(filename)
                read_stream += stream

        # increment the time interval
        start += time_interval


    outfile = 'overview_%s_%s_%s.png' % (channel, start_time.strftime("%Y-%m-%d"), end_time.strftime("%Y-%m-%d"))
    outfile = os.path.join(pdf_dir,outfile)

    

    # if original_channel is not None:
    #     for tr in read_stream:
    #         tr.stats.channel = 'ASI'
    #         if len(tr) > 2:
    #             tr.data = np.diff(tr.data)
    #         else:
    #             read_stream.remove(tr)
    #
    #     channel = original_channel

    # print(read_stream)
    plot_overview_from_stream(
      read_stream,
      channel=channel,
      stations=stations,
      ylim=ylim,
      ylim_masking=ylim_masking,
      start_time=start_time,
      end_time=end_time,
      plot_solar_presence=plot_solar_presence,
      solar_presence_dir=solar_presence_dir,
      # plot_time_delays=plot_time_delays,
      # time_delays_dir=time_delays_dir,
      merge_locations=merge_locations,
      save_pdf=save_pdf,
      outfile=outfile)

# # used 17 October
# def stream_from_directory(
#     top_level_dir,
#     start_time,
#     stations,
#     channels,
#     merge_locations=True,
#     end_time=None):
#
#     # helper method to read files from the directory
#     # and return a stream
#
#     if end_time is None:
#         end_time = start_time + timedelta(hours=3)
#
#     time_interval = timedelta(hours=24)
#
#     return_stream = Stream()
#
#     start = start_time
#     while start < end_time:
#         for station in stations:
#             for channel in channels:
#                 filename = find_filename_date(station, '*', channel, start)
#                 filename = "{}.gz".format(filename)
#                 directory = find_dir(top_level_dir,start.year,station,channel)
#                 filename = os.path.join(directory,filename)
#
#                 stream = read(filename)
#                 stream.trim(starttime=start_time, endtime=end_time)
#
#                 if merge_locations:
#                     for tr in stream:
#                         tr.stats.location = ''
#
#                 return_stream += stream
#
#         # increment the time interval
#         start += time_interval
#     return_stream.merge()
#     return_stream = return_stream.sort(keys=['channel'])
#
#     return return_stream

# changed to new structure 25-08-2021 
# note potential errors with upper/lower and .gz 
def stream_from_directory_new(
    top_level_dir,
    start_time,
    stations,
    channels,
    merge_locations=False,
    end_time=None):

    # helper method to read files from the directory
    # and return a stream

    if end_time is None:
        # as a default, retrieve the whole day 
        end_time = start_time + timedelta(hours=24)

    time_interval = timedelta(hours=24)

    return_stream = Stream()

    start = start_time
    while start < end_time:
        for station in stations:
            for channel in channels:
                filename = find_filename_date_upper(station.upper(), '*', channel.upper(), start)
                # filename = "{}.gz".format(filename)
                # print(filename)
                directory = find_dir_upper(top_level_dir,station,start)
                # print(directory)
                filename = os.path.join(directory,filename)
                if glob.glob(filename):
                    print(filename)


                    stream = read(filename)
                    # stream.trim(starttime=start_time, endtime=end_time)

                    if merge_locations:
                        for tr in stream:
                            tr.stats.location = ''

                        # print('temporary measure !!!')    
                        # for tr in stream:
                        #     tr.stats.delta = DELTA*4

                    return_stream += stream

                else:
                    print('File not found ', filename)



        # increment the time interval
        start += time_interval

    print(return_stream)
    for tr in return_stream:
        print(tr.stats.sampling_rate, tr.stats.delta)

    return_stream.merge()
    return_stream = return_stream.sort(keys=['channel'])

    return_stream.trim(starttime=start_time,endtime=end_time)

    return return_stream

def find_filename_date_lower(station, location, channel, datetime):
    # xa.s12..afr.1973.324.0.mseed
    # for writing, an empty option is possible
    # for reading, an empty option or a wildcard (*) is possible
    return '%s.%s.%s.%s.%s.%03d.*.mseed' % ('xa',station, location, channel,
      str(datetime.year), datetime.julday)

def find_filename_date_upper(station, location, channel, datetime):
    # xa.s12..afr.1973.324.0.mseed
    # for writing, an empty option is possible
    # for reading, an empty option or a wildcard (*) is possible
    return '%s.%s.%s.%s.%s.%03d.*.MSEED' % ('XA',station, location, channel,
      str(datetime.year), datetime.julday)

def find_dir_lower(top_level_dir,station,starttime):
    # /Users/cnunn/lunar_data/PDART_CONTINUOUS_MAIN_TAPES/s16/1972/366
    year = str(starttime.year)
    day = str('{:03}'.format(starttime.julday))
    return os.path.join(top_level_dir,station,year,day)

def find_dir_upper(top_level_dir,station,starttime):
    # /Users/cnunn/lunar_data/PDART_CONTINUOUS_MAIN_TAPES/s16/1972/366
    year = str(starttime.year)
    day = str('{:03}'.format(starttime.julday))
    return os.path.join(top_level_dir,station,year,day)

def save_availability_coarse(
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z'),
    channels=['MHZ','MH1','MH2','SHZ'],
    stations=['S11','S12','S14','S15','S16'],
    file_dir='.',
    pdf_dir='.',
    merge_locations=True,
    read_gzip=True,
    save_pdf=False,
    outfile=None):
    '''
    (Coarsely) save the availability for the all the stations and channels.

    See also save_availability() which contains more detail. 
    Then call plot_availability()
    For the figure, and without much detail. 
    '''

    # stations=['S12']
    # channels=['MH2']
    # channels=['ASI']
    # stations=['S12']

    time_interval = timedelta(hours=24)

    # stations = [sta.lower() for sta in stations]
    # channels = [cha.lower() for cha in channels]

    segs = []


    start = start_time
    while start < end_time:

        # year = start.year
        # julday = start.julday

        filename = find_filename_date_upper('*', '*', '*', start)
        directory = find_dir_upper(file_dir,'*',start)
        filename = os.path.join(directory,filename)

        if read_gzip:
            filename = '%s.gz' % (filename)

        # check that the file pattern exists
        if not glob.glob(filename):
            msg = 'view.py cannot find file: {}'.format(filename)
            print(msg)
            logging.info(msg)
            # increment the time interval
            start += time_interval
            continue

        stream = read(filename)

        for station in stations:
            station_upper = station.upper()
            station_stream = stream.select(station=station)
            # continue
            for channel in channels:
                channel_upper = channel.upper()
                channel_stream = station_stream.select(channel=channel)
                for location in ['00','01','']: 
                    location_stream = channel_stream.select(location=location)
                    location_stream.merge()
                    for tr in location_stream:  
                        seg = tr.stats.starttime,tr.stats.endtime,station_upper,channel_upper,location
                        segs.append(seg)


        # increment the time interval
        start += time_interval



    if outfile is None:
        'availability1.csv'
    with open(outfile,'w') as out:
        csv_out=csv.writer(out)
        csv_out.writerow(['t1','t2','station','channel','location'])
        for row in segs:
            csv_out.writerow(row)            


def save_availability(
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z'),
    channels=['MHZ','MH1','MH2','SHZ'],
    stations=['S11','S12','S14','S15','S16'],
    file_dir='.',
    pdf_dir='.',
    merge_locations=True,
    read_gzip=True,
    save_pdf=False,
    outfile=None):
    '''
    Save the availability for the all the stations and channels.

    See also save_availability_coarse() which has less fine detail. 
    Then call plot_availability()
    For the figure, and without much detail. 
    '''

    # stations=['S12']
    # channels=['MH2']
    # channels=['ASI']
    # stations=['S12']

    time_interval = timedelta(hours=24)

    stations = [sta.lower() for sta in stations]
    channels = [cha.lower() for cha in channels]

    segs = []


    start = start_time
    while start < end_time:

        # year = start.year
        # julday = start.julday

        filename = find_filename_date_lower('*', '*', '*', start)
        directory = find_dir_lower(file_dir,'*',start)
        filename = os.path.join(directory,filename)

        if read_gzip:
            filename = '%s.gz' % (filename)

        # check that the file pattern exists
        if not glob.glob(filename):
            msg = 'view.py cannot find file: {}'.format(filename)
            print(msg)
            logging.info(msg)
            # increment the time interval
            start += time_interval
            continue

        stream = read(filename)

        for station in stations:
            station_upper = station.upper()
            station_stream = stream.select(station=station)
            # continue
            for channel in channels:
                channel_upper = channel.upper()
                channel_stream = station_stream.select(channel=channel)
                for location in ['00','01','']: 
                    location_stream = channel_stream.select(location=location)
                    for tr in location_stream:
                        # replace any gaps less than one minute
                        # make a pandas data series
                        data_series = pd.Series(tr.data)
                        # replace the -1 values with a null
                        data_series.replace(-1.0, pd.NA, inplace=True)
                        # make a linear interpolation to fill gaps up to about one minute
                        if channel == 'SHZ':
                            limit = 3180
                        else:
                            limit = 398
                        
                        data_series.interpolate(method='linear', axis=0, limit=limit, 
                          inplace=True, limit_direction=None, limit_area='inside', 
                          downcast=None)
                        # also get rid of the first gap (if there is one)
                        data_series.interpolate(method='backfill', axis=0, limit=1, 
                          inplace=True, limit_direction=None, limit_area=None, 
                          downcast=None)

                        # replace any nans - be careful with the datatype again
                        data_series.fillna(-1, inplace=True)
                        tr.data=data_series.to_numpy(dtype='int')
                        tr.data = np.ma.masked_where(tr.data == -1, tr.data)
                    # location_stream.print_gaps()

                    # need to check that actual gaps are recorded 
                    location_stream = location_stream.split()
                    for tr in location_stream:  
                        seg = tr.stats.starttime,tr.stats.endtime,station_upper,channel_upper,location
                        segs.append(seg)


        # increment the time interval
        start += time_interval



    if outfile is None:
        'availability1.csv'
    with open(outfile,'w') as out:
        csv_out=csv.writer(out)
        csv_out.writerow(['t1','t2','station','channel','location'])
        for row in segs:
            csv_out.writerow(row)            


# def save_availability(
#     start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
#     end_time=UTCDateTime('1977-09-30T21:00:00.000000Z'),
#     channels=['MHZ','MH1','MH2','SHZ'],
#     stations=['S11','S12','S14','S15','S16'],
#     # locations=['','00','01'],
#     locations=['01'],
#     file_dir='.',
#     pdf_dir='.',
#     merge_locations=True,
#     read_gzip=True,
#     save_pdf=False,
#     outfile=None):
#     '''
#     Save the availability for the all the stations and channels.
#     Then call plot_availability()
#     For the figure, and without much detail. 
#     '''
# 
#     # stations=['S12']
#     # channels=['MH2']
#     # channels=['ASI']
#     # stations=['S12']
# 
#     time_interval = timedelta(hours=24)
# 
#     stations = [sta.lower() for sta in stations]
#     channels = [cha.lower() for cha in channels]
# 
#     segs = []
# 
#     for channel in channels:
#         channel1 = channel.upper()
#         for station in stations:
#             station1 = station.upper()
#             for location in locations:
#                 start = start_time
#                 while start < end_time:
# 
#                     # year = start.year
#                     # julday = start.julday
# 

# 
#                     filename = find_filename_date_lower(station, location, channel, start)
#                     # print(' this is what I look for ', filename)
#                     directory = find_dir_lower(file_dir,station,start)
#                     filename = os.path.join(directory,filename)
# 
#                     if read_gzip:
#                         filename = '%s.gz' % (filename)
# 
#                     # check that the file pattern exists
#                     if not glob.glob(filename):
#                         print(channel, station)
#                         msg = 'view.py cannot find file: {}'.format(filename)
#                         print(msg)
#                         logging.info(msg)
#                         # increment the time interval
#                         start += time_interval
#                         continue
# 
# 
#                     stream = read(filename)
# 
#                     channel = stream[0].stats.channel
#                     # print(channel)
# 
#                     for tr in stream:
#                         # replace any gaps less than one minute
#                         # make a pandas data series
#                         data_series = pd.Series(tr.data)
#                         # replace the -1 values with a null
#                         data_series.replace(-1.0, pd.NA, inplace=True)
#                         # print('replace 1')
#                         # print(data_series[0:90].to_string())
#                         # make a linear interpolation to fill gaps up to about one minute
#                         if channel == 'SHZ':
#                             limit = 3180
#                         else:
#                             limit = 398
# 
#                         data_series.interpolate(method='linear', axis=0, limit=limit, 
#                           inplace=True, limit_direction=None, limit_area='inside', 
#                           downcast=None)
#                         # also get rid of the first gap (if there is one)
#                         data_series.interpolate(method='backfill', axis=0, limit=1, 
#                           inplace=True, limit_direction=None, limit_area=None, 
#                           downcast=None)
# 
#                         # print('interpolate')
#                         # print(data_series[0:90].to_string())
#                         # replace any nans - be careful with the datatype again
#                         data_series.fillna(-1, inplace=True)
#                         # print('fillna')
#                         # print(data_series[0:90].to_string())
#                         tr.data=data_series.to_numpy(dtype='int')
#                         # print('tr data')
#                         # print(tr.data[0:90])
#                         tr.data = np.ma.masked_where(tr.data == -1, tr.data)
#                         # print('tr data end')
#                         # print(tr.data[0:90])

#                     stream.print_gaps()
#                     # 

# 
#                     # need to check that actual gaps are recorded 
#                     stream = stream.split()
#                     for tr in stream:        
#                         seg = tr.stats.starttime,tr.stats.endtime,station1,channel1,location
#                         print(seg)
#                         segs.append(seg)
# 
#                     # increment the time interval
#                     start += time_interval
# 
# 
# 
#     if outfile is None:
#         'availability1.csv'
#     with open(outfile,'w') as out:
#         csv_out=csv.writer(out)
#         csv_out.writerow(['t1','t2','station','channel','location'])
#         for row in segs:
#             csv_out.writerow(row)            

            

                

    # # call_ALL_plot_availability()
    # 
    # 
    # # call_shz_save_availability_work_tapes()
    # st = read()
    # st = st[0:1]
    # tr_orig = st[0].copy()
    # for tr in st:
    #     for i in range(0,500,20):
    #         if i < len(tr.data):
    #             # no gaps after interpolation
    #             tr.data[i:i+1] = -1
    # 
    #     # tr.data= np.ma.masked_where(tr.data==-1, tr.data)
    #     data_series = pd.Series(tr.data)
    #     # data_series = pd.Series(tr.data).astype('int64')
    #     # data_series = pd.Series(data_series).astype('Int64')
    #     # be careful with the data type here 
    #     data_series.replace(-1.0, pd.NA, inplace=True)
    #     data_series.interpolate(method='linear', axis=0, limit=1, inplace=True, limit_direction=None, limit_area='inside', downcast=None)
    #     # replace any nans - be careful with the datatype again
    #     data_series.fillna(-1.0, inplace=True)
    #     tr.data=data_series.to_numpy()
    #     tr.data = np.ma.masked_where(tr.data == -1.0, tr.data)
    # 
    # 
    # # st += tr_orig
    # st.plot()
    # 
    # # put a gap in 
    # # interpolate it 
    # # get the trace back to a masked array 
    # # check that only one missing sample is interpolated 
    # 
    # 
    # 
    # # are they 
    # st.print_gaps()  

                
                # for tr in stream:
                #     if -1 in tr.data:
                #         print('-1 found')
                #     print(len(tr))
                # 
                # if channel not in 'SHZ':
                # 
                #     stream.merge(method=1, fill_value='interpolate', interpolation_samples=60)
                #     stream.print_gaps(1)
                    

            

                # print('before')
                # print(len(stream))
                # # interpolate if there are any small gaps 
                # stream.merge(method=1, fill_value='interpolate', interpolation_samples=60)
                # print('after')
                # print(len(stream))
                # # 
                # # split the stream
                # stream = stream.split()
                # 
                # print('after splitting')
                # print(len(stream))
                # 

    



def plot_availability(
    filenames=[],
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z'),
    channels=['MHZ','MH1','MH2','SHZ'],
    stations=['S11','S12','S14','S15','S16'],
    pdf_dir='.',
    merge_locations=True,
    read_gzip=True,
    save_pdf=False,
    outfile=None):
    '''
    Save the availability for the all the stations and channels.
    Then call plot_availability()
    For the figure, and without much detail. 
    '''

    import pandas as pd

    # From the Tol color palette
    # https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
    tol_dark_blue ='#341B88' 
    tol_green = '#357932'
    tol_blue_green = '#4CAA9A'
    tol_pale_blue = '#82C0E1'

    colors = [tol_green, tol_pale_blue, tol_dark_blue]
    labels = ['peaked mid-period','flat mid-period', 'short-period']

    list_line_segments = []

    df_all = None
    for filename in filenames:
        df1 = pd.read_csv(filename, dtype=str)

        if df_all is None:
            df_all = df1
        else:
            df_all = pd.concat([df_all, df1], ignore_index = False)


    print('read in files')

    
    df_all['t1'] = df_all['t1'].astype('datetime64[ns, UTC]')
    df_all['t2'] = df_all['t2'].astype('datetime64[ns, UTC]')

    df_all = df_all[df_all.station.isin(stations)]
    df_all = df_all[df_all.channel.isin(channels)]

    df_all['location'].fillna('02', inplace=True)

    for i, location in enumerate(['00','01','02']):
        color = colors[i]
        
    

        df = df_all[df_all.location == location].copy()

        t1a =  df['t1'].to_numpy()
        t1b = date2num(t1a)
        df['t1_num'] = t1b

        t2a =  df['t2'].to_numpy()
        t2b = date2num(t2a)
        df['t2_num'] = t2b

        df['y'] = pd.NA

        df['y'] = np.where((df.channel == 'SHZ') & (df.station == 'S16'), 1, df['y'])
        df['y'] = np.where((df.channel == 'MH1') & (df.station == 'S16'), 2, df['y'])
        df['y'] = np.where((df.channel == 'MH2') & (df.station == 'S16'), 3, df['y'])
        df['y'] = np.where((df.channel == 'MHZ') & (df.station == 'S16'), 4, df['y'])

        df['y'] = np.where((df.channel == 'SHZ') & (df.station == 'S15'), 6, df['y'])
        df['y'] = np.where((df.channel == 'MH1') & (df.station == 'S15'), 7, df['y'])
        df['y'] = np.where((df.channel == 'MH2') & (df.station == 'S15'), 8, df['y'])
        df['y'] = np.where((df.channel == 'MHZ') & (df.station == 'S15'), 9, df['y'])

        df['y'] = np.where((df.channel == 'SHZ') & (df.station == 'S14'), 11, df['y'])
        df['y'] = np.where((df.channel == 'MH1') & (df.station == 'S14'), 12, df['y'])
        df['y'] = np.where((df.channel == 'MH2') & (df.station == 'S14'), 13, df['y'])
        df['y'] = np.where((df.channel == 'MHZ') & (df.station == 'S14'), 14, df['y'])

        # SHZ at S12 didn't work
        df['y'] = np.where((df.channel == 'SHZ') & (df.station == 'S12'), -10, df['y'])
        df['y'] = np.where((df.channel == 'MH1') & (df.station == 'S12'), 17, df['y'])
        df['y'] = np.where((df.channel == 'MH2') & (df.station == 'S12'), 18, df['y'])
        df['y'] = np.where((df.channel == 'MHZ') & (df.station == 'S12'), 19, df['y'])

        df['y'] = np.where((df.channel == 'SHZ') & (df.station == 'S11'), 21, df['y'])
        df['y'] = np.where((df.channel == 'MH1') & (df.station == 'S11'), 22, df['y'])
        df['y'] = np.where((df.channel == 'MH2') & (df.station == 'S11'), 23, df['y'])
        df['y'] = np.where((df.channel == 'MHZ') & (df.station == 'S11'), 24, df['y'])


        df['seg_a'] = list(zip(df.t1_num, df.y))
        df['seg_b'] = list(zip(df.t2_num, df.y))
        df['seg'] = list(zip(df.seg_a, df.seg_b))
        
        segs = list(df['seg'])


        line_segments = LineCollection(segs, linewidths=4,
                                       colors=color, linestyle='solid', label=labels[i])

        list_line_segments.append(line_segments)

        


    print('finished reading files')
    print('made line segments')

    fig = plt.figure(figsize=(15, 4))
    # gs = gridspec.GridSpec(5, 1, hspace=0.001)
    # 
    # ax0 = plt.subplot(gs[0])

    plt.title('Data Availability', fontsize=16)
    ylim_min = 0
    ylim_max = 25
    plt.yticks([2.5,7.5,12.5,17.5,22.5], ['S16', 'S15', 'S14', 'S12', 'S11'], fontsize=14)
    # limits need to set after defining the ticks
    plt.ylim(ylim_min,ylim_max)
    xlim_min=date2num(start_time.datetime)
    xlim_max=date2num(end_time.datetime)
    plt.xlim(xlim_min,xlim_max)

    axes = fig.get_axes()
    axes[0].add_collection(list_line_segments[0])
    axes[0].add_collection(list_line_segments[1])
    axes[0].add_collection(list_line_segments[2])

    for i, ax in enumerate(axes):

        # customise the tick locations
        locator = AutoDateLocator(minticks=3, maxticks=9)
        ax.xaxis.set_major_formatter(ObsPyAutoDateFormatter(locator))
        ax.xaxis.set_major_locator(locator)

        # customise the minor tick locations
        locator = MonthLocator()
        ax.xaxis.set_minor_locator(locator)
        plt.setp(ax.get_xticklabels(), fontsize=14)

    an_fs = 11

    plt.annotate(xy=(1.01,1/ylim_max), s='SHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,2/ylim_max), s='MH2', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,3/ylim_max), s='MH1', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,4/ylim_max), s='MHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    plt.annotate(xy=(1.01,6/ylim_max), s='SHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,7/ylim_max), s='MH2', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,8/ylim_max), s='MH1', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,9/ylim_max), s='MHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    plt.annotate(xy=(1.01,11/ylim_max), s='SHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,12/ylim_max), s='MH2', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,13/ylim_max), s='MH1', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,14/ylim_max), s='MHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    plt.annotate(xy=(1.01,16/ylim_max), s='SHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,17/ylim_max), s='MH2', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,18/ylim_max), s='MH1', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,19/ylim_max), s='MHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    plt.annotate(xy=(1.01,21/ylim_max), s='SHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,22/ylim_max), s='MH2', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,23/ylim_max), s='MH1', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    plt.annotate(xy=(1.01,24/ylim_max), s='MHZ', fontsize=an_fs,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    # for S11 if required
    # plt.annotate(xy=(1.01,21/ylim_max), text='SHZ', fontsize=14,
    #   xycoords="axes fraction", horizontalalignment='left',
    #   verticalalignment='center')
    # plt.annotate(xy=(1.01,22/ylim_max), text='MH2', fontsize=14,
    #   xycoords="axes fraction", horizontalalignment='left',
    #   verticalalignment='center')
    # plt.annotate(xy=(1.01,23/ylim_max), text='MH1', fontsize=14,
    #   xycoords="axes fraction", horizontalalignment='left',
    #   verticalalignment='center')
    # plt.annotate(xy=(1.01,24/ylim_max), text='MHZ', fontsize=14,
    #   xycoords="axes fraction", horizontalalignment='left',
    #   verticalalignment='center')


    plt.legend(loc='lower left', fontsize=12, handletextpad=0.1)
    plt.subplots_adjust(left=0.04, right=0.96, top=0.9, bottom=0.12)

    if outfile is None: 
        outfile = 'AvailabilityX.png'

    if save_pdf:
        filename=os.path.join(pdf_dir,outfile)
        print(filename)
        plt.savefig(filename)
    else:
        plt.show()




def calc_average_sampling_interval(
    stations=['S11','S12','S14','S15','S16'],
    top_level_dir='.',
    av_dir='.',
    start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z'),
    read_gzip=True,
    write_gzip=True,
    sampling_average=360,
    # calc_moving_average=True
    ):
    '''
    Calculates the average framecounts and saves as new files.
    '''

    time_interval = timedelta(hours=24)

    # return_stream = Stream()

    # check that the overall directory exists
    if not os.path.exists(av_dir):
        msg = ("The directory {} doesn't exist".format(av_dir))
        raise IOError(msg)

    start = start_time
    while start < end_time:
        for station in stations:
            # # Note that location 00 is missing
            # for location in np.arange(1,14):
            #     location = "%02d" % (location,)
            #     filename = find_filename_date(station, location, 'ATT', start)
            #     if read_gzip:
            #         filename = "{}.gz".format(filename)
            #     directory = find_dir(file_dir,start.year,station,'ATT')
            # 
            #     filename = os.path.join(directory, filename)
            # 
            # 
            # 
            #     # check that the file pattern exists
            #     if not glob.glob(filename):
            #         msg = 'view.py cannot find file: {}'.format(filename)
            #         # print(msg)
            #         logging.info(msg)
            #         continue
            # 
            #     stream = read(filename)
            stream = stream_from_directory_new(
                top_level_dir=top_level_dir,
                start_time=start,
                stations=[station],
                channels=['ATT'],
                merge_locations=False,
                end_time=None)


            # print('temporary measure !!!')    
            # for tr in stream:
            #     print('before')
            #     print(tr.stats)
            #     tr.stats.delta = 0.6036
            #     print(tr.stats)
            #     print('end')      
            

            stream.merge()
            # print(len(stream))
            # 
            # print(stream)
            # exit()

            stream = stream.sort(keys=['starttime'])

            remove_negative_ones(stream)
            

            for i, tr in enumerate(stream):
                # slice the array with a step of the sampling_average
                # e.g. if sampling_average is the same as the
                # samples in a frame (360) then every 360th point
                # this gives a timestamp at the beginning of each frame

                # print(tr.stats)
                # 
                # for a in range(-10,0):
                #     print(UTCDateTime(tr.data[a]))
                # 

                # tr.plot(method='full')
                
                tr.data = tr.data[0::sampling_average]


                # it is only valid to 3 decimal places
                # tr.data = np.ma.round(tr.data,3)
                # diff function to calculate consecutive values
                # note that we now have one value fewer than we did before
                # print(tr.data[-1])

                if np.ma.isMaskedArray(tr.data):

                    # print('yes masked')
                    tr.data = np.ma.diff(tr.data)
                    # tr.plot(method='full')
                else: 
                    tr.data = np.diff(tr.data)


                if tr.stats.npts > 0:
                    tr.data = np.ma.round(tr.data,3)

                    tr.data = tr.data / sampling_average

                    # if calc_moving_average:
                    #     tr.data = moving_average(tr.data)

                    # # if any of the values are zero, make a masked trace
                    # if 0 in tr.data:
                    #     tr.data = ma.masked_equal(tr.data, 0)

                    # correct the starttime - half the sampling_average for the averaging,
                    tr.stats.starttime = tr.stats.starttime + (sampling_average/2)*DELTA*4
                    # print('delta ', tr.stats.delta)
                    # print(DELTA)
                    tr.stats.delta = DELTA*4 * sampling_average
                    tr.stats.mseed.encoding = 'FLOAT64'
                    tr.stats.channel = 'ASI'

                    if tr.stats.npts < 1:
                        stream.remove(tr)
                    elif isinstance(tr.data,np.ma.MaskedArray) and tr.data.count() < 1:
                        stream.remove(tr)

                else:
                    stream.remove(tr)

                for tr in stream:
                    if np.ma.isMaskedArray(tr.data):
                        # mask the ASI trace which has gaps with a fill value
                        tr.data = tr.data.filled(fill_value=INVALID_99)

                # XXXX

                for tr in stream: 
                    # print('what do I have now? - should be one day, one station')
                    # print(stream)
                    # stream.plot()
                    location = tr.stats.location

                    av_filename = find_filename_date_upper(station, location, 'ASI', start)
                    

                    # make the subdirectory with the station name
                    # av_directory = make_dir(av_dir,start.year,station,'ASI')
                    av_directory = make_dir(av_dir,station,start,lower=False)
                    av_filename = os.path.join(av_directory, av_filename)

# /Users/nunn/lunar_data/PDART_ELAPSED_TIME/1973/S12/ASI/S12.XA..ASI.1973.360.gz


                    
                    # save
                    # stream = stream.split()


                    tr.write(av_filename, 'MSEED')
                    print(tr)
                    print(av_filename)
                    


                    av_filename_gzip = '%s.gz' % (av_filename)
                    # # this is slow
                    if write_gzip:
                        with open(av_filename, 'rb') as f_in, gzip.open(av_filename_gzip, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                        os.unlink(av_filename)
                        print('writing file ', av_filename)

                

        # increment the time interval
        start += time_interval
