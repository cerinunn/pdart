#!/usr/bin/env python

'''
Plot an overview of the seismgrams - can choose different channels.

'''

from __future__ import print_function
import numpy as np
import os

import matplotlib  
from matplotlib import pyplot as plt

from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import stream_from_directory_new
from obspy.imaging.util import (_set_xaxis_obspy_dates, _id_key, _timestring)
from matplotlib.dates import date2num
from pdart.util import relative_timing_trace, relative_timing_stream
from matplotlib.dates import AutoDateFormatter, AutoDateLocator, HOURLY
from obspy.core import Stream
import pandas as pd
from pdart.util import linear_interpolation

    # start_time = UTCDateTime('1971-02-07T00:45:00')

SECONDS_PER_DAY = 3600.0 * 24.0
DELTA = 0.1509433962


# From the Tol color palette
# https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
tol_dark_blue ='#341B88' 
tol_green = '#357932'
tol_blue_green = '#4CAA9A'
tol_pale_blue = '#82C0E1'


def plot_timing_multipage(top_level_dir_good=None, start_time=None, end_time=None, stations=['S12','S14','S15','S16'], channels=['MH1','MHZ'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False):

    if end_time is None or start_time is None or top_level_dir_good is None:
        print('start_time, end_time and top_level_dir_good are all required')
        return

    start_time_local = start_time

    while start_time_local < end_time:

        out_filename = 'timing_dir/checkingchannels_{}.png'.format(start_time_local.strftime('%Y.%j'))
        out_filename = os.path.join(out_dir,out_filename)

        fig = plt.figure(figsize=(8.5, 11))
        # gs = gridspec.GridSpec(7, 2, hspace=0.001)
        gs = gridspec.GridSpec(5, 1, hspace=0.1)

        fig.suptitle('Checking Channels', fontsize=16)

        for local_count in range(0,5,1):

            end_time_local = start_time_local+3600*24
        
            # print('Times ', local_count, start_time_local, end_time_local)
            stream = Stream()
            if end_time_local <= end_time:
                stream = stream_from_directory_new(
                  top_level_dir=top_level_dir_good,
                  start_time=start_time_local,
                  stations=stations,
                  channels=channels,
                  end_time=end_time_local)

                # stream=stream.select(channel=channel).copy()
                stream.trim(starttime=start_time,endtime=end_time)

                # print(stream)

                # relative_timing_stream(stream)
                for tr in stream:
                    linear_interpolation(tr,interpolation_limit=1)
                
                xmin = time_to_xvalue(start_time_local)
                xmax = time_to_xvalue(start_time_local + 3600*24)

            # x_annotate = xmax + (xmax-xmin)*.01
            x_annotate = 86400

            #######

            ax0 = plt.subplot(gs[local_count])

            # y_max = 1
            # y_min = -1

            for trace in stream:
                label = trace.stats.channel
                if label == 'MH1':
                    color=tol_dark_blue
                elif label == 'MH2':
                    color=tol_green    
                elif label == 'MHZ':
                    color=tol_blue_green  
                elif label == 'SHZ':
                    color=tol_pale_blue

        #         x1 = trace.times()[0]
        #         x2 = trace.times()[-1]
        #         y1 = trace.data[0]
        #         y2 = trace.data[-1]
        # 
        # #         y = m * x + b
        #         m = (y1 - y2) / (x1 - x2)
        #         b = (x1 * y2 - x2 * y1) / (x1 - x2)
        # 
        #         straight_line = m*trace.times() + b
        #         trace_straight = trace.copy()
        #         trace_straight.data = straight_line
                
                # plot_trace(ax0,trace_straight,type='relative',color='r',linewidth=1)

                plot_trace(ax0,trace,type='relative',color=color,linewidth=2)
                # print(x_annotate,trace.data[-1])
                # ax0.annotate(label, xy=(x_annotate,trace.data[-1]),
                #   fontsize=13, horizontalalignment="left", verticalalignment="center",
                #   xycoords="data", color=color, annotation_clip=False)

                # max1 = trace.data.max()
                # if max1 > y_max:
                #     y_max = max1
                # min1 = trace.data.min()
                # if min1 < y_min:
                #     y_min = min1

            if len(stream) == 0:
                plt.yticks([])  
 
            plt.xticks([])
            ax0.set_xlim(0, 86400)
            ax0.set_ylim(450, 550)
            # print('local count',local_count, y_min, y_max )

            if local_count == 4:
                plt.xticks([0,6*3600,12*3600,18*3600,24*3600],['00:00', '06:00', '12:00', '18:00', '00:00'])
                plt.xlabel(channels)

            # plt.xticks(np.arange(5), ['0', '6', '12', '18', '24'])

            ax0.annotate(start_time_local.strftime('%Y-%m-%d %Y.%-j'), xy=(0.01,0.99),
              fontsize=13, horizontalalignment="left", verticalalignment="top",
              xycoords="axes fraction", color='k', annotation_clip=False)

            start_time_local += 3600*24

    


        plt.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.12)

        # plt.xlabels(['0', '6', '12', '18', '24'])

        if save_fig:
            print('Writing file ', out_filename)
            plt.savefig(out_filename)
        if plot_fig:
            plt.show()
        plt.close()

# def times_to_minutes(times_in_seconds):
#     return ((times_in_seconds / 60) - 2)

def time_to_xvalue(t):
    return date2num(t.datetime)

def plot_trace(ax,trace,type='normal', color='k',linewidth=1,linestyle='solid'):

    if type == 'relative':
        # use seconds of relative sample times and shift by trace's
        # start time, which was set relative to `reftime`.
        reftime = UTCDateTime(year=trace.stats.starttime.year,
          julday=trace.stats.starttime.julday)
        x_values = (
            trace.times() + (trace.stats.starttime - reftime))
    else:
        # convert seconds of relative sample times to days and add
        # start time of trace.
        x_values = ((trace.times() / SECONDS_PER_DAY) +
                    date2num(trace.stats.starttime.datetime))
    ax.plot(x_values, trace.data, color=color,
            linewidth=linewidth, linestyle=linestyle)

def plot_change(ax,trace,idx_list_neg,idx_list_pos,type='normal',neg_color='r',pos_color='g'):
    if len(idx_list_neg) > 0: 
        x_values = []
        y_values = []
        for idx in idx_list_neg:
            if type == 'relative':
                # use seconds of relative sample times and shift by trace's
                # start time, which was set relative to `reftime`.
                reftime = UTCDateTime(year=trace.stats.starttime.year,
                  julday=trace.stats.starttime.julday)
                x_value = (
                    trace.times()[idx] + (trace.stats.starttime - reftime))
                x_values.append(x_value)
                
            else:
                # convert seconds of relative sample times to days and add
                # start time of trace.
                x_value = ((trace.times()[idx] / SECONDS_PER_DAY) +
                            date2num(trace.stats.starttime.datetime))
                x_values.append(x_value)
            y_values.append(trace.data[idx])
            
        # print('x ', x_values)
        # print('y ', y_values)
        ax.scatter(x_values, y_values, color=neg_color, s=100)

    if len(idx_list_pos) > 0: 
        x_values = []
        y_values = []
        for idx in idx_list_pos:
            if type == 'relative':
                # use seconds of relative sample times and shift by trace's
                # start time, which was set relative to `reftime`.
                reftime = UTCDateTime(year=trace.stats.starttime.year,
                  julday=trace.stats.starttime.julday)
                x_value = (
                    trace.times()[idx] + (trace.stats.starttime - reftime))
                x_values.append(x_value)
                
            else:
                # convert seconds of relative sample times to days and add
                # start time of trace.
                x_value = ((trace.times()[idx] / SECONDS_PER_DAY) +
                            date2num(trace.stats.starttime.datetime))
                x_values.append(x_value)
            y_values.append(trace.data[idx])
            
        # print('x ', x_values)
        # print('y ', y_values)
        ax.scatter(x_values, y_values, color=pos_color, s=100)

        # print('end of plot change')


def plot_set_x_ticks(type='normal', *args, **kwargs):  # @UnusedVariable
    """
    Goes through all axes in pyplot and sets time ticks on the x axis.
    """
    from matplotlib import ticker
# axes = plt.gcf().get_axes()
    # plt.gcf().subplots_adjust(hspace=0)
    # Loop over all but last axes.
    axes = plt.gcf().get_axes()
    for ax in axes[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)
    # set bottom most axes:
    ax = axes[-1]
    if type == "relative":
        locator = ticker.MaxNLocator(5)
        ax.xaxis.set_major_locator(locator)
    else:
        _set_xaxis_obspy_dates(ax)
    plt.setp(ax.get_xticklabels(), fontsize='medium',
             rotation=0)

def times_to_seconds(times_in_seconds):
    return (times_in_seconds - 120)

if __name__ == "__main__":


    # view channel MHZ and MH1 from 1976-11-15T00:00:00.000000Z to 1976-11-20T00:00:00.000000Zs
    plot_timing_multipage(top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime('1976-11-15T00:00:00.000000Z'), end_time=UTCDateTime('1976-11-20T00:00:00.000000Z'), stations=['S14'], channels=['MHZ','MH1'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

