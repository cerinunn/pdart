#!/usr/bin/env python

'''
Plot timing divergence (difference between sample time and timestamp)

Used for the timing divergence plot in the paper. 
Extensively used for checking and correcting the timing - (e.g. checking 
for the software clock error and the tapehead error where the 
work tapes were sometimes 400 or 800 ms late.) Zipped versions of these
plots are available (timing_<YEAR>.zip in pdart/Electronic_Supplement/)
'''

from __future__ import print_function
import numpy as np
import os

# Qt5Agg seems to work best on Mac - try 'TkAgg' if that works for you
# put this after the other imports, otherwise it can be overridden
import matplotlib  
# matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import stream_from_directory_new
from obspy.imaging.util import (_set_xaxis_obspy_dates, _id_key, _timestring)
from matplotlib.dates import date2num
from pdart.util import relative_timing_trace, relative_timing_stream, remove_negative_ones
from matplotlib.dates import AutoDateFormatter, AutoDateLocator, HOURLY
from obspy.core import Stream
import pandas as pd

    # start_time = UTCDateTime('1971-02-07T00:45:00')

SECONDS_PER_DAY = 3600.0 * 24.0
DELTA = 0.1509433962
SHORT_DELTA = 0.1509


# From the Tol color palette
# https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
tol_dark_blue ='#341B88' 
tol_green = '#357932'
tol_blue_green = '#4CAA9A'
tol_pale_blue = '#82C0E1'
midnight_blue = '#191970'


def plot_timing_multipage(top_level_dir_bad=None, top_level_dir_good=None, start_time=None, end_time=None, stations=['S12','S14','S15','S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False):

    if end_time is None or start_time is None or top_level_dir_bad is None or top_level_dir_good is None:
        print('start_time, end_time, top_level_dir_bad and top_level_dir_good are all required')
        return

    start_time_local = start_time

    while start_time_local < end_time:


        if len(stations) > 1:
            out_filename = 'timing_dir/timing_divergence_{}_all_XX.png'.format(start_time_local.strftime('%Y.%j'))
        else:
            out_filename = 'timing_dir/timing_divergence_{}_{}XX.png'.format(start_time_local.strftime('%Y.%j'),stations[0])
        out_filename = os.path.join(out_dir,out_filename)

        fig = plt.figure(figsize=(8.5, 11))
        # gs = gridspec.GridSpec(7, 2, hspace=0.001)
        gs = gridspec.GridSpec(7, 2, hspace=0.1)

        fig.suptitle('Divergence between Sample Time and Timestamp', fontsize=16)

        for local_count in range(0,10,2):

            end_time_local = start_time_local+3600*24

                
            # print('Times ', local_count, start_time_local, end_time_local)
            stream = Stream()
            if end_time_local <= end_time:
                stream = stream_from_directory_new(
                  top_level_dir=top_level_dir_bad,
                  start_time=start_time_local,
                  stations=stations,
                  channels=['ATT'],
                  end_time=end_time_local)

                stream=stream.select(channel='ATT').copy()
                stream.trim(starttime=start_time,endtime=end_time)

                # remove_negative_ones(stream,channels=['ATT'])

                # print(stream)

                relative_timing_stream(stream)
                
                xmin = time_to_xvalue(start_time_local)
                xmax = time_to_xvalue(start_time_local + 3600*24)

            # x_annotate = xmax + (xmax-xmin)*.01
            x_annotate = 86400

            #######

            ax0 = plt.subplot(gs[local_count])

            y_max = 1
            y_min = -1

            for trace in stream:
                label = trace.stats.station
                if label == 'S12':
                    color=tol_dark_blue
                elif label == 'S14':
                    color=tol_green    
                elif label == 'S15':
                    color=tol_blue_green  
                elif label == 'S16':
                    color=tol_pale_blue
                elif label == 'S11':
                    color=midnight_blue
                

                x1 = trace.times()[0]
                x2 = trace.times()[-1]
                y1 = trace.data[0]
                y2 = trace.data[-1]
                
        #         y = m * x + b
                m = (y1 - y2) / (x1 - x2)
                b = (x1 * y2 - x2 * y1) / (x1 - x2)

                straight_line = m*trace.times() + b
                trace_straight = trace.copy()
                trace_straight.data = straight_line
                
                plot_trace(ax0,trace_straight,type='relative',color='r',linewidth=1)

                plot_trace(ax0,trace,type='relative',color=color,linewidth=2)
                # print(x_annotate,trace.data[-1])
                ax0.annotate(label, xy=(x_annotate,trace.data[-1]),
                  fontsize=13, horizontalalignment="left", verticalalignment="center",
                  xycoords="data", color=color, annotation_clip=False)

                max1 = trace.data.max()
                if max1 > y_max:
                    y_max = max1
                min1 = trace.data.min()
                if min1 < y_min:
                    y_min = min1

            if len(stream) == 0:
                plt.yticks([])  
 
            plt.xticks([])
            ax0.set_xlim(0, 86400)
            ax0.set_ylim(y_min, y_max)
            # print('local count',local_count, y_min, y_max )

            if local_count == 8:
                plt.xticks([0,6*3600,12*3600,18*3600,24*3600],['00:00', '06:00', '12:00', '18:00', '00:00'])
                plt.xlabel('ATT relative (before)')

            # plt.xticks(np.arange(5), ['0', '6', '12', '18', '24'])

            ax0.annotate(start_time_local.strftime('%Y-%m-%d %Y.%-j'), xy=(0.01,0.99),
              fontsize=13, horizontalalignment="left", verticalalignment="top",
              xycoords="axes fraction", color='k', annotation_clip=False)

            # print(local_count)


            ####### Get the traces for the right-hand side of the plot

            # get it again - will get the updated version if necessary

            stream = Stream()
            if end_time_local <= end_time:
                stream = stream_from_directory_new(
                  top_level_dir=top_level_dir_good,
                  start_time=start_time_local,
                  stations=stations,
                  channels=['ATT'],
                  end_time=end_time_local)

                stream=stream.select(channel='ATT').copy()
                stream.trim(starttime=start_time,endtime=end_time)

                # print(stream)

                relative_timing_stream(stream)
                
                xmin = time_to_xvalue(start_time_local)
                xmax = time_to_xvalue(start_time_local + 3600*24)

            # x_annotate = xmax + (xmax-xmin)*.01
            x_annotate = 86400

            # plot on the right
            ax0 = plt.subplot(gs[local_count+1])

            y_max = 1
            y_min = 0

            for trace in stream:
                label = trace.stats.station
                if label == 'S12':
                    color=tol_dark_blue
                elif label == 'S14':
                    color=tol_green    
                elif label == 'S15':
                    color=tol_blue_green  
                elif label == 'S16':
                    color=tol_pale_blue
                elif label == 'S11':
                    color=midnight_blue

                x1 = trace.times()[0]
                x2 = trace.times()[-1]
                y1 = trace.data[0]
                y2 = trace.data[-1]
                
        #         y = m * x + b
                m = (y1 - y2) / (x1 - x2)
                b = (x1 * y2 - x2 * y1) / (x1 - x2)

                straight_line = m*trace.times() + b
                trace_straight = trace.copy()
                trace_straight.data = straight_line

                plot_trace(ax0,trace_straight,type='relative',color='r',linewidth=1)

                plot_trace(ax0,trace,type='relative',color=color,linewidth=2)
                # print(x_annotate,trace.data[-1])
                ax0.annotate(label, xy=(x_annotate,trace.data[-1]),
                  fontsize=13, horizontalalignment="left", verticalalignment="center",
                  xycoords="data", color=color, annotation_clip=False)

                max1 = trace.data.max()
                if max1 > y_max:
                    y_max = max1
                min1 = trace.data.min()
                if min1 < y_min:
                    y_min = min1
    
            if len(stream) == 0:
                plt.yticks([])  
                
            plt.xticks([])
            ax0.set_xlim(0, 86400)
            ax0.set_ylim(y_min-0.75, y_max+0.75)
            # print('local count',local_count+1, y_min, y_max )

            # plt.xticks(np.arange(5), ['0', '6', '12', '18', '24'])

            ax0.annotate(start_time_local.strftime('%Y-%m-%d %Y.%-j'), xy=(0.01,0.99),
              fontsize=13, horizontalalignment="left", verticalalignment="top",
              xycoords="axes fraction", color='k', annotation_clip=False)


            if local_count == 8:
                plt.xticks([0,6*3600,12*3600,18*3600,24*3600],['00:00', '06:00', '12:00', '18:00', '00:00'])
                plt.xlabel('ATT relative (after)')

            start_time_local += 3600*24

    
        print('at the end')

        plt.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.12)

        # plt.xlabels(['0', '6', '12', '18', '24'])

        if save_fig:
            print('Writing file ', out_filename)
            plt.savefig(out_filename)
        if plot_fig:
            plt.show()
        plt.close()


def plot_timing(stream=None, start_time=None, timedelta_hours=24, end_time=None, stations=['S12','S14','S15','S16'], include_line=False, out_dir='../extra_plots_output', out_filename=None, save_fig=True, plot_fig=True):

    if end_time is None:
        end_time = start_time + timedelta(hours=timedelta_hours)

    

    stream=stream.select(channel='ATT').copy()

    # print('temp')
    # stream.trim(starttime=UTCDateTime('1976-03-01T09:38:34'),endtime=UTCDateTime('1976-03-01T09:39'))
    # fig = stream[0].plot(handle=True, show=False,method='full')
    # plt.ylim(-2,2)
    # plt.show()

    # 194521114.185 194521114.789 194521115.393 -- -- -- -- -- 97260559.3095


    # remove_negative_ones(stream,channels=['ATT'])
    # 
    # stream.trim(starttime=UTCDateTime('1976-03-01T09:38'),endtime=UTCDateTime('1976-03-01T09:39'))
    # fig = stream[0].plot(handle=True, show=False,method='full')
    # plt.ylim(-2,2)
    # plt.show()

    # TODO - fix the end time
    stream.trim(starttime=start_time,endtime=end_time)

    # XXXX
    
    xmin = time_to_xvalue(start_time)
    xmax = time_to_xvalue(end_time)

    x_annotate = xmax + (xmax-xmin)*.01

    ########

    fig = plt.figure(figsize=(6.5, 3))
    gs = gridspec.GridSpec(1, 1, hspace=0.001)

    ax0 = plt.subplot(gs[0])

    for trace in stream:
        label = trace.stats.station
        if label == 'S12':
            color=tol_dark_blue
        elif label == 'S14':
            color=tol_green    
        elif label == 'S15':
            color=tol_blue_green  
        elif label == 'S16':
            color=tol_pale_blue
        elif label == 'S11':
            color=midnight_blue

        trace_divergence = trace.copy() 
        relative_timing_trace(trace_divergence)

        # plot the line based on the divergence trace 
        if include_line:
            x1 = trace_divergence.times()[0]
            x2 = trace_divergence.times()[-1]
            y1 = trace_divergence.data[0]
            y2 = trace_divergence.data[-1]

    #         y = m * x + b
            m = (y1 - y2) / (x1 - x2)
            b = (x1 * y2 - x2 * y1) / (x1 - x2)

            straight_line = m*trace_divergence.times() + b
            trace_straight = trace_divergence.copy()
            trace_straight.data = straight_line

            # print(trace_straight.data.min())
            # print(trace_straight.data.max())
            # print(trace_straight.times()[0])    
            # print(trace_straight.times()[-1])  

            plot_trace(ax0,trace_straight,color='r',linewidth=1) 

        # plot the relative trace
        plot_trace(ax0,trace_divergence,color=color,linewidth=3)

        ax0.annotate(label, xy=(x_annotate,trace_divergence.data[-1]),
          fontsize=13, horizontalalignment="left", verticalalignment="center",
          xycoords="data", color=color, annotation_clip=False)

        ax0.annotate(start_time.strftime('%Y-%m-%d %Y.%-j'), xy=(0.01,0.99),
          fontsize=13, horizontalalignment="left", verticalalignment="top",
          xycoords="axes fraction", color='k', annotation_clip=False)


  

            # # find the shifts based on the divergence trace 
            # idx_list_neg, idx_list_pos = correct_shifts(trace_divergence,trace)
            # 
            # # plot the shifts based on the divergent trace
            # plot_change(ax0,trace_divergence,idx_list_neg, idx_list_pos)

            # 
            # XXXX
        

                
          

    ax0.set_title('Divergence between Sample Time and Timestamp', fontsize=16)
    # ax0.set_ylim(-1,5.5)
    ax0.set_ylabel('Sample Time\nminus Timestamp [s]', rotation=90,fontsize=14, labelpad=-2)

    plt.subplots_adjust(left=0.11, right=0.93, top=0.9, bottom=0.1)

    plot_set_x_ticks()
    ax0.set_xlim(xmin, xmax)


    if out_filename is None: 
        out_filename = 'timing_divergence_{}XX.png'.format(start_time.strftime('%Y.%j'))

    out_filename = os.path.join(out_dir,out_filename)


    if save_fig:
        plt.savefig(out_filename)
    if plot_fig:
        plt.show()
    plt.close()


def plot_timing_from_dir(top_level_dir=None, start_time=None, timedelta_hours=24, end_time=None, stations=['S12','S14','S15','S16'], include_line=False, out_dir='../extra_plots_output', out_filename=None, save_fig=True, plot_fig=True ):

    if end_time is None:
        end_time = start_time + timedelta(hours=timedelta_hours)

    # print('Currently from tmp directory'

    stream = stream_from_directory_new(
      top_level_dir=top_level_dir,
      start_time=start_time,
      stations=stations,
      channels=['ATT'],
      end_time=end_time)

    # print(stream)
    # print(stream[0].stats.starttime)
    # print(UTCDateTime(stream[0].data[0]))
    # exit()

    plot_timing(stream=stream, start_time=start_time, timedelta_hours=timedelta_hours, end_time=end_time, stations=stations, include_line=include_line, out_dir=out_dir, out_filename=out_filename, save_fig=save_fig, plot_fig=plot_fig)

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

    # plot for the paper
    plot_timing_from_dir(top_level_dir='/Users/cnunn/lunar_data/PDART_V2',start_time=UTCDateTime('1973-06-30T00:00:00.00000Z'))

    exit()

    # plot for the Electronic Supplement - all days!!!

    # 1969 - these look really bad!!
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1969,julday=201), end_time=UTCDateTime(year=1969,julday=241), stations=['S11'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    # 1969 - S12 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1969,julday=321), end_time=UTCDateTime(year=1970,julday=1), stations=['S12'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    # 1970 - S12 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1970,julday=1), end_time=UTCDateTime(year=1971,julday=1), stations=['S12'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    # 1971 - S12 / S14 / S15
    # plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1971,julday=1), end_time=UTCDateTime(year=1972,julday=1), stations=['S12','S14','S15'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False)
 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1971,julday=1), end_time=UTCDateTime(year=1972,julday=1), stations=['S12'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1971,julday=1), end_time=UTCDateTime(year=1972,julday=1), stations=['S14'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1971,julday=1), end_time=UTCDateTime(year=1972,julday=1), stations=['S15'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 


    # 1972 - all
    # plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1972,julday=1), end_time=UTCDateTime(year=1973,julday=1), stations=['S12','S14','S15','S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1972,julday=1), end_time=UTCDateTime(year=1973,julday=1), stations=['S12'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1972,julday=1), end_time=UTCDateTime(year=1973,julday=1), stations=['S14'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1972,julday=1), end_time=UTCDateTime(year=1973,julday=1), stations=['S15'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1972,julday=1), end_time=UTCDateTime(year=1973,julday=1), stations=['S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 


    # 1973 - all
    # plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1973,julday=1), end_time=UTCDateTime(year=1974,julday=1), stations=['S12','S14','S15','S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1973,julday=1), end_time=UTCDateTime(year=1974,julday=1), stations=['S12'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1973,julday=1), end_time=UTCDateTime(year=1974,julday=1), stations=['S14'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1973,julday=1), end_time=UTCDateTime(year=1974,julday=1), stations=['S15'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1973,julday=1), end_time=UTCDateTime(year=1974,julday=1), stations=['S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    # 1974
    # plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1974,julday=1), end_time=UTCDateTime(year=1975,julday=1), stations=['S12','S14','S15','S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1974,julday=1), end_time=UTCDateTime(year=1975,julday=1), stations=['S12'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1974,julday=1), end_time=UTCDateTime(year=1975,julday=1), stations=['S14'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1974,julday=1), end_time=UTCDateTime(year=1975,julday=1), stations=['S15'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1974,julday=1), end_time=UTCDateTime(year=1975,julday=1), stations=['S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 


    # 1975
    # plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1975,julday=1), end_time=UTCDateTime(year=1976,julday=1), stations=['S12','S14','S15','S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1975,julday=1), end_time=UTCDateTime(year=1976,julday=1), stations=['S12'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1975,julday=1), end_time=UTCDateTime(year=1976,julday=1), stations=['S14'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1975,julday=1), end_time=UTCDateTime(year=1976,julday=1), stations=['S15'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1975,julday=1), end_time=UTCDateTime(year=1976,julday=1), stations=['S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    # 1976
    # plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1976,julday=1), end_time=UTCDateTime(year=1977,julday=1), stations=['S12','S14','S15','S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1976,julday=1), end_time=UTCDateTime(year=1977,julday=1), stations=['S12'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1976,julday=1), end_time=UTCDateTime(year=1977,julday=1), stations=['S14'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1976,julday=1), end_time=UTCDateTime(year=1977,julday=1), stations=['S15'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1976,julday=1), end_time=UTCDateTime(year=1977,julday=1), stations=['S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    # 1977
    # plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1977,julday=1), end_time=UTCDateTime(year=1977,julday=276), stations=['S12','S14','S15','S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1977,julday=1), end_time=UTCDateTime(year=1977,julday=276), stations=['S12'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1977,julday=1), end_time=UTCDateTime(year=1977,julday=276), stations=['S14'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1977,julday=1), end_time=UTCDateTime(year=1977,julday=276), stations=['S15'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 
    plot_timing_multipage(top_level_dir_bad='/Users/cnunn/lunar_data/PDART', top_level_dir_good='/Users/cnunn/lunar_data/PDART_V2', start_time=UTCDateTime(year=1977,julday=1), end_time=UTCDateTime(year=1977,julday=276), stations=['S16'], out_dir='../extra_plots_output', save_fig=True, plot_fig=False) 

