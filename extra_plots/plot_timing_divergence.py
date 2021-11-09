#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import os

# Qt5Agg seems to work best on Mac - try 'TkAgg' if that works for you
# put this after the other imports, otherwise it can be overridden
import matplotlib  
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import stream_from_directory_new
from obspy.imaging.util import (_set_xaxis_obspy_dates, _id_key, _timestring)
from matplotlib.dates import date2num
from pdart.util import relative_timing_trace

    # start_time = UTCDateTime('1971-02-07T00:45:00')

SECONDS_PER_DAY = 3600.0 * 24.0
DELTA = 0.1509433962


def plot_timing(stream=None, start_time=None, timedelta_hours=24, end_time=None, stations=['S12','S14','S15','S16'], out_dir='../extra_plots_output', out_filename=None, save_fig=True, plot_fig=True):

    if end_time is None:
        end_time = start_time + timedelta(hours=timedelta_hours)

    stream=stream.select(channel='ATT').copy()
    stream.trim(starttime=start_time,endtime=end_time)

    relative_timing_trace(stream)
    
    xmin = time_to_xvalue(start_time)
    xmax = time_to_xvalue(end_time)

    x_annotate = xmax + (xmax-xmin)*.01

    ########

    fig = plt.figure(figsize=(15, 4))
    gs = gridspec.GridSpec(1, 1, hspace=0.001)

    ax0 = plt.subplot(gs[0])

    # From the Tol color palette
    # https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
    tol_dark_blue ='#341B88' 
    tol_green = '#357932'
    tol_blue_green = '#4CAA9A'
    tol_pale_blue = '#82C0E1'

    stream_S12 = stream.select(station='S12')
    if len(stream_S12) > 0:
        trace_S12 = stream_S12[0]
    
        plot_trace(ax0,trace_S12,color=tol_dark_blue,linewidth=3)
        ax0.annotate('S12', xy=(x_annotate,trace_S12.data[-1]),
          fontsize=13, horizontalalignment="left", verticalalignment="center",
          xycoords="data", color=tol_dark_blue, annotation_clip=False)

        ax0.set_title('Divergence between Sample Time and Timestamp', fontsize=16)
        # ax0.set_ylim(-1,5.5)
        ax0.set_ylabel('Sample Time\nminus Timestamp [s]', rotation=90,fontsize=14)

    stream_S14 = stream.select(station='S14')
    if len(stream_S14) > 0:
        trace_S14 = stream_S14[0]

        plot_trace(ax0,trace_S14,color=tol_green,linewidth=3)
        ax0.annotate('S14', xy=(x_annotate,trace_S14.data[-1]),
          fontsize=13, horizontalalignment="left", verticalalignment="center", 
          xycoords="data", color=tol_green, annotation_clip=False)

    stream_S15 = stream.select(station='S15')
    if len(stream_S15) > 0:
        trace_S15 = stream_S15[0]

        plot_trace(ax0,trace_S15,color=tol_blue_green,linewidth=3)
        ax0.annotate('S15', xy=(x_annotate,trace_S15.data[-1]),
          fontsize=13, horizontalalignment="left", verticalalignment="center",
          xycoords="data", color=tol_blue_green, annotation_clip=False)

    stream_S16 = stream.select(station='S16')
    if len(stream_S16) > 0:
        trace_S16 = stream_S16[0]

        plot_trace(ax0,trace_S16,color=tol_pale_blue,linewidth=3)
        ax0.annotate('S16', xy=(x_annotate,trace_S16.data[-1]),
          fontsize=13, horizontalalignment="left", verticalalignment="center",
          xycoords="data", color=tol_pale_blue, annotation_clip=False)

    plt.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.12)

    plot_set_x_ticks()
    ax0.set_xlim(xmin, xmax)


    if out_filename is None: 
        out_filename = 'timing_divergence_{}XX.png'.format(start_time.strftime('%Y.%j'))

    out_filename = os.path.join(out_dir,out_filename)


    if save_fig:
        plt.savefig(out_filename)
    if plot_fig:
        plt.show()


def plot_timing_from_dir(top_level_dir=None, start_time=None, timedelta_hours=24, end_time=None, stations=['S12','S14','S15','S16'], out_dir='../extra_plots_output', out_filename=None, save_fig=True, plot_fig=True ):

    if end_time is None:
        end_time = start_time + timedelta(hours=timedelta_hours)

    # print('Currently from tmp directory'

    stream = stream_from_directory_new(
      top_level_dir=top_level_dir,
      start_time=start_time,
      stations=stations,
      channels=['ATT'],
      end_time=end_time)

    plot_timing(stream=stream, start_time=start_time, timedelta_hours=timedelta_hours, end_time=end_time, stations=stations, out_dir=out_dir, out_filename=out_filename, save_fig=save_fig, plot_fig=plot_fig)

# def times_to_minutes(times_in_seconds):
#     return ((times_in_seconds / 60) - 2)

def time_to_xvalue(t):
    return date2num(t.datetime)

def plot_trace(ax,trace,type='normal', color='k',linewidth=1,linestyle='solid'):

    if type == 'relative':
        # use seconds of relative sample times and shift by trace's
        # start time, which was set relative to `reftime`.
        x_values = (
            trace.times() + (trace.stats.starttime - self.reftime))
    else:
        # convert seconds of relative sample times to days and add
        # start time of trace.
        x_values = ((trace.times() / SECONDS_PER_DAY) +
                    date2num(trace.stats.starttime.datetime))
    ax.plot(x_values, trace.data, color=color,
            linewidth=linewidth, linestyle=linestyle)


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
    plot_timing_from_dir(top_level_dir='/Users/cnunn/lunar_data/PDART',start_time=UTCDateTime('1973-06-30T00:00:00.00000Z'))

