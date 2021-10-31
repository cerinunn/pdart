#!/usr/bin/env python

from __future__ import print_function
import numpy as np
# Qt5Agg seems to work best on Mac - try 'TkAgg' if that works for you
# put this after the other imports, otherwise it can be overridden
import matplotlib  
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import stream_from_directory_new
from pdart.snippet import remove_negative_ones, linear_interpolation



    # start_time = UTCDateTime('1971-02-07T00:45:00')

# remove_negative_ones(stream_after)
#     for tr in stream_after:
#         linear_interpolation(tr)

def five_tracks():

    onset = UTCDateTime('1971-02-07T00:45:45.600000Z')
    start_time = onset - timedelta(minutes=2)

    stations = ['S14']
    channels = ['MH1', 'MH2', 'MHZ', 'SHZ', 'ATT']
    end_time = start_time + timedelta(minutes=(2+12))

    stream = stream_from_directory_new(
      top_level_dir='/Users/cnunn/lunar_data/PDART_CONTINUOUS_MAIN_TAPES',
      start_time=start_time,
      stations=stations,
      channels=channels,
      end_time=end_time)

    # print(stream)
    print('Applying a linear interpolation across any single sample gaps')
    for tr in stream:
        linear_interpolation(tr)

    # stream.plot(type='relative')
    # exit()


    ########

    fig = plt.figure(figsize=(15, 4))
    gs = gridspec.GridSpec(5, 1, hspace=0.001)

    ax0 = plt.subplot(gs[0])
# 
# 510
# 
# 509
# 
# 512
# 
# 531

    trace_MHZ = stream.select(channel='MHZ')[0]
    ax0.plot(times_to_seconds(trace_MHZ.times()), trace_MHZ.data, color='k')
    ax0.set_xlim(-2*60, 12*60)
    ax0.set_xticks([0,180,360,540,720])
    ax0.set_xticks(np.arange(-120,720,60), minor=True)
    ax0.set_title('5 Tracks', fontsize=20)

    ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)
    ax0.set_ylabel('DU', fontsize=14)
    ax0.annotate(xy=(0.01,0.9), text=onset.strftime("%Y-%m-%d %H:%M:%S"),
      fontsize=13, horizontalalignment="left", verticalalignment="top",
      xycoords="axes fraction")

    ax1 = plt.subplot(gs[1], sharex=ax0)
    trace_MH1 = stream.select(channel='MH1')[0]
    ax1.plot(times_to_seconds(trace_MH1.times()), trace_MH1.data, color='k')

    yticks = ax1.set_yticks([430,510,590])
    # yticks[0].label1.set_visible(False)
    # yticks[-1].label1.set_visible(False)
    # ax1.set_ylim(510-140, 510+140)
    # ax1.set_yticks(np.arange(400, 1000, 200))
    ax1.set_ylabel('DU', fontsize=14)

    ax2 = plt.subplot(gs[2], sharex=ax0)
    trace_MH2 = stream.select(channel='MH2')[0]
    ax2.plot(times_to_seconds(trace_MH2.times()), trace_MH2.data, color='k')

    # ax2.set_ylim(510-60, 510+60)
    ax2.set_ylim(280, 1000)
    ax2.set_yticks([310,510,710])
    ax2.set_ylabel('DU', fontsize=14)

    ax3 = plt.subplot(gs[3], sharex=ax0)
    trace_SHZ = stream.select(channel='SHZ')[0]
    ax3.plot(times_to_seconds(trace_SHZ.times()), trace_SHZ.data, color='k')

    ax3.set_yticks([330,530,730])
    ax3.set_ylim(530-300,530+300)
    ax3.set_ylabel('DU', fontsize=14)

    ax4 = plt.subplot(gs[4], sharex=ax0)
    trace_ATT = stream.select(channel='ATT')[0]
    ax4.plot(times_to_seconds(trace_ATT.times()),
      times_to_seconds(trace_ATT.data-trace_ATT.stats.starttime.timestamp), color='k')

    ax4.set_ylim(-180, 820)
    ax4.set_yticks([0,360,720])
    ax4.set_ylabel('s',rotation=90, fontsize=14)

    ax4.annotate(xy=(0.01,0.9), text='+{} s'.format(onset.timestamp),
      fontsize=13, horizontalalignment="left", verticalalignment="top",
      xycoords="axes fraction")
    ax4.annotate(xy=(0.99,0.1), text=trace_ATT.stats.station,
      fontsize=16, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")
    ax4.set_xlabel('Seconds after arrival time', fontsize=14)

    xticklabels = (ax0.get_xticklabels() + ax1.get_xticklabels() +
      ax2.get_xticklabels() + ax3.get_xticklabels())
    plt.setp(xticklabels, visible=False)

    ax0.tick_params(length=6, width=1, which='minor')
    ax1.tick_params(length=6, width=1, which='minor')
    ax2.tick_params(length=6, width=1, which='minor')
    ax3.tick_params(length=6, width=1, which='minor')
    ax4.tick_params(length=6, width=1, which='minor')


    ax0.yaxis.set_label_coords(-0.04, 0.5)
    ax1.yaxis.set_label_coords(-0.04, 0.5)
    ax2.yaxis.set_label_coords(-0.04, 0.5)
    ax3.yaxis.set_label_coords(-0.04, 0.5)
    ax4.yaxis.set_label_coords(-0.04, 0.5)

    ax0.annotate(xy=(1.01,0.5), text='MHZ', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax1.annotate(xy=(1.01,0.5), text='MH1', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax2.annotate(xy=(1.01,0.5), text='MH2', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax3.annotate(xy=(1.01,0.5), text='SHZ', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax4.annotate(xy=(1.01,0.5), text='ATT', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    plt.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.12)
    plt.savefig('5tracks_incl_SHZ_XXXX.png')
    plt.show()


# updated 25-08-21
def five_tracks_original():

    onset = UTCDateTime('1971-02-07T00:45:46.800000Z')
    start_time = onset - timedelta(minutes=2)

    stations = ['S12']
    channels = ['MH1', 'MH2', 'MHZ', 'ATT', 'AFR']
    end_time = start_time + timedelta(minutes=15)

    stream = stream_from_directory_new(
      top_level_dir='/Users/cnunn/lunar_data/PDART_CONTINUOUS_MAIN_TAPES',
      start_time=start_time,
      stations=stations,
      channels=channels,
      end_time=end_time)

    ########

    fig = plt.figure(figsize=(15, 4))
    gs = gridspec.GridSpec(5, 1, hspace=0.001)

    ax0 = plt.subplot(gs[0])

    trace_MHZ = stream.select(channel='MHZ')[0]
    ax0.plot(times_to_seconds(trace_MHZ.times()), trace_MHZ.data, color='k')
    ax0.set_xlim(-2*60, 15*60)
    ax0.set_xticks(np.arange(0,15*60,6*50))
    ax0.set_xticks(np.arange(-180,15*60,60), minor=True)
    ax0.set_title('5 Tracks', fontsize=20)

    ax0.set_yticks(np.arange(480, 560, 20))
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    ax0.set_ylim(510-50, 510+50)
    ax0.set_ylabel('DU', fontsize=14)
    ax0.annotate(xy=(0.01,0.9), text=onset.strftime("%Y-%m-%d %H:%M:%S"),
      fontsize=13, horizontalalignment="left", verticalalignment="top",
      xycoords="axes fraction")

    ax1 = plt.subplot(gs[1], sharex=ax0)
    trace_MH1 = stream.select(channel='MH1')[0]
    ax1.plot(times_to_seconds(trace_MH1.times()), trace_MH1.data, color='k')

    # yticks = ax1.set_yticks(np.arange(400, 1100, 200))
    # yticks[0].label1.set_visible(False)
    # yticks[-1].label1.set_visible(False)
    ax1.set_ylim((200, 1000))
    ax1.set_yticks(np.arange(400, 1000, 200))
    ax1.set_ylabel('DU', fontsize=14)

    ax2 = plt.subplot(gs[2], sharex=ax0)
    trace_MH2 = stream.select(channel='MH2')[0]
    ax2.plot(times_to_seconds(trace_MH2.times()), trace_MH2.data, color='k')

    ax2.set_ylim((380, 580))
    ax2.set_yticks(np.arange(430, 580, 50))
    ax2.set_ylabel('DU', fontsize=14)

    ax3 = plt.subplot(gs[3], sharex=ax0)
    trace_AFR = stream.select(channel='AFR')[0]
    ax3.plot(times_to_seconds(trace_AFR.times()), trace_AFR.data, color='k')

    ax3.set_yticks(np.arange(0, 100, 45))
    ax3.set_ylim((-10, 100))
    ax3.set_ylabel('Frame-\ncount',rotation=90, fontsize=14)

    ax4 = plt.subplot(gs[4], sharex=ax0)
    trace_ATT = stream.select(channel='ATT')[0]
    ax4.plot(times_to_seconds(trace_ATT.times()),
      times_to_seconds(trace_ATT.data-trace_ATT.stats.starttime.timestamp), color='k')

    ax4.set_ylim((-200, 900))
    ax4.set_yticks(np.arange(-200, 800, 200))
    ax4.set_ylabel('s',rotation=90, fontsize=14)

    ax4.annotate(xy=(0.01,0.9), text='+{} s'.format(onset.timestamp),
      fontsize=13, horizontalalignment="left", verticalalignment="top",
      xycoords="axes fraction")
    ax4.set_xlabel('Seconds after arrival time', fontsize=14)

    xticklabels = (ax0.get_xticklabels() + ax1.get_xticklabels() +
      ax2.get_xticklabels() + ax3.get_xticklabels())
    plt.setp(xticklabels, visible=False)

    ax0.tick_params(length=6, width=1, which='minor')
    ax1.tick_params(length=6, width=1, which='minor')
    ax2.tick_params(length=6, width=1, which='minor')
    ax3.tick_params(length=6, width=1, which='minor')
    ax4.tick_params(length=6, width=1, which='minor')


    ax0.yaxis.set_label_coords(-0.04, 0.5)
    ax1.yaxis.set_label_coords(-0.04, 0.5)
    ax2.yaxis.set_label_coords(-0.04, 0.5)
    ax3.yaxis.set_label_coords(-0.03, 0.5)
    ax4.yaxis.set_label_coords(-0.04, 0.5)

    ax0.annotate(xy=(1.01,0.5), text='MHZ', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax1.annotate(xy=(1.01,0.5), text='MH1', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax2.annotate(xy=(1.01,0.5), text='MH2', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax3.annotate(xy=(1.01,0.5), text='AFR', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax4.annotate(xy=(1.01,0.5), text='ATT', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    plt.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.12)
    plt.savefig('../extra_plots_output/5tracks_XXXX.png')
    plt.show()

# def times_to_minutes(times_in_seconds):
#     return ((times_in_seconds / 60) - 2)

def times_to_seconds(times_in_seconds):
    return (times_in_seconds - 120)

if __name__ == "__main__":
    five_tracks()
