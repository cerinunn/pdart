#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import stream_from_directory_new

    # start_time = UTCDateTime('1971-02-07T00:45:00')

# 1970-03-24T15:59:28.208000Z

def calibration_example():

    onset = UTCDateTime('1970-03-26T02:04:48.295000Z')
    start_time = onset - timedelta(minutes=1)

    stations = ['S12']
    channels = ['MH1', 'MH2', 'MHZ']
    end_time = start_time + timedelta(minutes=4.5)

    stream = stream_from_directory_new(
      top_level_dir='/Users/cnunn/lunar_data/PDART_CORRECTED',
      start_time=start_time,
      stations=stations,
      channels=channels,
      end_time=end_time)

    onset2 = UTCDateTime('1976-05-18T04:29:04.511000Z')
    start_time = onset2 - timedelta(minutes=1)

    stations = ['S12']
    channels = ['MH1', 'MH2', 'MHZ']
    end_time = start_time + timedelta(minutes=9)

    stream2 = stream_from_directory_new(
      top_level_dir='/Users/cnunn/lunar_data/PDART_CORRECTED',
      start_time=start_time,
      stations=stations,
      channels=channels,
      end_time=end_time)

    # stream2.plot()
    # exit()

    ########

    fig = plt.figure(figsize=(7.5, 4.5))
    gs = gridspec.GridSpec(7, 1, hspace=0.001)

    ax0 = plt.subplot(gs[0])

    trace_MHZ = stream.select(channel='MHZ')[0]
    ax0.plot(times_to_minutes(trace_MHZ.times()), trace_MHZ.data, color='k')
    ax0.set_xticks(np.arange(0,6,1))
    ax0.set_xticks(np.arange(-1,6,1/4), minor=True)
    ax0.set_title('Example Calibration Pulses', fontsize=20)


    ax0.set_yticks(np.arange(510,540,10))
    ax0.set_ylim(500, 540)
    # ax0.annotate(xy=(0.01,0.9), s=onset.strftime("%Y-%m-%d %H:%M:%S"),
    #   fontsize=13, horizontalalignment="left", verticalalignment="top",
    #   xycoords="axes fraction")

    ax1 = plt.subplot(gs[1], sharex=ax0)
    trace_MH1 = stream.select(channel='MH1')[0]
    ax1.plot(times_to_minutes(trace_MH1.times()), trace_MH1.data, color='k')

    # yticks = ax1.set_yticks(np.arange(400, 1100, 200))
    # yticks[0].label1.set_visible(False)
    # yticks[-1].label1.set_visible(False)
    ax1.set_ylim((500, 540))
    ax1.set_yticks(np.arange(510,540,10))
    ax1.set_ylabel('DU', fontsize=14, rotation=0)

    ax2 = plt.subplot(gs[2], sharex=ax0)
    trace_MH2 = stream.select(channel='MH2')[0]
    ax2.plot(times_to_minutes(trace_MH2.times()), trace_MH2.data, color='k')

    ax2.set_ylim(((500, 540)))
    ax2.set_yticks(np.arange(510,540,10))
    ax2.annotate(xy=(0.99,0.4), s=onset.strftime("%Y-%m-%d %H:%M:%S"),
      fontsize=13, horizontalalignment="right", verticalalignment="top",
      xycoords="axes fraction")


    ax3 = plt.subplot(gs[4])
    trace_MHZ = stream2.select(channel='MHZ')[0]
    ax3.plot(times_to_minutes(trace_MHZ.times()), trace_MHZ.data, color='k')
    ax3.set_ylim(((200, 890)))
    ax3.set_yticks(np.arange(300,900,200))

    ax4 = plt.subplot(gs[5], sharex=ax3)
    trace_MH1 = stream2.select(channel='MH1')[0]
    ax4.plot(times_to_minutes(trace_MH1.times()), trace_MH1.data, color='k')
    ax4.set_ylim(((200, 890)))
    ax4.set_yticks(np.arange(300,900,200))
    ax4.set_ylabel('DU', fontsize=14, rotation=0)

    ax5 = plt.subplot(gs[6], sharex=ax3)
    trace_MH2 = stream2.select(channel='MH2')[0]
    ax5.plot(times_to_minutes(trace_MH2.times()), trace_MH2.data, color='k')
    ax5.set_ylim(((200, 890)))
    ax5.set_yticks(np.arange(300,900,200))

    ax5.annotate(xy=(0.99,0.4), s=onset2.strftime("%Y-%m-%d %H:%M:%S"),
      fontsize=13, horizontalalignment="right", verticalalignment="top",
      xycoords="axes fraction")



    ax3.set_xticks(np.arange(-1,10,1/4), minor=True)
    ax3.set_xlim(-1,8)
    ax5.set_xlabel('Minutes after first calibration pulse', fontsize=14)


    ax0.tick_params(length=6, width=1, which='minor')
    ax1.tick_params(length=6, width=1, which='minor')
    ax2.tick_params(length=6, width=1, which='minor')
    ax3.tick_params(length=6, width=1, which='minor')
    ax4.tick_params(length=6, width=1, which='minor')
    ax5.tick_params(length=6, width=1, which='minor')

    ax1.yaxis.set_label_coords(-0.09, 0.5)
    ax4.yaxis.set_label_coords(-0.09, 0.5)

    ax0.annotate(xy=(1.01,0.5), s='MHZ', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax1.annotate(xy=(1.01,0.5), s='MH1', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax2.annotate(xy=(1.01,0.5), s='MH2', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    ax3.annotate(xy=(1.01,0.5), s='MHZ', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax4.annotate(xy=(1.01,0.5), s='MH1', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')
    ax5.annotate(xy=(1.01,0.5), s='MH2', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    xticklabels = (ax0.get_xticklabels() + ax1.get_xticklabels()
      + ax3.get_xticklabels() + ax4.get_xticklabels())
    plt.setp(xticklabels, visible=False)



    plt.subplots_adjust(left=0.1, right=0.92, top=0.9, bottom=0.12)
    plt.savefig('../extra_plots_output/calibration_examples.png')
    plt.show()

def times_to_minutes(times_in_seconds):
    return ((times_in_seconds - 60) / 60)

if __name__ == "__main__":
    calibration_example()
