#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import stream_from_directory
from obspy import read_inventory
import os
from obspy.core import read
    # start_time = UTCDateTime('1971-02-07T00:45:00')

from pdart.diffusion.view_single_seismogram import remove_response

# update 25-08-21
def single_seismogram(title):

# 1969-11-20T22:17:17.7
    # onset is 42.4 s Lognonne 2003
    onset = UTCDateTime('1969-11-20TT22:17:17.700000Z') 
# +42.4
    # onset = UTCDateTime('1969-11-20T22:17:700000Z')
    start_time = onset - timedelta(minutes=2)

    stations = ['S12']
    channels = ['MHZ']
    end_time = start_time + timedelta(minutes=15)

    stream = stream_from_directory(
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
    ax0.set_title(title, fontsize=20)

    ax0.set_yticks(np.arange(480, 560, 20))
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    ax0.set_ylim(510-50, 510+50)
    ax0.set_ylabel('DU', fontsize=14)
    ax0.annotate(xy=(0.01,0.9), text=onset.strftime("%Y-%m-%d %H:%M:%S"),
      fontsize=13, horizontalalignment="left", verticalalignment="top",
      xycoords="axes fraction")

    xticklabels = (ax0.get_xticklabels())
    plt.setp(xticklabels, visible=False)

    ax0.tick_params(length=6, width=1, which='minor')

    ax0.yaxis.set_label_coords(-0.04, 0.5)

    ax0.annotate(xy=(1.01,0.5), text='MHZ', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    plt.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.12)
    plt.savefig('Apollo12_LM_impact_XXXX.png')
    plt.show()

# def single_seismogram_remove_response_short(title):
# 
#     # peaked mode
#     inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
#     # onset is 42.4 s Lognonne 2003
#     onset = UTCDateTime('1969-11-20TT22:17:17.700000Z')
#     start_time = onset - timedelta(minutes=2)
#     # XXXX
#     station = 'S12'
#     channel = 'MHZ'
#     pre_filt = [0.1, 0.3,0.7,1]
#     # end_time = UTCDateTime('1971:02:07T02:35.25')
# 
#     end_time = onset + timedelta(minutes=60)
# 
# # 1969-11-20T22:17:17.7
# 
#     stream =  remove_response_from_seismogram(inv_name=inv_name,
#           start_time=start_time,
#           station=station,
#           channel=channel,
#           pre_filt=pre_filt,
#           water_level=None,
#           end_time=end_time,
#           plot=False)
# 
#     ########
# 
#     fig = plt.figure(figsize=(15, 4))
#     gs = gridspec.GridSpec(1, 1, hspace=0.001)
# 
#     ax0 = plt.subplot(gs[0])
# 
#     trace_MHZ = stream.select(channel='MHZ')[0]
#     ax0.plot(times_to_seconds(trace_MHZ.times()), trace_MHZ.data, color='k')
#     ax0.set_xlim(-2*60, 5*60)
#     # print('short')
#     ax0.set_xticks(np.arange(0,6*60,6*50),minor=False)
#     ax0.set_xticks(np.arange(-180,6*60,60), minor=True)
#     ax0.set_title(title, fontsize=20)
# 
#     # ax0.set_yticks(np.arange(480, 560, 20))
#     # ax0.set_yticks(np.arange(460,580,20), minor=True)
#     ax0.set_ylim(-1.1e-8, 1.01e-8)
#     ax0.set_ylabel('Displacement [m]', fontsize=14)
#     ax0.annotate(xy=(0.01,0.9), text=onset.strftime("%Y-%m-%d %H:%M:%S"),
#       fontsize=13, horizontalalignment="left", verticalalignment="top",
#       xycoords="axes fraction")
# 
#     # xticklabels = (ax0.get_xticklabels())
#     # plt.setp(xticklabels, visible=False)
# 
#     ax0.tick_params(length=3, width=1, which='minor')
#     ax0.tick_params(length=6, width=1, which='major')
# 
#     ax0.yaxis.set_label_coords(-0.04, 0.5)
# 
#     ax0.annotate(xy=(1.01,0.5), text='MHZ', fontsize=16,
#       xycoords="axes fraction", horizontalalignment='left',
#       verticalalignment='center')
# 
#     ax0.set_xlabel('Time after impact [s]', fontsize=14)
#     ax0.yaxis.set_label_coords(-0.04, 0.5)
# 
#     plt.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.12)
#     plt.savefig('Apollo12_LM_impact_XXXX.png')
#     plt.show()

def single_seismogram_remove_response(title,onset,pick=None):

    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    onset = UTCDateTime('1969-11-20TT22:17:17.700000Z')
    start_time = onset - timedelta(minutes=2)
    # XXXX
    station = 'S12'
    channel = 'MHZ'
    pre_filt = [0.1, 0.3,0.7,1]
    # end_time = UTCDateTime('1971:02:07T02:35.25')

    end_time = onset + timedelta(minutes=60)

# 1969-11-20T22:17:17.7

    # reset the timing 
    

    # make a correction
    # find actual time of onset 
    # print(onset.time)

    stream =  remove_response_from_seismogram(inv_name=inv_name,
          start_time=start_time,
          station=station,
          channel=channel,
          pre_filt=pre_filt,
          water_level=None,
          end_time=end_time,
          plot=False)

    ########

    fig = plt.figure(figsize=(15, 4))
    gs = gridspec.GridSpec(1, 1, hspace=0.001)

    ax0 = plt.subplot(gs[0])

    trace_MHZ = stream.select(channel='MHZ')[0]
    ax0.plot(times_to_seconds(trace_MHZ.times()), trace_MHZ.data, color='k')
    ax0.set_xlim(-2*60, 60*60)
    ax0.set_xticks(np.arange(0,61*60,6*50),minor=False)
    ax0.set_xticks(np.arange(-180,61*60,60), minor=True)

    # pick_markP = pick - onset
    # plt.gca().axvline(x=pick_markP, 
    #                       color='r', linewidth=2)

    ax0.set_title(title, fontsize=20)

    # ax0.set_yticks(np.arange(480, 560, 20))
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    ax0.set_ylim(-1.1e-8, 1.01e-8)
    ax0.set_ylabel('Displacement [m]', fontsize=14)
    ax0.annotate(xy=(0.01,0.9), text=onset.strftime("%Y-%m-%d %H:%M:%S"),
      fontsize=13, horizontalalignment="left", verticalalignment="top",
      xycoords="axes fraction")

    # xticklabels = (ax0.get_xticklabels())
    # plt.setp(xticklabels, visible=False)

    ax0.tick_params(length=3, width=1, which='minor')
    ax0.tick_params(length=6, width=1, which='major')

    ax0.yaxis.set_label_coords(-0.04, 0.5)

    ax0.annotate(xy=(1.01,0.5), text='MHZ', fontsize=16,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    ax0.set_xlabel('Time after impact [s]', fontsize=14)
    ax0.yaxis.set_label_coords(-0.04, 0.5)

    ax0.plot(times_to_seconds(trace_MHZ.times()),
      times_to_seconds(trace_MHZ.data-trace_MHZ.stats.starttime.timestamp), color='k')

    plt.subplots_adjust(left=0.06, right=0.95, top=0.9, bottom=0.13)
    plt.savefig('../extra_plots_output/Apollo12_LM_impact_XXXX.png')
    plt.show()

# def times_to_minutes(times_in_seconds):
#     return ((times_in_seconds / 60) - 2)



# copied from /Users/cnunn/python_packages/pdart/extra_plots/view_response.py
def remove_response_from_seismogram(
    inv_name,
    start_time,
    station,
    channel,
    pre_filt,
    end_time=None,
    outfile=None,
    output='DISP',
    water_level=None,
    plot=True):

    # read the response file
    inv = read_inventory(inv_name)

    if end_time is None:
        time_interval = timedelta(hours=3)
        end_time = start_time + time_interval
# xa.s12..att.1969.324.0.mseed

    filename = '%s.%s.*.%s.%s.%03d.0.mseed' % ('xa',station.lower(), channel.lower(),
       str(start_time.year), start_time.julday)
    filename = os.path.join('/Users/cnunn/lunar_data/PDART_CONTINUOUS_MAIN_TAPES',station.lower(),str(start_time.year),str(start_time.julday),filename)

    stream = read(filename)
    stream = stream.select(channel=channel)
    stream.trim(starttime=start_time, endtime=end_time)

    # remove location (ground station)
    for tr in stream:
        tr.stats.location = ''

    # detrend
    stream.detrend('linear')

    # taper the edges
    # if there are gaps in the seismogram - EVERY short trace will be tapered
    # this is required to remove the response later
    # stream.taper(max_percentage=0.05, type='cosine')

    # experiment with tapering? not tapering preserves the overall shape better
    # but it may required

    # merge the streams
    stream.merge()
    if stream.count() > 1:
        print('Too many streams - exiting')

    # find the gaps in the trace
    if isinstance(stream[0].data,np.ma.MaskedArray):
        mask = np.ma.getmask(stream[0].data)
    else:
        mask = None

    # split the stream, then refill it with zeros on the gaps
    stream = stream.split()
    stream = stream.merge(fill_value=0)

    # for i, n in enumerate(stream[0].times()):
    #     # print(n)
    #     stream[0].data[i]=np.sin(2*np.pi*(1/25)*n)

    stream.attach_response(inv)
    # print('here')

    # zero_mean=False - because the trace can be asymmetric - remove the mean ourselves
    # do not taper here - it doesn't work well with the masked arrays - often required
    # when there are gaps - if necessary taper first
    # water level - this probably doesn't have much impact - because we are pre filtering
    # stream.remove_response(pre_filt=pre_filt,output="DISP",water_level=30,zero_mean=False,taper=False,plot=True,fig=outfile)
    for tr in stream:
        remove_response(tr, pre_filt=pre_filt,output=output,water_level=water_level,zero_mean=False,taper=False,plot=plot,fig=outfile)

    for tr in stream:
        tr.stats.location = 'changed'

    if mask is not None:
        stream[0].data = np.ma.array(stream[0].data, mask = mask)

    print(stream)

    return stream

def times_to_seconds(times_in_seconds):
    return (times_in_seconds - 120)

if __name__ == "__main__":
    # single_seismogram(title='Impact of Apollo 12 Lunar Ascent Module')
    onset = UTCDateTime('1969-11-20TT22:17:17.700000Z')
      # <pick publicID="smi:nunn19/pick/00001/lognonne03/S12/P">
    pick = UTCDateTime('1969-11-20T22:17:42.400000Z')
    arrival_time = pick - onset
    print(arrival_time)

    single_seismogram_remove_response(title='Impact of Apollo 12 Lunar Ascent Module',onset=onset,pick=pick)

    # single_seismogram_remove_response_short(title='Impact of Apollo 12 Lunar Ascent Module')


