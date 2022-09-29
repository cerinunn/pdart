#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import stream_from_directory_new
from obspy import read_inventory
import os
from obspy.core import read
from pdart.util import linear_interpolation, remove_negative_ones_trace
    # start_time = UTCDateTime('1971-02-07T00:45:00')

# from pdart.diffusion.view_single_seismogram import remove_response

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

    stream = stream_from_directory_new(
      top_level_dir=top_level_dir,
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
    ax0.annotate(onset.strftime("%Y-%m-%d %H:%M:%S"), xy=(0.01,0.9),
      fontsize=13, horizontalalignment="left", verticalalignment="top",
      xycoords="axes fraction")

    xticklabels = (ax0.get_xticklabels())
    plt.setp(xticklabels, visible=False)

    ax0.tick_params(length=6, width=1, which='minor')

    ax0.yaxis.set_label_coords(-0.04, 0.5)

    ax0.annotate('MHZ', xy=(1.01,0.5), fontsize=16,
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



def single_seismogram_remove_response(top_level_dir,title,onset,pick=None,plot_response=False):

    seconds_before_arrival_time = 2*60
    # seconds_before_arrival_time = 0
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0.xml"
    onset = UTCDateTime('1969-11-20TT22:17:17.700000Z')
    start_time = onset - timedelta(seconds=seconds_before_arrival_time)
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

    stream =  remove_response_from_seismogram(
          top_level_dir=top_level_dir,
          inv_name=inv_name,
          start_time=start_time,
          station=station,
          channel=channel,
          pre_filt=pre_filt,
          water_level=None,
          end_time=end_time,
          plot_response=plot_response,
          plot=False)

    # stream.plot()
    # exit()

    ########

    fig = plt.figure(figsize=(6.5, 2.5))
    gs = gridspec.GridSpec(1, 1, hspace=0.001)

    ax0 = plt.subplot(gs[0])

    # seconds_before_arrival_time = 0
    trace_MHZ = stream.select(channel='MHZ')[0]
    ax0.plot((trace_MHZ.times() - seconds_before_arrival_time), trace_MHZ.data, color='k')

    
    ax0.set_xticks(np.arange(0,61*60,300),minor=False)
    ax0.set_xticks(np.arange(-180,61*60,60), minor=True)

    ax0.set_xlim(-seconds_before_arrival_time, 60*60)

    # pick_markP = pick - onset
    # plt.gca().axvline(x=pick_markP, 
    #                       color='r', linewidth=2)

    ax0.set_title(title, fontsize=12)

    # ax0.set_yticks(np.arange(480, 560, 20))
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    ax0.set_ylim(-0.51e-8, 0.51e-8)
    ax0.set_ylabel('Displacement [m]', fontsize=12, labelpad=0)
    ax0.annotate(onset.strftime("%Y-%m-%d %H:%M:%S"), xy=(0.01,0.97),
      fontsize=12, horizontalalignment="left", verticalalignment="top",
      xycoords="axes fraction")

    # xticklabels = (ax0.get_xticklabels())
    # plt.setp(xticklabels, visible=False)

    ax0.tick_params(length=3, width=1, which='minor')
    ax0.tick_params(length=6, width=1, which='major')

    # ax0.yaxis.set_label_coords(10, 0.5)

    ax0.annotate('MHZ', xy=(1.01,0.5), fontsize=12,
      xycoords="axes fraction", horizontalalignment='left',
      verticalalignment='center')

    ax0.set_xlabel('Time after impact [s]', fontsize=12, labelpad=0)
    # ax0.yaxis.set_label_coords(-0.04, 0.5)

    plt.subplots_adjust(left=0.08, right=0.93, top=0.9, bottom=0.17)
    plt.savefig('../extra_plots_output/Apollo12_LM_impact_XXXX.png')
    plt.show()

# def times_to_minutes(times_in_seconds):
#     return ((times_in_seconds / 60) - 2)



# copied from /Users/cnunn/python_packages/pdart/extra_plots/view_response.py
def remove_response_from_seismogram(
    top_level_dir,
    inv_name,
    start_time,
    station,
    channel,
    pre_filt,
    end_time=None,
    outfile=None,
    output='DISP',
    water_level=None,
    plot_response=False,
    plot=True):

    # read the response file
    inv = read_inventory(inv_name)


    if end_time is None:
        time_interval = timedelta(hours=3)
        end_time = start_time + time_interval
# xa.s12..att.1969.324.0.mseed

    # filename = '%s.%s.*.%s.%s.%03d.0.MSEED' % ('XA',station.upper(), channel.upper(),
    #    str(start_time.year), start_time.julday)
    # filename = os.path.join(top_level_dir,station.upper(),str(start_time.year),str(start_time.julday),filename)



    print(start_time, end_time)


    # include in a two minute buffer 
    stream = stream_from_directory_new(
      top_level_dir=top_level_dir,
      start_time=start_time-2*60,
      stations=[station],
      channels=[channel],
      end_time=end_time+ 2*60)

    for tr in stream:
        remove_negative_ones_trace(tr)
        # interpolate across the gaps of one sample 
        linear_interpolation(tr,interpolation_limit=1)
    stream.merge()


    for tr in stream:
        # optionally interpolate across any gap 
        # for removing the instrument response from a seimogram, 
        # it is useful to get a mask, then interpolate across the gaps, 
        # then mask the trace again. 
        if tr.stats.channel in ['MH1', 'MH2', 'MHZ']:

            # add linear interpolation but keep the original mask
            original_mask = linear_interpolation(tr,interpolation_limit=None)
            # remove the instrument response
            pre_filt = [0.1,0.3,0.9,1.1]
            tr.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP",
                       water_level=None, plot=plot_response)
            if plot_response:
                plt.show()
            # apply the mask back to the trace 
            tr.data = np.ma.masked_array(tr, mask=original_mask)

        elif tr.stats.channel in ['SHZ']:

            # add linear interpolation but keep the original mask
            original_mask = linear_interpolation(tr,interpolation_limit=None)
            # remove the instrument response
            pre_filt = [1,2,11,13] 
            tr.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP",
                       water_level=None, plot=plot_response)
            if plot_response:
                plt.show()
            
            # apply the mask back to the trace 
            tr.data = np.ma.masked_array(tr, mask=original_mask)

            # if mask is not None:
            #     stream[0].data = np.ma.array(stream[0].data, mask = mask)

    # remove the two minute buffer 
    stream.trim(start_time,end_time)

    # stream = read(filename)
    # stream = stream.select(channel=channel)
    # stream.trim(starttime=start_time, endtime=end_time)

    # remove location (ground station)
    # for tr in stream:
    #     tr.stats.location = ''

    # detrend
    # stream.detrend('linear')

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

    # # find the gaps in the trace
    # if isinstance(stream[0].data,np.ma.MaskedArray):
    #     mask = np.ma.getmask(stream[0].data)
    # else:
    #     mask = None
    # 
    # # split the stream, then refill it with zeros on the gaps
    # stream = stream.split()
    # stream = stream.merge(fill_value=0)

    # for i, n in enumerate(stream[0].times()):
    #     # print(n)
    #     stream[0].data[i]=np.sin(2*np.pi*(1/25)*n)

    # stream.attach_response(inv)
    # # print('here')
    # 
    # # zero_mean=False - because the trace can be asymmetric - remove the mean ourselves
    # # do not taper here - it doesn't work well with the masked arrays - often required
    # # when there are gaps - if necessary taper first
    # # water level - this probably doesn't have much impact - because we are pre filtering
    # # stream.remove_response(pre_filt=pre_filt,output="DISP",water_level=30,zero_mean=False,taper=False,plot=True,fig=outfile)
    # for tr in stream:
    #     tr.remove_response(pre_filt=pre_filt,output=output,water_level=water_level,zero_mean=False,taper=False,plot=plot,fig=outfile)

    # for tr in stream:
    #     tr.stats.location = 'changed'



    # print(stream)

    return stream

# def times_to_seconds(times_in_seconds, seconds_before_arrival_time):
#     return (times_in_seconds - time_delta(seconds=seconds_before_arrival_time))

if __name__ == "__main__":
    top_level_dir='/Users/cnunn/lunar_data/PDART_V2'
    # single_seismogram(title='Impact of Apollo 12 Lunar Ascent Module')
    onset = UTCDateTime('1969-11-20TT22:17:17.700000Z')
      # <pick publicID="smi:nunn19/pick/00001/lognonne03/S12/P">
    pick = UTCDateTime('1969-11-20T22:17:42.400000Z')
    arrival_time = pick - onset
    print(arrival_time)

    single_seismogram_remove_response(top_level_dir=top_level_dir,title='Impact of Apollo 12 Lunar Ascent Module',onset=onset,pick=pick)

    # single_seismogram_remove_response_short(title='Impact of Apollo 12 Lunar Ascent Module')


