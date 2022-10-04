#!/usr/bin/env python

'''
Plot example moonquakes with instrument removed.


'''

from __future__ import print_function
import numpy as np
import os

import matplotlib  
from matplotlib import pyplot as plt

from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.util import linear_interpolation, timing_correction
from pdart.view import stream_from_directory_new
from obspy.core.inventory import read_inventory
from obspy.core import read

SECONDS_PER_DAY = 3600.0 * 24.0

# From the Tol color palette
# https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
tol_dark_blue ='#341B88' 
midnight_blue = '#191970'

def view_Apollo(starttime= UTCDateTime('1973-03-13T07:30:00.0'),endtime = UTCDateTime('1973-03-13T09:30:00.0'),
  network='XA',station='S14',channel='MH1',location='*',plot_seismogram=True,plot_response=False):
    """Get a stream, remove the instrument response, and return the stream. 
    
    """

    inv = read_inventory('/Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0.xml')

    stream = stream_from_directory_new(
        top_level_dir = '/Users/cnunn/lunar_data/PDART_V2',
        start_time=starttime,
        stations=[station],
        channels=[channel],
        merge_locations=False,
        end_time=endtime)

    # # get the response file (wildcards allowed)
    # inv = client.get_stations(starttime=starttime, endtime=endtime,
    #     network=network, sta=station, loc=location, channel=channel,
    #     level="response")
    
    for tr in stream:
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

    if plot_seismogram:
        stream.plot(equal_scale=False,size=(1000,600),method='full')

    return stream 



def plot_examples():


#  <text>deep moonquake</text>
# 1973-06-05T11
# 
# <type>meteorite</type>
# 1972-05-11T13
# 
#  <text>shallow moonquake</text>
# 1973-03-13T08
# 
# <text>S-IVB impact</text>
# 1971-02-04T07


    onset_deep=UTCDateTime("1973-06-05T11:12:18")
    onset_meteoroid=UTCDateTime("1972-05-11T13:34:37")
    onset_shallow=UTCDateTime("1973-03-13T08:01:33")
    onset_artifical=UTCDateTime("1971-02-04T07:41:32")
    
    
    deep = view_Apollo(starttime=onset_deep-700,endtime=onset_deep+7300,channel='MHZ',station='S12',plot_seismogram=False)
    deep = deep.trim(starttime=onset_deep-600,endtime=onset_deep+7200)
    shallow = view_Apollo(starttime=onset_shallow-700,endtime=onset_shallow+7300,channel='MHZ',station='S12',plot_seismogram=False)
    shallow = shallow.trim(starttime=onset_shallow-600,endtime=onset_shallow+7200)
    meteoroid = view_Apollo(starttime=onset_meteoroid-700,endtime=onset_meteoroid+7300,channel='MHZ',station='S12',plot_seismogram=False)
    meteoroid = meteoroid.trim(starttime=onset_meteoroid-600,endtime=onset_meteoroid+7200)
    artificial = view_Apollo(starttime=onset_artifical-700,endtime=onset_artifical+7300,channel='MHZ',station='S12',plot_seismogram=False)
    artificial = artificial.trim(starttime=onset_artifical-600,endtime=onset_artifical+7200)

    fig = plt.figure(figsize=(6.2, 4))
    gs = gridspec.GridSpec(2, 2, hspace=0.2,wspace=0.2)

    # top left - Deep Moonquake 
    ax0 = plt.subplot(gs[0])

    ax0.plot(times_to_minutes(deep[0].times()), deep[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-10, 120)
    ax0.set_xticks([0,60,120])
    ax0.set_xticks(np.arange(-10,120,10), minor=True)
    ax0.axes.xaxis.set_ticklabels([])
    ax0.set_title('Deep Moonquake', fontsize=13, pad=-1)

    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)
    ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    ax0.annotate(onset_deep.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    # top right - Meteoroid  
    ax0 = plt.subplot(gs[1])

    ax0.plot(times_to_minutes(meteoroid[0].times()), meteoroid[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-10, 120)
    ax0.set_xticks([0,60,120])
    ax0.set_xticks(np.arange(-10,120,10), minor=True)
    ax0.axes.xaxis.set_ticklabels([])
    ax0.set_title('Meteoroid Impact', fontsize=13, pad=-1)

    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)

    ax0.annotate(onset_meteoroid.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    # bottom left - Shallow 
    ax0 = plt.subplot(gs[2])

    ax0.plot(times_to_minutes(shallow[0].times()), shallow[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-10, 120)
    ax0.set_xticks([0,60,120])
    ax0.set_xticks(np.arange(-10,120,10), minor=True)
    ax0.set_title('Shallow Moonquake', fontsize=13, pad=-1)

    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)
    ax0.set_xlabel('Minutes afer arrival time', fontsize=10)
    ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    ax0.annotate(onset_shallow.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    # bottom right - Artificial  
    ax0 = plt.subplot(gs[3])

    ax0.plot(times_to_minutes(artificial[0].times()), artificial[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-10, 120)
    ax0.set_xticks([0,60,120])
    ax0.set_xticks(np.arange(-10,120,10), minor=True)
    ax0.set_title('Artifical Impact', fontsize=13, pad=-1)

    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)
    ax0.set_xlabel('Minutes after arrival time', fontsize=10)
    ax0.annotate(onset_artifical.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    plt.subplots_adjust(left=0.1, right=0.975, top=0.93, bottom=0.11)
    plt.savefig('../extra_plots_output/moonquake_examples_MHZ.png', dpi=300)
    plt.show()

def plot_deep():


    onset_deep=UTCDateTime("1973-06-05T11:12:18")
    onset_meteoroid=UTCDateTime("1972-05-11T13:34:37")
    onset_shallow=UTCDateTime("1973-03-13T08:01:33")
    onset_artifical=UTCDateTime("1971-02-04T07:41:32")
    
    
    deep = view_Apollo(starttime=onset_deep-700,endtime=onset_deep+7300,channel='MHZ',station='S12',plot_seismogram=False)
    deep = deep.trim(starttime=onset_deep-600,endtime=onset_deep+7200)
    # shallow = view_Apollo(starttime=onset_shallow-700,endtime=onset_shallow+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # shallow = shallow.trim(starttime=onset_shallow-600,endtime=onset_shallow+7200)
    # meteoroid = view_Apollo(starttime=onset_meteoroid-700,endtime=onset_meteoroid+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # meteoroid = meteoroid.trim(starttime=onset_meteoroid-600,endtime=onset_meteoroid+7200)
    # artificial = view_Apollo(starttime=onset_artifical-700,endtime=onset_artifical+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # artificial = artificial.trim(starttime=onset_artifical-600,endtime=onset_artifical+7200)

    fig = plt.figure(figsize=(6.2, 4))
    gs = gridspec.GridSpec(1, 1, hspace=0.2,wspace=0.2)

    # top left - Deep Moonquake 
    ax0 = plt.subplot(gs[0])

    ax0.plot(times_to_minutes(deep[0].times()), deep[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-10, 120)
    ax0.set_xticks([0,60,120])
    ax0.set_xticks(np.arange(-10,120,10), minor=True)
    ax0.set_title('Deep Moonquake', fontsize=13, pad=-1)

    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)
    ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    ax0.annotate(onset_deep.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    # # top right - Shallow Moonquake 
    # ax0 = plt.subplot(gs[1])
    # 
    # ax0.plot(times_to_minutes(shallow[0].times()), shallow[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.axes.xaxis.set_ticklabels([])
    # ax0.set_title('Shallow Moonquake', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # 
    # ax0.annotate(onset_shallow.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # bottom left - Meteoroid 
    # ax0 = plt.subplot(gs[2])
    # 
    # ax0.plot(times_to_minutes(meteoroid[0].times()), meteoroid[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Meteoroid Impact', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_xlabel('Minutes afer arrival time', fontsize=10)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # ax0.annotate(onset_meteoroid.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # # bottom right - Artificial  
    # ax0 = plt.subplot(gs[3])
    # 
    # ax0.plot(times_to_minutes(artificial[0].times()), artificial[0].data*10**9, color=tol_dark_blue,linewidth=1)

    # ax0.set_title('Artifical Impact', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)

    ax0.set_xlabel('Minutes after arrival time', fontsize=10)
    # ax0.annotate(onset_artifical.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    plt.subplots_adjust(left=0.12, right=0.975, top=0.93, bottom=0.11)
    plt.savefig('../extra_plots_output/moonquake_deep_MHZ_XXXX.png')
    # plt.show()

def plot_shallow():


    onset_deep=UTCDateTime("1973-06-05T11:12:18")
    onset_meteoroid=UTCDateTime("1972-05-11T13:34:37")
    onset_shallow=UTCDateTime("1973-03-13T08:01:33")
    onset_artifical=UTCDateTime("1971-02-04T07:41:32")
    
    
    # deep = view_Apollo(starttime=onset_deep-700,endtime=onset_deep+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # deep = deep.trim(starttime=onset_deep-600,endtime=onset_deep+7200)
    shallow = view_Apollo(starttime=onset_shallow-700,endtime=onset_shallow+7300,channel='MHZ',station='S12',plot_seismogram=False)
    shallow = shallow.trim(starttime=onset_shallow-600,endtime=onset_shallow+7200)
    # meteoroid = view_Apollo(starttime=onset_meteoroid-700,endtime=onset_meteoroid+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # meteoroid = meteoroid.trim(starttime=onset_meteoroid-600,endtime=onset_meteoroid+7200)
    # artificial = view_Apollo(starttime=onset_artifical-700,endtime=onset_artifical+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # artificial = artificial.trim(starttime=onset_artifical-600,endtime=onset_artifical+7200)

    fig = plt.figure(figsize=(6.2, 4))
    gs = gridspec.GridSpec(1, 1, hspace=0.2,wspace=0.2)

    # # top left - Deep Moonquake 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(deep[0].times()), deep[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Deep Moonquake', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # ax0.annotate(onset_deep.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # top right - Shallow Moonquake 
    ax0 = plt.subplot(gs[0])
    
    ax0.plot(times_to_minutes(shallow[0].times()), shallow[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-10, 120)
    ax0.set_xticks([0,60,120])
    ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.axes.xaxis.set_ticklabels([])
    ax0.set_title('Shallow Moonquake', fontsize=13, pad=-1)
    # 
    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)
    ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    
    ax0.annotate(onset_shallow.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    # bottom left - Meteoroid 
    # ax0 = plt.subplot(gs[2])
    # 
    # ax0.plot(times_to_minutes(meteoroid[0].times()), meteoroid[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Meteoroid Impact', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_xlabel('Minutes afer arrival time', fontsize=10)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # ax0.annotate(onset_meteoroid.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # # bottom right - Artificial  
    # ax0 = plt.subplot(gs[3])
    # 
    # ax0.plot(times_to_minutes(artificial[0].times()), artificial[0].data*10**9, color=tol_dark_blue,linewidth=1)

    # ax0.set_title('Artifical Impact', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)

    ax0.set_xlabel('Minutes after arrival time', fontsize=10)
    # ax0.annotate(onset_artifical.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # plt.subplots_adjust(left=0.12, right=0.975, top=0.93, bottom=0.11)
    plt.savefig('../extra_plots_output/moonquake_shallow_MHZ_XXXX.png')
    # plt.show()


def plot_meteoroid():


    onset_deep=UTCDateTime("1973-06-05T11:12:18")
    onset_meteoroid=UTCDateTime("1972-05-11T13:34:37")
    onset_shallow=UTCDateTime("1973-03-13T08:01:33")
    onset_artifical=UTCDateTime("1971-02-04T07:41:32")
    
    
    # deep = view_Apollo(starttime=onset_deep-700,endtime=onset_deep+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # deep = deep.trim(starttime=onset_deep-600,endtime=onset_deep+7200)
    # shallow = view_Apollo(starttime=onset_shallow-700,endtime=onset_shallow+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # shallow = shallow.trim(starttime=onset_shallow-600,endtime=onset_shallow+7200)
    meteoroid = view_Apollo(starttime=onset_meteoroid-700,endtime=onset_meteoroid+7300,channel='MHZ',station='S12',plot_seismogram=False)
    meteoroid = meteoroid.trim(starttime=onset_meteoroid-600,endtime=onset_meteoroid+7200)
    # artificial = view_Apollo(starttime=onset_artifical-700,endtime=onset_artifical+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # artificial = artificial.trim(starttime=onset_artifical-600,endtime=onset_artifical+7200)

    fig = plt.figure(figsize=(6.2, 4))
    gs = gridspec.GridSpec(1, 1, hspace=0.2,wspace=0.2)

    # # top left - Deep Moonquake 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(deep[0].times()), deep[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Deep Moonquake', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # ax0.annotate(onset_deep.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # # top right - Shallow Moonquake 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(shallow[0].times()), shallow[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # # ax0.axes.xaxis.set_ticklabels([])
    # ax0.set_title('Shallow Moonquake', fontsize=13, pad=-1)
    # # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # 
    # ax0.annotate(onset_shallow.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # bottom left - Meteoroid 
    ax0 = plt.subplot(gs[0])
    
    ax0.plot(times_to_minutes(meteoroid[0].times()), meteoroid[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-10, 120)
    ax0.set_xticks([0,60,120])
    ax0.set_xticks(np.arange(-10,120,10), minor=True)
    ax0.set_title('Meteoroid Impact', fontsize=13, pad=-1)
    
    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)
    ax0.set_xlabel('Minutes afer arrival time', fontsize=10)
    ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    ax0.annotate(onset_meteoroid.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    # # bottom right - Artificial  
    # ax0 = plt.subplot(gs[3])
    # 
    # ax0.plot(times_to_minutes(artificial[0].times()), artificial[0].data*10**9, color=tol_dark_blue,linewidth=1)

    # ax0.set_title('Artifical Impact', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)

    ax0.set_xlabel('Minutes after arrival time', fontsize=10)
    # ax0.annotate(onset_artifical.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # plt.subplots_adjust(left=0.12, right=0.975, top=0.93, bottom=0.11)
    plt.savefig('../extra_plots_output/moonquake_meteoroid_MHZ_XXXX.png')
    # plt.show()


def plot_artificial():


    onset_deep=UTCDateTime("1973-06-05T11:12:18")
    onset_meteoroid=UTCDateTime("1972-05-11T13:34:37")
    onset_shallow=UTCDateTime("1973-03-13T08:01:33")
    onset_artifical=UTCDateTime("1971-02-04T07:41:32")
    
    
    # deep = view_Apollo(starttime=onset_deep-700,endtime=onset_deep+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # deep = deep.trim(starttime=onset_deep-600,endtime=onset_deep+7200)
    # shallow = view_Apollo(starttime=onset_shallow-700,endtime=onset_shallow+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # shallow = shallow.trim(starttime=onset_shallow-600,endtime=onset_shallow+7200)
    # meteoroid = view_Apollo(starttime=onset_meteoroid-700,endtime=onset_meteoroid+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # meteoroid = meteoroid.trim(starttime=onset_meteoroid-600,endtime=onset_meteoroid+7200)
    artificial = view_Apollo(starttime=onset_artifical-700,endtime=onset_artifical+7300,channel='MHZ',station='S12',plot_seismogram=False)
    artificial = artificial.trim(starttime=onset_artifical-600,endtime=onset_artifical+7200)

    fig = plt.figure(figsize=(6.2, 4))
    gs = gridspec.GridSpec(1, 1, hspace=0.2,wspace=0.2)

    # # top left - Deep Moonquake 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(deep[0].times()), deep[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Deep Moonquake', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # ax0.annotate(onset_deep.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # # top right - Shallow Moonquake 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(shallow[0].times()), shallow[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # # ax0.axes.xaxis.set_ticklabels([])
    # ax0.set_title('Shallow Moonquake', fontsize=13, pad=-1)
    # # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # 
    # ax0.annotate(onset_shallow.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # # bottom left - Meteoroid 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(meteoroid[0].times()), meteoroid[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Meteoroid Impact', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_xlabel('Minutes afer arrival time', fontsize=10)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # ax0.annotate(onset_meteoroid.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # bottom right - Artificial  
    ax0 = plt.subplot(gs[0])
    
    ax0.plot(times_to_minutes(artificial[0].times()), artificial[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-10, 120)
    ax0.set_xticks([0,60,120])
    ax0.set_xticks(np.arange(-10,120,10), minor=True)

    ax0.set_title('Artifical Impact', fontsize=13, pad=-1)
    
    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)

    ax0.set_xlabel('Minutes after arrival time', fontsize=10)
    ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    ax0.annotate(onset_artifical.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    # plt.subplots_adjust(left=0.12, right=0.975, top=0.93, bottom=0.11)
    plt.savefig('../extra_plots_output/moonquake_artificial_MHZ_XXXX.png')
    # plt.show()


def plot_artificial_ten_minutes():


    onset_deep=UTCDateTime("1973-06-05T11:12:18")
    onset_meteoroid=UTCDateTime("1972-05-11T13:34:37")
    onset_shallow=UTCDateTime("1973-03-13T08:01:33")
    onset_artifical=UTCDateTime("1971-02-04T07:41:32")
    
    
    # deep = view_Apollo(starttime=onset_deep-700,endtime=onset_deep+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # deep = deep.trim(starttime=onset_deep-600,endtime=onset_deep+7200)
    # shallow = view_Apollo(starttime=onset_shallow-700,endtime=onset_shallow+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # shallow = shallow.trim(starttime=onset_shallow-600,endtime=onset_shallow+7200)
    # meteoroid = view_Apollo(starttime=onset_meteoroid-700,endtime=onset_meteoroid+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # meteoroid = meteoroid.trim(starttime=onset_meteoroid-600,endtime=onset_meteoroid+7200)
    artificial = view_Apollo(starttime=onset_artifical-700,endtime=onset_artifical+7300,channel='MHZ',station='S12',plot_seismogram=False)
    artificial = artificial.trim(starttime=onset_artifical-600,endtime=onset_artifical+7200)

    fig = plt.figure(figsize=(6.2, 4))
    gs = gridspec.GridSpec(1, 1, hspace=0.2,wspace=0.2)

    # # top left - Deep Moonquake 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(deep[0].times()), deep[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Deep Moonquake', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # ax0.annotate(onset_deep.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # # top right - Shallow Moonquake 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(shallow[0].times()), shallow[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # # ax0.axes.xaxis.set_ticklabels([])
    # ax0.set_title('Shallow Moonquake', fontsize=13, pad=-1)
    # # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # 
    # ax0.annotate(onset_shallow.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # # bottom left - Meteoroid 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(meteoroid[0].times()), meteoroid[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Meteoroid Impact', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_xlabel('Minutes afer arrival time', fontsize=10)

    # ax0.annotate(onset_meteoroid.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # bottom right - Artificial  
    ax0 = plt.subplot(gs[0])
    
    ax0.plot(times_to_minutes(artificial[0].times()), artificial[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-1, 10)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)

    ax0.set_title('Artifical Impact', fontsize=13, pad=-1)
    
    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)

    ax0.set_xlabel('Minutes after arrival time', fontsize=10)
    ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    ax0.annotate(onset_artifical.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    # plt.subplots_adjust(left=0.12, right=0.975, top=0.93, bottom=0.11)
    plt.savefig('../extra_plots_output/moonquake_artificial_ten_minutes_MHZ_XXXX.png')
    # plt.show()


def plot_artificial_one_minute():


    onset_deep=UTCDateTime("1973-06-05T11:12:18")
    onset_meteoroid=UTCDateTime("1972-05-11T13:34:37")
    onset_shallow=UTCDateTime("1973-03-13T08:01:33")
    onset_artifical=UTCDateTime("1971-02-04T07:41:32")
    
    
    # deep = view_Apollo(starttime=onset_deep-700,endtime=onset_deep+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # deep = deep.trim(starttime=onset_deep-600,endtime=onset_deep+7200)
    # shallow = view_Apollo(starttime=onset_shallow-700,endtime=onset_shallow+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # shallow = shallow.trim(starttime=onset_shallow-600,endtime=onset_shallow+7200)
    # meteoroid = view_Apollo(starttime=onset_meteoroid-700,endtime=onset_meteoroid+7300,channel='MHZ',station='S12',plot_seismogram=False)
    # meteoroid = meteoroid.trim(starttime=onset_meteoroid-600,endtime=onset_meteoroid+7200)
    artificial = view_Apollo(starttime=onset_artifical-700,endtime=onset_artifical+7300,channel='MHZ',station='S12',plot_seismogram=False)
    artificial = artificial.trim(starttime=onset_artifical-600,endtime=onset_artifical+60)

    fig = plt.figure(figsize=(6.2, 4))
    gs = gridspec.GridSpec(1, 1, hspace=0.2,wspace=0.2)

    # # top left - Deep Moonquake 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(deep[0].times()), deep[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Deep Moonquake', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)

    # ax0.annotate(onset_deep.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # # top right - Shallow Moonquake 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(shallow[0].times()), shallow[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # # ax0.axes.xaxis.set_ticklabels([])
    # ax0.set_title('Shallow Moonquake', fontsize=13, pad=-1)
    # # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    # 
    # ax0.annotate(onset_shallow.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # # bottom left - Meteoroid 
    # ax0 = plt.subplot(gs[0])
    # 
    # ax0.plot(times_to_minutes(meteoroid[0].times()), meteoroid[0].data*10**9, color=tol_dark_blue,linewidth=1)
    # ax0.set_xlim(-10, 120)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)
    # ax0.set_title('Meteoroid Impact', fontsize=13, pad=-1)
    # 
    # # ax0.set_yticks([470,510,550])
    # # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # # ax0.set_ylim(510-80, 510+80)
    # ax0.set_xlabel('Minutes afer arrival time', fontsize=10)

    # ax0.annotate(onset_meteoroid.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
    #   fontsize=10, horizontalalignment="right", verticalalignment="bottom",
    #   xycoords="axes fraction")

    # bottom right - Artificial  
    ax0 = plt.subplot(gs[0])
    
    ax0.plot(times_to_minutes(artificial[0].times()), artificial[0].data*10**9, color=tol_dark_blue,linewidth=1)
    ax0.set_xlim(-0.5, 1)
    # ax0.set_xticks([0,60,120])
    # ax0.set_xticks(np.arange(-10,120,10), minor=True)


    ax0.set_title('Artifical Impact', fontsize=13, pad=-1)
    
    # ax0.set_yticks([470,510,550])
    # ax0.set_yticks(np.arange(460,580,20), minor=True)
    # ax0.set_ylim(510-80, 510+80)

    ax0.set_xlabel('Minutes after arrival time', fontsize=10)
    ax0.set_ylabel(r'x10$^9$ m', fontsize=11, labelpad=-2)
    ax0.annotate(onset_artifical.strftime("%Y-%m-%d %H:%M:%S") ,xy=(0.99,0.01),
      fontsize=10, horizontalalignment="right", verticalalignment="bottom",
      xycoords="axes fraction")

    # plt.subplots_adjust(left=0.12, right=0.975, top=0.93, bottom=0.11)
    plt.savefig('../extra_plots_output/moonquake_artificial_one_minute_MHZ_XXXX.png')
    # plt.show()


def times_to_minutes(times_in_seconds):
    return ((times_in_seconds / 60) - 10)

if __name__ == "__main__":

    # plot_examples()
    plot_deep()
    plot_shallow()
    plot_meteoroid()
    plot_artificial()
    plot_artificial_ten_minutes()
    plot_artificial_one_minute()