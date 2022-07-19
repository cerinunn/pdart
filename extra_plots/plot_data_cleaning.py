#!/usr/bin/env python
'''
Plot data cleaning and a digital spike. No longer used. 


'''


from __future__ import print_function
import numpy as np
import os

# Qt5Agg seems to work best on Mac - try 'TkAgg' if that works for you
# put this after the other imports, otherwise it can be overridden
import matplotlib  
# matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
# 
# from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
# from matplotlib import gridspec
# from pdart.view import stream_from_directory_new
# from obspy.imaging.util import (_set_xaxis_obspy_dates, _id_key, _timestring)
# from matplotlib.dates import date2num
# from pdart.util import relative_timing_trace, relative_timing_stream, remove_negative_ones
# from matplotlib.dates import AutoDateFormatter, AutoDateLocator, HOURLY
from obspy.core import Stream, read
# import pandas as pd
# import random
# # from pdart.csv_join_work_tapes import stream_import, initial_cleanup, merge_channel_stream, despike3, loose_frame_diff, read_file, process_list
# from pdart.csv_join_work_tapes import despike3, to_Int64, stream_import
# from pdart.csv_check_work_tapes import calculate_gaps, add_or_minus_frame
import pdart.config as config
import pdart.auth as auth
from obspy.clients.fdsn.client import Client
from pdart.util import linear_interpolation


def view_Apollo(stream=None,starttime= UTCDateTime('1973-03-13T07:30:00.0'),endtime = UTCDateTime('1973-03-13T09:30:00.0'),
  network='XA',station='S14',channel='MH1',location='*',plot_seismogram=True,plot_response=False):
    """Snippet to read in raw seismogram and remove the instrument response for Apollo.
    
    """
    
    user=auth.user
    auth_password=auth.auth_password

    if user == '' or auth_password == '':
        print('Set user and auth_password in auth.py')
        return
    
    client = Client("IRIS",user=user,password=auth_password)

    # get the response file (wildcards allowed)
    inv = client.get_stations(starttime=starttime, endtime=endtime,
        network=network, sta=station, loc=location, channel=channel,
        level="response")

    if stream is None:
        stream = client.get_waveforms(network=network, station=station, channel=channel, location=location, starttime=starttime, endtime=endtime)

    else:
        stream.trim(starttime=starttime,endtime=endtime)
        
    
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
    


# 1976-12-04T10:45:00 1976-12-04T10:50:00

def plot_data_spike():

    # print(UTCDateTime('1976-12-04T10:45:00').julday)
    # exit()

    stream = read('/Users/cnunn/lunar_data/PDART_CORRECTED/s12/1976/339/XA.S12.00.MH1.1976.339.0.MSEED')
    stream += read('/Users/cnunn/lunar_data/GEOSCOPE_lunar_data/1976/S12/MH1/S12.XA..MH1.1976.339')
    # stream.trim(starttime= UTCDateTime('1976-12-04T10:45:00'),endtime = UTCDateTime('1976-12-04T10:50:00'))
    stream.trim(starttime= UTCDateTime('1976-12-04T10:00:00'),endtime = UTCDateTime('1976-12-04T11:15:00'))
    for tr in stream:
        # interpolate across the gaps of one sample 
        linear_interpolation(tr,interpolation_limit=1)
    stream.merge()
    # stream.plot(equal_scale=False,size=(1000,600),method='full')

    fig = stream.plot(handle=True, show=False,size=(1200,600), method='full')
    for i, ax in enumerate(fig.get_axes()):
        if i == 0:
            pass
            # ax.set_ylim(-50,50)
        elif i == 1:
            ax.set_ylim(500,550)

    plt.show()

def plot_data_cleaning():
    # print(UTCDateTime('1976-03-04T08:09:59').julday)
    stream = read('/Users/cnunn/lunar_data/GEOSCOPE_lunar_data/1976/S12/MH1/S12.XA..MH1.1976.339')
    for tr in stream:
        tr.data = tr.data + 512
        tr.stats.starttime = tr.stats.starttime + 2.1


    stream += read('/Users/cnunn/lunar_data/PDART_CORRECTED/s12/1976/339/XA.S12.00.MH1.1976.339.0.MSEED')
    # stream += read('/Users/cnunn/lunar_data/GEOSCOPE_lunar_data/1976/S12/MH1/S12.XA..MH1.1976.64')
    # stream.trim(starttime= UTCDateTime('1976-03-04T08:09:59'),endtime = UTCDateTime('1976-03-04T14:19:59'))
    stream.trim(starttime= UTCDateTime('1976-12-04T10:00:00'),endtime = UTCDateTime('1976-12-04T11:15:00'))
    for tr in stream:
        # interpolate across the gaps of one sample 
        linear_interpolation(tr,interpolation_limit=1)
    stream.merge()
    # stream.plot(equal_scale=False,size=(1000,600),method='full')

    fig = stream.plot(handle=True, show=False,size=(1200,600), method='full')
    # for i, ax in enumerate(fig.get_axes()):
        # if i == 0:
        #     ax.set_ylim(500,550)
        # elif i == 1:
        #     ax.set_ylim(500,550)

    plt.show()


# def plot_data_cleaning():
#     # print(UTCDateTime('1976-03-04T08:09:59').julday)
#     stream = read('/Users/cnunn/lunar_data/GEOSCOPE_lunar_data/1976/S12/MH1/S12.XA..MH1.1976.64')
#     for tr in stream:
#         tr.data = tr.data + 512
# 
#     stream += read('/Users/cnunn/lunar_data/PDART_CORRECTED/s12/1976/064/XA.S12.00.MH1.1976.064.0.MSEED')
#     # stream += read('/Users/cnunn/lunar_data/GEOSCOPE_lunar_data/1976/S12/MH1/S12.XA..MH1.1976.64')
#     # stream.trim(starttime= UTCDateTime('1976-03-04T08:09:59'),endtime = UTCDateTime('1976-03-04T14:19:59'))
#     stream.trim(starttime= UTCDateTime('1976-12-04T10:00:00'),endtime = UTCDateTime('1976-12-04T11:15:00'))
#     for tr in stream:
#         # interpolate across the gaps of one sample 
#         linear_interpolation(tr,interpolation_limit=1)
#     stream.merge()
#     # stream.plot(equal_scale=False,size=(1000,600),method='full')
# 
#     fig = stream.plot(handle=True, show=False,size=(1200,600), method='full')
#     for i, ax in enumerate(fig.get_axes()):
#         print(i)
#         if i == 0:
#             ax.set_ylim(500,550)
#         elif i == 1:
#             ax.set_ylim(500,550)
# 
#     plt.show()

    # 
    # stream = view_Apollo(starttime= UTCDateTime('1976-03-04T08:09:59'),endtime = UTCDateTime('1976-03-04T14:19:59'),
    #   network='XA',station='S16',channel='MH1',location='*',plot_seismogram=True)
    # 



if __name__ == "__main__":
    plot_data_cleaning()
    # plot_data_spike()