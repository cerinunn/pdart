#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

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
import operator
import fnmatch
import shutil
import pandas as pd
from random import choice, seed
from itertools import groupby
import sys, traceback
from collections import OrderedDict
from pdart.csv_check_work_tapes import minus_frame, add_or_minus_frame, frame_diff

from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, Trace, Stats, read

# from pdart.save_24_hours import update_starttimes
import pdart.config as config
# from pdart.csv_import_work_tapes import find_output_dir, make_output_dir, make_filelist
# import matplotlib.pyplot as plt

# Qt5Agg seems to work best on Mac - try 'TkAgg' if that works for you
# put this after the other imports, otherwise it can be overridden
import matplotlib  
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

from pdart.util import maximize_plot

global total
total = 0
global gst_total
gst_total = 0

# global DELTA
DELTA = 0.1509433962
# slightly shorter delta, so that the timing doesn't run over 24 hours 
# SHORT_DELTA = 0.1509000000
INTERVAL = '603774us'

# INVALID = -99999
# INVALID = 0
INVALID = -1
# X should point north, Y east, but this is not always the case
# so we rename LPX to MH1, and LPY to MH2
NETWORK='XA'
TEST_SHIFT = 20 

# max percent error is based on 12 s in 24 hours (the maximum 
# time lag seen so far)
MAX_PERCENT_ERROR= 0.014

# global first_record
# global last_record

# consecutive frames
# FRAMES = list(range(0,90))


'''
# read in each of the csv files for the station

# read csv file and make a dataframe

# make a list of the dataframes

# process the list of dataframes with the process_list() method 
# it's this method which controls the start time of each of the streams,
# so this method is very important 

# sort the list by starttime

# sometimes there's a bit of overlap, but if we delete a couple of 
# samples, we can continue from the previous trace  

# sometimes, there's more than one record that could be the start record 

# import the stream for every record, using the same 
# sample_time0 and index0 using the stream_import() method 

# if required, clean the spikes (set config.clean_spikes = True)

# if required, merge the streams (set config.combine_ground_stations = True)

# shift the traces to reflect the different sampling time for each channel

# finally, save the streams with save_stream()
# merged streams will be saved in the following format s12/1976/229/xa.s12..mh1.1976.229.0.mseed
# unmerged streams will be saved in the following format   s12/1976/229/xa.s12.3.mh1.1976.068.0.mseed where the ground station is also included 
'''

def call_csv_join_work_tapes(
    processed_dir='.',
    join_dir='.',
    log_dir='.',
    wildcard_style='*',
    year_start=None,
    year_end=None,
    day_start=None,
    day_end=None,
    stations=['S11','S12','S14','S15','S16'],
    # logging_level=logging.DEBUG
    logging_level=logging.INFO,
    test=False, # these parameters are used for testing ONLY
    test_start=None, # these parameters are used for testing ONLY
    test_end=None, # these parameters are used for testing ONLY
    test_framejump=40, # these parameters are used for testing ONLY
    test_cutout_start=None,  # these parameters are used for testing ONLY
    test_cutout_end=None,  # these parameters are used for testing ONLY    
    ):

    '''
    Make files which are joined into the right days
    Calls csv_import_work_tapes()
    '''

    # TODO -fix these date ranges - they are not quite right 
    # if no filenames have been passed, then look for them 
    for year in range(year_start,year_end+1):
        for day in range(day_start,day_end+1):

            if config.combine_ground_stations == True:
                log_filename = 'join.{}.{}.log'.format(year,day)
            else:
                log_filename = 'ground_stations.{}.{}.log'.format(year,day)    
            log_filename = os.path.join(log_dir,log_filename)
            logging.basicConfig(filename=log_filename, filemode='w', 
              level=logging_level,format='%(message)s')

            print('Processing {}.{}'.format(year,day))

            logging.info('############################################')
            logging.info('Processing {}.{}'.format(year,day))
            logging.info('############################################')
    # /Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.6.S12.6.1976.063.00_00_00_199000.csv.gz
    # /Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.6.S12.601.1976.063.00_32_59_285000.csv.gz

            logging.info('#########')                                   
            for station in stations:
                df_list = []
                # wildcard_filename = 'pse.a16.1.83.S16.12.1972.195.00_00_00_159000*XXX.csv.gz'
                # wildcard_filename = 'pse.a15.2.168.S15.12.1972.195.00_00_00_144000.csv.gz'
                wildcard_filename = '{}.*.*.{}.*.{}.{:03}.*.csv.gz'.format(wildcard_style,station,year,day)
                print('wildcard filename ', processed_dir, wildcard_filename)
                # read in each of the csv files for the station
                for filename in glob.glob(os.path.join(processed_dir,wildcard_filename)):
                    # print(filename)
                    
                    # if filename not in ('/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.6.S12.6.1976.063.00_00_00_199000.csv.gz',
                    #     '/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.6.S12.601.1976.063.00_32_59_285000.csv.gz',
                    #     '/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.6.S12.9.1976.063.05_35_56_339000.csv.gz'):

                    # if filename not in (
                    #     '/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.19.S12.6.1976.069.05_08_48_484000.csv.gz',
                    #     '/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.11.21.S12.4.1976.069.06_20_12_532000.csv.gz'
                    #     ):
                    #     continue
                    # if filename in (
                    #     '/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.19.S12.6.1976.069.05_08_48_484000.csv.gz'
                    #     ):
                    #     continue
                    # if os.path.basename(filename) not in ['wtn.3.18.S16.5.1976.128.03_19_59_792000.csv.gz','wtn.3.18.S16.8.1976.128.04_08_31_629000.csv.gz']:
                    #     continue
                    if 'dropped' not in filename:
                        try: 
                            gzip_filename = filename
                            logging.info(gzip_filename)
                            df_list = read_file(filename,df_list,logging_level,join_dir,log_dir)
                        except Exception as e:
                            logging.info(traceback.format_exc())
                            logging.info(traceback.print_exc(file=sys.stdout))
                            # reset_config()
                            print('Warning, continuing')
                            # print('Not continuing')

                # XXXX
                logging.info('#########')
                if len(df_list) > 0:
                    # if test:
                    #     print('Testing dealing with problems with the trace. Start by moving the frame count by 40.')
                    #     for df_list1 in df_list:
                    #         df = df_list1[5]
                    #         if test_start is not None and test_end is not None: 
                    #             idx_list =  df[(df['corr_timestamp'] >= test_start) & (df['corr_timestamp'] <= test_end)].index.tolist()
                    #             # df['new_frame'] = df['frame']
                    #             for idx in idx_list:
                    #                 frame = df['frame'].iloc[idx]
                    #                 df.at[idx,'frame'] = add_or_minus_frame(frame,test_framejump)

                            # if test_cutout_start is not None and test_cutout_end is not None: 
                            #     idx_list =  df[(df['corr_timestamp'] >= test_cutout_start) & (df['corr_timestamp'] <= test_cutout_end)].index.tolist()
                            #     # df['new_frame'] = df['frame']
                            #     df.drop(df.index(to_drop)], inplace=True)
                            #     df.reset_index(inplace=True,drop=True)
                            #     for idx in idx_list:
                            #         frame = df['frame'].iloc[idx]
                            #         df.at[idx,'frame'] = add_or_minus_frame(frame,test_framejump)

                    
                    
                    try: 
                        # process the list of dataframes
                        process_list(df_list,join_dir)
                    except Exception as e:
                        logging.info(traceback.format_exc())
                        logging.info(traceback.print_exc(file=sys.stdout))
                        # reset_config()
                        print('Warning, error in process_list, continuing')
    

            # close the log file so that we can open a new one
            logger = logging.getLogger()
            handlers = logger.handlers[:]
            for handler in handlers:
                handler.close()
                logger.removeHandler(handler)

# def reset_config():
#     config.last_frame_S12=None
#     config.last_timestamp_S12=None
#     config.last_ground_station_S12=None
#     config.last_frame_S14=None
#     config.last_timestamp_S14=None
#     config.last_ground_station_S14=None
#     config.last_frame_S15=None
#     config.last_timestamp_S15=None
#     config.last_ground_station_S15=None
#     config.last_frame_S16=None
#     config.last_timestamp_S16=None
#     config.last_ground_station_S16=None

def read_file(gzip_filename,df_list,logging_level=logging.INFO,join_dir='.',log_dir='.',
  df=None # only used for test purposes
):

    # unless we are testing, df is None, so read in the file
    if df is None:
        # read csv file and make a dataframe
        df = pd.read_csv(gzip_filename, dtype=str, comment='#')

    if len(df) < 3:
        logging.debug('WARNING: Only {} record(s), not importing'.format(len(df)))
        return df_list
    else:
        logging.debug('{} record(s)'.format(len(df)))

    df = initial_cleanup(df)

    # get some details about the data 
    starttime = UTCDateTime(df.corr_timestamp.iloc[0])
    endtime = UTCDateTime(df.corr_timestamp.iloc[-1])
    orig_station = df.orig_station.iloc[0]
    orig_ground_station = df.orig_ground_station.iloc[0]
    corr_ground_station = df.corr_ground_station.iloc[0]
    attempt_merge = True
    
    # make a list of the dataframes
    df_list.append([starttime, endtime, orig_station,orig_ground_station, corr_ground_station, df, attempt_merge])

    return df_list

def process_list(df_list,join_dir):

    # for each station, find the first anchor
    # run through the rest of the dataframes for the station


    orig_station = df_list[0][2]
    # logging.info('This is the orig station {}'.format(orig_station))

    if orig_station == 'S12':
        # last_frame_join = config.last_frame_S12
        last_timestamp_join = config.last_timestamp_S12
        # last_ground_station_join = config.last_ground_station_S12
        # mean_S12_MH1 = config.mean_S12_MH1
        # mean_S12_MH2 = config.mean_S12_MH2
        # mean_S12_MHZ = config.mean_S12_MHZ
    elif orig_station == 'S14':
        # last_frame_join = config.last_frame_S14
        last_timestamp_join = config.last_timestamp_S14
        # last_ground_station_join = config.last_ground_station_S14
        # mean_join_MH1 = config.mean_S14_MH1
        # mean_join_MH2 = config.mean_S14_MH2
        # mean_join_MHZ = config.mean_S14_MHZ
    elif orig_station == 'S15':
        # last_frame_join = config.last_frame_S15
        last_timestamp_join = config.last_timestamp_S15
        # last_ground_station_join = config.last_ground_station_S15
        # mean_join_MH1 = config.mean_S15_MH1
        # mean_join_MH2 = config.mean_S15_MH2
        # mean_join_MHZ = config.mean_S15_MHZ
    elif orig_station == 'S16':
        # last_frame_join = config.last_frame_S16
        last_timestamp_join = config.last_timestamp_S16
        # last_ground_station_join = config.last_ground_station_S16
        # mean_join_MH1 = config.mean_S16_MH1
        # mean_join_MH2 = config.mean_S16_MH2
        # mean_join_MHZ = config.mean_S16_MHZ

    stream = Stream()
    # stations = ['S12','S14','S15','S16']
    # 
    # # split by station
    # for station in stations:
    #     df_list_station = []
    #     for df_list1 in df_list:
    #         orig_station = df_list1[2]
    #         if orig_station == station:
    #             df_list_station.append(df_list1)

    # sort by starttime
    sorted_df_list = sorted(df_list, key=operator.itemgetter(0))

    continue_idx = None
    # if last_ground_station_join is not None:
    #     # make sure we are continuing with the record 
    #     # only perform this check if the endtime was in the last second
    #     # before midnight
    #     last_timestamp_join = UTCDateTime(last_timestamp_join)
    #     day = last_timestamp_join.day
    #     if (last_timestamp_join + 1).day != day:
    #         for i, df_list1 in enumerate(sorted_df_list):
    #             starttime = df_list1[0]
    #             corr_ground_station = df_list1[4]
    #             if ((corr_ground_station == last_ground_station_join) and 
    #                 (starttime < (last_timestamp_join + 1))):
    #                 continue_idx = i
    #                 logging.info('INFO: Continuing with trace from ground station {} index={}'.format(corr_ground_station,i))
    #                 break

    if continue_idx is not None: 
        if continue_idx != 0:
            # logging.info('EXCEPTION: Raise Exception to test this strange case.')
            # raise Exception

            # sometimes there's a bit of overlap, but if we delete a couple of 
            # samples, we can continue from the previous trace  
            for i, df_list1 in enumerate(sorted_df_list):
                if i < continue_idx:
                    starttime = df_list1[0]

                    logging.info('WARNING: Deleting two records from trace {}, to continue from trace {}'.format(i, continue_idx))

                    # these starttimes are close, so delete the first couple of 
                    # records in this trace 
                    df = df_list1[5]
                    df = df.iloc[2:]
                    # reindex
                    df.reset_index(inplace=True,drop=True)
                    new_starttime = UTCDateTime(df.corr_timestamp.iloc[0])
                    sorted_df_list[i][0] = new_starttime
                    sorted_df_list[i][5] = df
                    
            # re-sort by starttime
            sorted_df_list = sorted(df_list, key=operator.itemgetter(0))
            # now it's been re-sorted, the continuing index will be first 
            logging.info('WARNING: Traces resorted')
            continue_idx = 0

    # sometimes, there's more than one record that could be the start record 
    for i, df_list1 in enumerate(sorted_df_list):
        if i == 0:
            starttime_i0 = df_list1[0]
        if i > 0:
            starttime = df_list1[0]
            if starttime < (starttime_i0+ 1):
                # these starttimes are close, so delete the first couple of 
                # records in this trace 
                df = df_list1[5]
                df = df.iloc[2:]
                logging.info('WARNING: Deleting two records from trace {}, to continue from initial trace {}'.format(i, 0))
                # reindex
                df.reset_index(inplace=True,drop=True)
                sorted_df_list[i][0] = df.corr_timestamp.iloc[0]
                sorted_df_list[i][5] = df
            else:
                break

    starttime0 = None
    index0 = None
    for i, df_list1 in enumerate(sorted_df_list):
        df = df_list1[5]
        # logging.info(i)
        # logging.info(df.head().to_string())
        # logging.info(df.tail().to_string())

        # for each station, set up the time_index
        if i==0:
            starttime0, sample_time0 = first_record(index_no=i,sorted_df_list=sorted_df_list)
        elif i > 0:
            later_records(index_no=i,sample_time0=sample_time0,sorted_df_list=sorted_df_list)
                
        df_list1 = sorted_df_list[i]
        df = df_list1[5]
        attempt_merge = df_list1[6]
        # logging.info('attempt_merge {} {}'.format(i, attempt_merge))

        gaps_7777 = (df['corr_gap_count'] == -7777).sum()
        if gaps_7777 > 0:
            orig_station = df_list1[2]
            orig_ground_station = df_list1[3]
    # df_list.append([starttime, endtime, orig_station,orig_ground_station, corr_ground_station, df])
            logging.info('# WARNING: Ground station/Station - {} {} {} frames are reset to zero (-7777 error)'.format(orig_ground_station,orig_station,gaps_7777))

        # 
        # logging.info('df?')
        # logging.info(df.head().to_string())
        # logging.info(df.dtypes.to_string())


        # import the stream for every record, using the same 
        # time_index0 and index0
        stream1 = stream_import(df,sample_time0,index0,attempt_merge)

        stream += stream1

    # if required, clean the spikes 
    if config.clean_spikes:
        for tr in stream:
            despike3(tr)

    # if required, merge the streams 
    if config.combine_ground_stations: 
        merged_stream = merge_streams(stream)
        stream = merged_stream

    # shift the traces to reflect the different sampling time for each channel
    timeshift_traces(stream)
    
    logging.debug(stream.__str__(extended=True))

    save_stream(stream,join_dir)
    
    # channel_stream = stream.select(channel='MHZ')
    # fig = channel_stream.plot(handle=True, show=False,size=(1200,600), method='full')
    # plt.show()
    # returns the stream - for test purposes only 
    return stream 

def merge_channel_stream(channel_stream,delta):

    logging.debug('merge_channel_stream')

    # channel_stream = channel_stream.merge(fill_value=INVALID)
    # this contains merged traces at this point 


    channels = set()
    for tr in channel_stream:
        channels.add (tr.stats.channel)

    if len(channels) > 1:
        print('EXCEPTION: merge_channel_stream() requires a single channel')
        logging.info('EXCEPTION: merge_channel_stream() requires a single channel')
        raise Exception

    # sort the stream
    channel_stream.sort(keys=['starttime', 'endtime'])


    # if there's only one stream, return it, because it doesn't require merging
    if len(channel_stream) < 2:
        return channel_stream

    new_channel_stream = Stream()


    for i, tr in enumerate(channel_stream.sort()):
        if i == 0:
            starttime0 = channel_stream[0].stats.starttime
            time_index0 = 0
            # make a basic trace 
            tr_new = channel_stream[0].copy()

            prev_tr = tr
            prev_time_index = np.arange(len(prev_tr.data))
            prev_d = {'data': prev_tr.data, 'time_index':prev_time_index }
            prev_df = pd.DataFrame(data=prev_d)
            prev_df['data'] = to_Int64_with_invalid(prev_df['data'])

        if i > 0:

            starttime = tr.stats.starttime
            time_int = round((starttime - starttime0)/delta)
            time_index = np.arange(len(tr.data)) + time_int
            d = {'data': tr.data, 'time_index':time_index}
            df = pd.DataFrame(data=d)
            df['data'] = to_Int64_with_invalid(df['data'])

            # merge the two data fields
            df_result = pd.merge(prev_df,df, how='outer', on='time_index')
            # logging.info(df_result.to_string())

            attempt_merge = tr.stats.attempt_merge

            if attempt_merge:

                # characterize the overlap
                count_diff =  len(df_result[~(df_result['data_x'].isna()) & ~(df_result['data_y'].isna()) & (df_result['data_y'] !=  df_result['data_x'])].index.tolist())
                count_overlap =  len(df_result[(~(df_result['data_x'].isna()) & ~(df_result['data_y'].isna()))].index.tolist())

                if count_overlap > 0:
                    percent_diff = 100*count_diff/count_overlap
                    if percent_diff > 5:
                        # if it's not overlapping properly, don't try to merge
                        attempt_merge = False
                        tr.stats.attempt_merge = False  
                        logging.info('INFO: Unable to merge because {:.1f}% of the records are different. {} {}'.format(percent_diff, tr.id, tr.stats.starttime))
                    
            df_result['data_x'] = to_Int64_with_invalid(df_result['data_x'])
            df_result['data_y'] = to_Int64_with_invalid(df_result['data_y'])        


            # if too many problems are found, no merge is made
            # instead, it just uses the later trace 
            if attempt_merge == False:
                logging.info('INFO: Adding trace without attempting a merge. {} {}'.format(tr.id, tr.stats.starttime))

                if config.view_corrected_traces:
                    stream_view_error = Stream()
                    # to view the error (later) if required
                    tr_data_x = channel_stream[0].copy()
                    tr_data_x.data = df_result['data_x'].to_numpy(dtype=np.int32, na_value=INVALID)
                    if INVALID in tr_data_x.data:
                        tr_data_x.data = ma.masked_where(tr_data_x.data == INVALID, tr_data_x.data)
                    stream_view_error.append(tr_data_x)
                    tr_data_y = channel_stream[0].copy()
                    tr_data_y.data = df_result['data_y'].to_numpy(dtype=np.int32, na_value=INVALID)
                    if INVALID in tr_data_y.data:
                        tr_data_y.data = ma.masked_where(tr_data_y.data == INVALID, tr_data_y.data)
                    stream_view_error.append(tr_data_y)

                # take the new data for the second part of the trace 
                df_result['data'] = np.where( (df_result['time_index'] >= time_int), 
                    df_result['data_y'], df_result['data_x'])



                if config.view_corrected_traces:
                    # view the error, if required
                    tr_data = channel_stream[0].copy()
                    df_result['data'] = to_Int64_with_invalid(df_result['data'])
                    tr_data.data = df_result['data'].to_numpy(dtype=np.int32, na_value=INVALID)

                    # either replace INVALID with gaps, to view (slow)
                    # if INVALID in tr_data.data:
                    #     tr_data.data = ma.masked_where(tr_data.data == INVALID, tr_data.data)
                    # or replace with a number near the middle
                    if INVALID in tr_data.data:
                        tr_data.data = np.where(tr_data.data == INVALID, 500, tr_data.data)
                    
                    # merge to a single record??

                    # how to view ??? 
                    
                    stream_view_error.append(tr_data)
                    # stream_view_error = stream_view_error.trim(starttime=stream_view_error[0].stats.starttime+0,endtime=stream_view_error[0].stats.starttime+300)
                    # stream_view_error = stream_view_error.trim(starttime=UTCDateTime('1977-08-20T21:25:42.441792Z')-30,endtime=UTCDateTime('1977-08-20T21:25:42.441792Z')+300)
        

                    stream_view_error.plot(method='full',size=(1200,600),automerge=True,handle=True,show=False)
                    maximize_plot()
                    plt.show()


            else: 

                # if tr.id == 'XA.S12..MH2':


                # if config.view_corrected_traces:
                #     stream_view_error = Stream()
                #     # to view the error (later) if required
                #     tr_data_x = channel_stream[0].copy()
                #     tr_data_x.data = df_result['data_x'].to_numpy(dtype=np.int32, na_value=INVALID)
                #     if INVALID in tr_data_x.data:
                #         tr_data_x.data = ma.masked_where(tr_data_x.data == INVALID, tr_data_x.data)
                #     stream_view_error.append(tr_data_x)
                #     tr_data_y = channel_stream[0].copy()
                #     tr_data_y.data = df_result['data_y'].to_numpy(dtype=np.int32, na_value=INVALID)
                #     if INVALID in tr_data_y.data:
                #         tr_data_y.data = ma.masked_where(tr_data_y.data == INVALID, tr_data_y.data)
                #     stream_view_error.append(tr_data_y)

                # if the first stream has nulls, replace with non nulls from the second stream
                df_result.loc[df_result['data_x'].isnull(),'data_x'] = df_result['data_y']
                # new column with all the data in 
                df_result['data'] = df_result['data_x']

                # find the mean
                # mean_x = df_result.data_x.mean()

                # logging.info(df_result.to_string())

                # idx_list =  df_result[~(df_result['data_x'].isna()) & ~(df_result['data_y'].isna()) & (df_result['data_y'] !=  df_result['data_x'])].index.tolist()
                # df_fil = df_result[~df_result['data_y'].isna()]

                # #################################
                # disply this if there are problems 
                # df_result['diff'] = np.where( (~(df_result['data_x'].isna()) & ~(df_result['data_y'].isna()) & ((df_result['data_y']) != df_result['data_x'])), 
                #     'yes_diff', 'no_diff')
                # logging.debug(df_result.to_string())

                


                # idx_list = df_result[~df_result['data_y'].isna()].index.tolist()
                # logging.info(df_fil.to_string())
                # idx_list =  ((~df_result['data_y'].isna())).index.tolist()
                # logging.info(df_result.iloc[idx_list].to_string())


                # TODO maybe make some counts
                # count the number of different records when both records are not
                # null
                # count1 = len(idx_list)


                # if tr.id == ''
                # logging.info(df_result.to_string())
                # logging.info(tr.id)
                # if tr.id == 'XA.S12..MH2':
                #     df_count = df_result.copy()
                #     df_count.dropna(inplace=True)
                #     count2 = (df_count['data_x'] != df_count['data_y']).sum()
                #     df_count3 = df_count[df_count['data_x'] != df_count['data_y']]              
                #     logging.info('This')
                #     logging.info(count1)
                #     logging.info(count2)
                #     logging.info(df_count3.to_string())


                # get a total count of errors 
                # df_count = df_result.copy()
                # df_count.dropna(inplace=True)
                # df_count['count'] = df_count[df_count['data_x'] != df_count['data_y']]
                    


                # copy the left side 

                # replace the value on the left with the value on the right,
                # if both values are not null and the value on the right is
                # closer to the mean 
                # df_result['data'] = np.where( (~(df_result['data_x'].isna()) & ~(df_result['data_y'].isna()) & ((abs(df_result['data_y'] - mean_x)) < abs(df_result['data_x'] - mean_x))), 
                #      df_result['data_y'], df_result['data_x'])

                # logging.info(df_result.to_string())
                # df_result['data'] = np.where( (df_result['data_x'].isna()) & ~(df_result['data_y'].isna()) , 
                #     df_result['data'], df_result['data_y'])

                # logging.info('stuff')
                # logging.info(df_result.iloc[4823:4824].to_string())
                # logging.info('stuff2')
                # logging.info(df_result.iloc[idx_list[0:10]].to_string())

                # logging.info(df_result.to_string())

                # count the number of different records between the data column and the data_x column when both records are not
                # null
                # count1 = np.count_nonzero( ~(df_result['data'].isna()) & ~(df_result['data_x'].isna()) & (df_result['data'] !=  df_result['data_x']))

            df_result.set_index('time_index', inplace=True)
            df_result = df_result.reindex(index = np.arange(df_result.index[0], df_result.index[-1] + 1), fill_value=pd.NA).reset_index()

            # make a new prev_df 
            prev_df = df_result[['data', 'time_index']].copy()


            # if attempt_merge:
            #     total = len(df)
            #     percent_error = 100 * count1/total
            # 
            #     if percent_error > 1: 
            #         print('SEVERE: Corrected {} data record(s) of {}({})'.format(count1,total,tr.id))
            #         error_starttime = starttime0 + idx_list[0]*delta
            #         error_endtime = starttime0 + idx_list[-1]*delta
            #         # my_error = starttime0 + 148620*delta
            # 
            #         # logging.info(str(tr))
            #         # logging.info('SEVERE: my time {}'.format(my_error))
            #         logging.info('SEVERE: Corrected {} data record(s) of {} ({}) from {} to {}'.format(count1,total,tr.id,error_starttime,error_endtime))
            #         print('SEVERE: Corrected {} data record(s) of {} ({}) from {} to {}'.format(count1,total,tr.id,error_starttime,error_endtime))
            #         if config.view_corrected_traces:
            #             # view the error, if required
            #             tr_data = channel_stream[0].copy()
            #             df_result['data'] = to_Int64_with_invalid(df_result['data'])
            #             tr_data.data = df_result['data'].to_numpy(dtype=np.int32, na_value=INVALID)
            # 
            #             # either replace INVALID with gaps, to view (slow)
            #             # if INVALID in tr_data.data:
            #             #     tr_data.data = ma.masked_where(tr_data.data == INVALID, tr_data.data)
            #             # or replace with a number near the middle
            #             if INVALID in tr_data.data:
            #                 tr_data.data = np.where(tr_data.data == INVALID, 500, tr_data.data)
            # 
            #             # merge to a single record??
            # 
            #             # how to view ??? 
            # 
            #             stream_view_error.append(tr_data)
            #             # stream_view_error = stream_view_error.trim(starttime=stream_view_error[0].stats.starttime+0,endtime=stream_view_error[0].stats.starttime+300)
            #             # stream_view_error = stream_view_error.trim(starttime=UTCDateTime('1977-08-20T21:25:42.441792Z')-30,endtime=UTCDateTime('1977-08-20T21:25:42.441792Z')+300)
            # 
            # 
            #             stream_view_error.plot(method='full',size=(1200,600),automerge=True,handle=True,show=False)
            #             maximize_plot()
            #             plt.show()
            #     elif count1 > 0: 
            #         logging.info('WARNING: Corrected {} data record(s) of {} ({})'.format(count1,total,tr.id))
            #         logging.info(starttime)

            
    df_result['data'] = to_Int64(df_result['data'])
    
    tr_new.data = df_result['data'].to_numpy(dtype=np.int32, na_value=INVALID)
    
    # if INVALID in tr_new.data:
    #     tr_new.data = ma.masked_where(tr_new.data == INVALID, tr_new.data)
    # 
    # check_trace_mean(tr_new)
    # put the single final trace into the new stream 

    new_channel_stream.append(tr_new)
    # for tr in new_channel_stream:
    #     logging.info('v end {}'.format(len(tr.data)))

    return new_channel_stream

def check_trace_mean(trace):

    station = trace.stats.station
    channel = trace.stats.channel
    if station == 'S12':
        if channel=='MH1':
            trace_mean = config.mean_S12_MH1
        elif channel=='MH2':
            trace_mean = config.mean_S12_MH2
        elif channel=='MHZ':
            trace_mean = config.mean_S12_MHZ
        elif channel=='SHZ':
            trace_mean = config.mean_S12_SHZ

    elif station == 'S14':
        if channel=='MH1':
            trace_mean = config.mean_S14_MH1
        elif channel=='MH2':
            trace_mean = config.mean_S14_MH2
        elif channel=='MHZ':
            trace_mean = config.mean_S14_MHZ
        elif channel=='SHZ':
            trace_mean = config.mean_S15_SHZ

    elif station == 'S15':
        if channel=='MH1':
            trace_mean = config.mean_S15_MH1
        elif channel=='MH2':
            trace_mean = config.mean_S15_MH2
        elif channel=='MHZ':
            trace_mean = config.mean_S15_MHZ
        elif channel=='SHZ':
            trace_mean = config.mean_S15_SHZ

    elif station == 'S16':
        if channel=='MH1':
            trace_mean = config.mean_S16_MH1
        elif channel=='MH2':
            trace_mean = config.mean_S16_MH2
        elif channel=='MHZ':
            trace_mean = config.mean_S16_MHZ
        elif channel=='SHZ':
            trace_mean = config.mean_S16_SHZ

    new_trace_mean = trace.data.mean()
    if trace_mean is not None: 
        if abs(new_trace_mean-trace_mean) > 3:  
                logging.info('WARNING: The mean of the trace was {:.2f} and is now {:.2f}'.format(trace_mean,new_trace_mean))     
                print('WARNING: The mean of the trace was {:.2f} and is now {:.2f}'.format(trace_mean,new_trace_mean))   
  
    if station == 'S12':
        if channel=='MH1':
            config.mean_S12_MH1 = new_trace_mean
        elif channel=='MH2':
            config.mean_S12_MH2 = new_trace_mean
        elif channel=='MHZ':
            config.mean_S12_MHZ = new_trace_mean
        elif channel=='SHZ':
            config.mean_S12_SHZ = new_trace_mean
    elif station == 'S14':
        if channel=='MH1':
            config.mean_S14_MH1 = new_trace_mean
        elif channel=='MH2':
            config.mean_S14_MH2 = new_trace_mean
        elif channel=='MHZ':
            config.mean_S14_MHZ = new_trace_mean
        elif channel=='SHZ':
            config.mean_S14_SHZ = new_trace_mean
    elif station == 'S15':
        if channel=='MH1':
            config.mean_S15_MH1 = new_trace_mean
        elif channel=='MH2':
            config.mean_S15_MH2 = new_trace_mean
        elif channel=='MHZ':
            config.mean_S15_MHZ = new_trace_mean
        elif channel=='SHZ':
            config.mean_S15_SHZ = new_trace_mean
    elif station == 'S16':
        if channel=='MH1':
            config.mean_S16_MH1 = new_trace_mean
        elif channel=='MH2':
            config.mean_S16_MH2 = new_trace_mean
        elif channel=='MHZ':
            config.mean_S16_MHZ = new_trace_mean
        elif channel=='SHZ':
            config.mean_S16_SHZ = new_trace_mean
        

def merge_streams(stream):

    merged_stream = Stream()

    for channel in ['SHZ']:
        channel_stream = stream.select(channel=channel)
        new_channel_stream = merge_channel_stream(channel_stream,delta=DELTA/8)
        merged_stream += new_channel_stream

    for channel in ['MH1', 'MH2', 'MHZ']:
        channel_stream = stream.select(channel=channel)
        new_channel_stream = merge_channel_stream(channel_stream,delta=DELTA)
        merged_stream += new_channel_stream

    for channel in ['ATT']:
        channel_stream = stream.select(channel=channel)
        channel_stream.merge(method=1)
        merged_stream += channel_stream

    if config.combine_ground_stations == False:
        for channel in ['AFR']:
            channel_stream = stream.select(channel=channel)
            channel_stream.merge(method=1)

        merged_stream += channel_stream

    # record the last entry
    station = merged_stream[0].stats.station
    if station == 'S12':
        # config.last_frame_S12=merged_stream.select(channel='AFR')[0].data[-1]
        config.last_timestamp_S12=merged_stream.select(channel='ATT')[0].data[-1]
        # config.last_ground_station_S12=merged_stream.select(channel='AGR')[0].data[-1]
        end_time = merged_stream.select(channel='ATT')[0].stats.endtime
        config.time_divergence_S12 = end_time - UTCDateTime(config.last_timestamp_S12)

        # logging.info('config.last_frame_S12={}'.format(config.last_frame_S12))
        logging.info('config.last_timestamp_S12={}'.format(config.last_timestamp_S12))
        # logging.info('config.last_ground_station_S12={}'.format(config.last_ground_station_S12))
        logging.info('config.time_divergence_S12={}'.format(config.time_divergence_S12))

        if abs(config.time_divergence_S12) > 20:
            logging.info('SEVERE: Time divergence is {} for S12'.format(config.time_divergence_S12))
        # 
        # try: 
        #     config.last_frame_S12_MH1=merged_stream.select(channel='MH1')[0].data.mean()
        # logging.info('config.mean_S12_MH1={:.2f}'.format(config.mean_S12_MH1))
        # logging.info('config.mean_S12_MH2={:.2f}'.format(config.mean_S12_MH2))
        # logging.info('config.mean_S12_MHZ={:.2f}'.format(config.mean_S12_MHZ))
        # logging.info('config.mean_S12_SHZ={:.2f}'.format(config.mean_S12_SHZ))

        # # except:
        # #     logging.info('#No MH1 for S12'.format(config.last_frame_S12_MH1))
        # 
        # try: 
        #     config.last_frame_S12_MH2=merged_stream.select(channel='MH2')[0].data.mean()
        #     logging.info('config.last_frame_S12_MH2={:.2f}'.format(config.last_frame_S12_MH2))
        # except:
        #     logging.info('#No MH2 for S12'.format(config.last_frame_S12_MH2))
        # 
        # try: 
        #     config.last_frame_S12_MHZ=merged_stream.select(channel='MHZ')[0].data.mean()
        #     logging.info('config.last_frame_S12_MH1={:.2f}'.format(config.last_frame_S12_MHZ))
        # except:
        #     logging.info('#No MHZ for S12'.format(config.last_frame_S12_MHZ))


    elif station == 'S14':
        # config.last_frame_S14=merged_stream.select(channel='AFR')[0].data[-1]
        config.last_timestamp_S14=merged_stream.select(channel='ATT')[0].data[-1]
        # config.last_ground_station_S14=merged_stream.select(channel='AGR')[0].data[-1]
        end_time = merged_stream.select(channel='ATT')[0].stats.endtime
        config.time_divergence_S14 = end_time - UTCDateTime(config.last_timestamp_S14)

        # logging.info('config.last_frame_S14={}'.format(config.last_frame_S14))
        logging.info('config.last_timestamp_S14={}'.format(config.last_timestamp_S14))
        # logging.info('config.last_ground_station_S14={}'.format(config.last_ground_station_S14))
        logging.info('config.time_divergence_S14={}'.format(config.time_divergence_S14))

        if abs(config.time_divergence_S14) > 20:
            logging.info('SEVERE: Time divergence is {} for S14'.format(config.time_divergence_S14))

        # try: 
        #     config.last_frame_S14_MH1=merged_stream.select(channel='MH1')[0].data.mean()
        #     logging.info('config.last_frame_S14_MH1={:.2f}'.format(config.last_frame_S14_MH1))
        # except:
        #     logging.info('#No MH1 for S14'.format(config.last_frame_S14_MH1))
        # 
        # try: 
        #     config.last_frame_S14_MH2=merged_stream.select(channel='MH2')[0].data.mean()
        #     logging.info('config.last_frame_S14_MH2={:.2f}'.format(config.last_frame_S14_MH2))
        # except:
        #     logging.info('#No MH2 for S14'.format(config.last_frame_S14_MH2))
        # 
        # try: 
        #     config.last_frame_S14_MHZ=merged_stream.select(channel='MHZ')[0].data.mean()
        #     logging.info('config.last_frame_S14_MH1={:.2f}'.format(config.last_frame_S14_MHZ))
        # except:
        #     logging.info('#No MHZ for S14'.format(config.last_frame_S14_MHZ))

        # logging.info('config.mean_S14_MH1={:.2f}'.format(config.mean_S14_MH1))
        # logging.info('config.mean_S14_MH2={:.2f}'.format(config.mean_S14_MH2))
        # logging.info('config.mean_S14_MHZ={:.2f}'.format(config.mean_S14_MHZ))
        # logging.info('config.mean_S14_SHZ={:.2f}'.format(config.mean_S14_SHZ))

    elif station == 'S15':
        # config.last_frame_S15=merged_stream.select(channel='AFR')[0].data[-1]
        config.last_timestamp_S15=merged_stream.select(channel='ATT')[0].data[-1]
        # config.last_ground_station_S15=merged_stream.select(channel='AGR')[0].data[-1]
        end_time = merged_stream.select(channel='ATT')[0].stats.endtime
        config.time_divergence_S15 = end_time - UTCDateTime(config.last_timestamp_S15)

        # logging.info('config.last_frame_S15={}'.format(config.last_frame_S15))
        logging.info('config.last_timestamp_S15={}'.format(config.last_timestamp_S15))
        # logging.info('config.last_ground_station_S15={}'.format(config.last_ground_station_S15))
        logging.info('config.time_divergence_S15={}'.format(config.time_divergence_S15))

        if abs(config.time_divergence_S15) > 20:
            logging.info('SEVERE: Time divergence is {} for S15'.format(config.time_divergence_S15))


        # logging.info('config.mean_S15_MH1={:.2f}'.format(config.mean_S15_MH1))
        # logging.info('config.mean_S15_MH2={:.2f}'.format(config.mean_S15_MH2))
        # logging.info('config.mean_S15_MHZ={:.2f}'.format(config.mean_S15_MHZ))
        # logging.info('config.mean_S15_SHZ={:.2f}'.format(config.mean_S15_SHZ))


        # try: 
        #     config.last_frame_S15_MH1=merged_stream.select(channel='MH1')[0].data.mean()
        #     logging.info('config.last_frame_S15_MH1={:.2f}'.format(config.last_frame_S15_MH1))
        # except:
        #     logging.info('#No MH1 for S15'.format(config.last_frame_S15_MH1))
        # 
        # try: 
        #     config.last_frame_S15_MH2=merged_stream.select(channel='MH2')[0].data.mean()
        #     logging.info('config.last_frame_S15_MH2={:.2f}'.format(config.last_frame_S15_MH2))
        # except:
        #     logging.info('#No MH2 for S15'.format(config.last_frame_S15_MH2))
        # 
        # try: 
        #     config.last_frame_S15_MHZ=merged_stream.select(channel='MHZ')[0].data.mean()
        #     logging.info('config.last_frame_S15_MH1={:.2f}'.format(config.last_frame_S15_MHZ))
        # except:
        #     logging.info('#No MHZ for S15'.format(config.last_frame_S15_MHZ))


    elif station == 'S16':
        # config.last_frame_S16=merged_stream.select(channel='AFR')[0].data[-1]
        config.last_timestamp_S16=merged_stream.select(channel='ATT')[0].data[-1]
        # config.last_ground_station_S16=merged_stream.select(channel='AGR')[0].data[-1]
        end_time = merged_stream.select(channel='ATT')[0].stats.endtime
        config.time_divergence_S16 = end_time - UTCDateTime(config.last_timestamp_S16)

        # logging.info('config.last_frame_S16={}'.format(config.last_frame_S16))
        logging.info('config.last_timestamp_S16={}'.format(config.last_timestamp_S16))
        # logging.info('config.last_ground_station_S16={}'.format(config.last_ground_station_S16))
        logging.info('config.time_divergence_S16={}'.format(config.time_divergence_S16))

        if abs(config.time_divergence_S16) > 20:
            logging.info('SEVERE: Time divergence is {} for S16'.format(config.time_divergence_S16))
        # 
        # logging.info('config.mean_S16_MH1={:.2f}'.format(config.mean_S16_MH1))
        # logging.info('config.mean_S16_MH2={:.2f}'.format(config.mean_S16_MH2))
        # logging.info('config.mean_S16_MHZ={:.2f}'.format(config.mean_S16_MHZ))
        # logging.info('config.mean_S16_SHZ={:.2f}'.format(config.mean_S16_SHZ))

        # try: 
        #     config.last_frame_S16_MH1=merged_stream.select(channel='MH1')[0].data.mean()
        #     logging.info('config.last_frame_S16_MH1={:.2f}'.format(config.last_frame_S16_MH1))
        # except:
        #     logging.info('#No MH1 for S16'.format(config.last_frame_S16_MH1))
        # 
        # try: 
        #     config.last_frame_S16_MH2=merged_stream.select(channel='MH2')[0].data.mean()
        #     logging.info('config.last_frame_S16_MH2={:.2f}'.format(config.last_frame_S16_MH2))
        # except:
        #     logging.info('#No MH2 for S16'.format(config.last_frame_S16_MH2))
        # 
        # try: 
        #     config.last_frame_S16_MHZ=merged_stream.select(channel='MHZ')[0].data.mean()
        #     logging.info('config.last_frame_S16_MH1={:.2f}'.format(config.last_frame_S16_MHZ))
        # except:
        #     logging.info('#No MHZ for S16'.format(config.last_frame_S16_MHZ))

    return merged_stream
    

# The old code kept the timing together until the end.
# This is probably a good idea 
def timeshift_traces(stream):
    '''
    We provide the MH1, MH2, MHZ and SHZ traces with a time-shift of 0.075, 0.094, 
    0.113 and 0.009 s
    relative to the ATT and AFR traces.
    '''

    for tr in stream:
        if tr.stats.channel == 'MH1':
            tr.stats.starttime = tr.stats.starttime + 0.075
        elif tr.stats.channel == 'MH2':
            tr.stats.starttime = tr.stats.starttime + 0.094
        elif tr.stats.channel == 'MHZ':
            tr.stats.starttime = tr.stats.starttime + 0.113
        elif tr.stats.channel == 'SHZ':
            tr.stats.starttime = tr.stats.starttime + 0.009


def first_record(index_no,sorted_df_list):

    df_list1 = sorted_df_list[index_no]
    starttime0 = df_list1[0]
    endtime = df_list1[1]
    orig_station = df_list1[2]
    # orig_ground_station = df_list1[3]
    corr_ground_station = df_list1[4]
    df = df_list1[5]

    # create a hash on the combined columns
    df['hash'] = pd.util.hash_pandas_object(df[[
        'orig_mh1_1','orig_mh2_1','orig_mhz_1',
        'orig_mh1_2','orig_mh2_2','orig_mhz_2',
        'orig_mh1_3','orig_mh2_3','orig_mhz_3',
        'orig_mh1_4','orig_mh2_4','orig_mhz_4',
        'frame']], 
        index=False)

    # logging.info('{} starttime0={} 1st frame={} endtime={} Last frame={} ground_station={}'.format(index_no, starttime0, df.frame.iloc[0], endtime, df.frame.iloc[-1],corr_ground_station))


    # first check that the starttime0 is near midnight
    midnight = UTCDateTime(year=starttime0.year,julday=starttime0.julday)
    test_time = midnight + 0.6091
    if starttime0 < test_time:
        df.time_index = np.arange(len(df))
        sample_time0 = starttime0
    else:
        logging.info('WARNING: Ground station: {} Station: {} Starttime not immediately after midnight {}'.format(corr_ground_station,orig_station,starttime0))
        delta4 = df.delta4.iloc[0]
        
        time_int = round((starttime0 - midnight)/delta4)

        # start the index at zero 
        df.time_index = np.arange(len(df))
        
        sample_time0 = midnight + time_int *DELTA*4

    logging.info('INFO: Zero index record. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={}'.format(index_no, orig_station, corr_ground_station, starttime0, sample_time0, starttime0-sample_time0))

    return starttime0, sample_time0

def later_records(index_no,sample_time0, sorted_df_list):

    # get the information from the sorted_df_list
    df_list1 = sorted_df_list[index_no]
    # starttime = df_list1[0]
    # endtime = df_list1[1]
    orig_station = df_list1[2]
    orig_ground_station = df_list1[3]
    corr_ground_station = df_list1[4]
    attempt_merge = df_list1[6]

    # find the current dataframe
    df = df_list1[5]
    frame = df.frame.iloc[0]

    # create a hash on the combined columns (without the ground station)
    df['hash'] = pd.util.hash_pandas_object(df[[
        'orig_mh1_1','orig_mh2_1','orig_mhz_1',
        'orig_mh1_2','orig_mh2_2','orig_mhz_2',
        'orig_mh1_3','orig_mh2_3','orig_mhz_3',
        'orig_mh1_4','orig_mh2_4','orig_mhz_4',
        'frame']], 
        index=False)

    start_timestamp = df.corr_timestamp.iloc[0]

    # find the previous dataframe
    df_list_prev = sorted_df_list[index_no-1]
    # prev_starttime = df_list_prev[0]
    # prev_endtime = df_list_prev[1]
    
    prev_df = df_list_prev[5]
    # find the last record in the last timeseries
    end_frame = prev_df.frame.iloc[-1]
    end_timestamp = prev_df.corr_timestamp.iloc[-1]
    
    end_time_index = prev_df.time_index.iloc[-1]

    timestamp_time_diff = (start_timestamp - end_timestamp).total_seconds()
    
    starttime= UTCDateTime(start_timestamp)

    found_index = []
    method=''

    
    # logging.info('{} Original starttime (data time)={} 1st frame={} Original endtime={} Last frame={} ground_station={}'.format(index_no, starttime, frame, endtime, df.frame.iloc[-1],corr_ground_station))

    overlap_found = False
    if start_timestamp <= (end_timestamp + pd.Timedelta(seconds=2)):
        start_timestamp_match = start_timestamp - pd.Timedelta(seconds=2)
        end_timestamp_match = start_timestamp + pd.Timedelta(seconds=2)
        overlap_found = True

        # logging.info('B')
        # found_index = (prev_df[(prev_df.corr_timestamp >= start_timestamp) & 
        #   (prev_df.corr_timestamp <= end_timestamp) & 
        #   (prev_df.frame == frame)].index)
        # logging.info(frame)
        # logging.info(prev_df.to_string())
        found_index_all = prev_df[prev_df['frame'] == frame].index

        # logging.info(df.iloc[0:1].to_string())

        # logging.info(found_index_all)
        # for idx in found_index_all:
        #     logging.info(prev_df.iloc[idx:idx+1].to_string())

        
        found_index = (prev_df[(prev_df.corr_timestamp >= start_timestamp_match) & 
          (prev_df.corr_timestamp <= end_timestamp_match) & 
          (prev_df.frame == frame)].index)

        

        if len(found_index) == 1:
            found_i = found_index[0]
            # make some tests 
            prev_hash = prev_df.hash.iloc[found_i]
            current_hash = df.hash.iloc[0]
            prev_timestamp = prev_df.corr_timestamp.iloc[found_i]

            if prev_hash != current_hash:
                logging.info('WARNING: Not a perfect match:')
                # TODO  make a bit more of this 
                logging.info(prev_df.iloc[found_i:found_i+1].to_string())
                logging.info(df.iloc[0:1].to_string())
            # if abs(prev_timestamp - starttime) > 1:
            #     logging.info('WARNING: Time difference greater than 1 s.')        
            
            # find the matching time index (probably won't be the LAST time index)
            prev_time_index = prev_df.time_index[found_i]
            # prev_corr_timestamp = prev_df.corr_timestamp[found_i]
            # prev_ground_station = prev_df.corr_ground_station[found_i]
            # note that time index is every 4 records (because there are 
            # four records per frame number)
            df.time_index = np.arange(len(df)) + prev_time_index
            sample_time = sample_time0 + df.iloc[0].time_index *DELTA*4
            # sample_time = starttime0 + prev_time_index*DELTA*4
            # corr_timestamp = df.corr_timestamp.iloc[0]
            # ground_station = df.corr_ground_station[0]
            # t_diff = UTCDateTime(corr_timestamp)  - UTCDateTime(prev_corr_timestamp)
            # 
            # # logging.info('Previous trace {} ({}), New trace {} ({}), Time diff={}'.format(prev_corr_timestamp, prev_ground_station, corr_timestamp, ground_station, t_diff))
            # # logging.info(tdiff={})
            # method='overlap'
            # logging.debug('method={} ground station time diff={} time_index={}'.format(method, t_diff, df.time_index.iloc[0]))
            logging.info('INFO: Overlapping traces found. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={}'.format(index_no, orig_station, orig_ground_station, UTCDateTime(starttime), sample_time, UTCDateTime(starttime)-sample_time))

            # 
            # logging.info('prev')
            # logging.info(prev_df.iloc[found_i:found_i+1].to_string())
            # 
            # logging.info('current')
            # logging.info(df.iloc[0:1].to_string())

        elif len(found_index) > 1:
            print('EXCEPTION: Too many found.')
            logging.info('EXCEPTION: Too many found.')
            raise Exception


    # if the item has not been found yet, try a different way 
    if len(found_index) == 0:

        delta4_1 = prev_df.delta4.iloc[-1]
        delta4_2 = df.delta4.iloc[0]

        gradient = (delta4_1 + delta4_2)/2
        # XXXX

        frame_gap_correct, actual_frame_gap, absolute_error, percent_error= loose_frame_diff(last_frame=end_frame,new_frame=frame,gap=timestamp_time_diff,delta4_1=delta4_1,delta4_2=delta4_2)
        if frame_gap_correct:
            index_diff = actual_frame_gap
        else: 
            # making a guess
            index_diff = round(timestamp_time_diff / gradient)

        # find the time interval to get the sample time 
        intervals = end_time_index + index_diff
        df.time_index = np.arange(len(df)) + intervals
        sample_time = sample_time0 + df.iloc[0].time_index *DELTA*4

        if frame_gap_correct:
            logging.info('INFO: Time gap with correct number of missing samples found. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time ))
        else: 
            if overlap_found: 
                logging.info('SEVERE: Overlap found and time gap does not agree with frame gap. Overwriting part of trace. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time))
                attempt_merge = False
            else: 
                logging.info('SEVERE: Time gap does not agree with frame gap. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time))

    # TODO how do we sort later - should the timing be updated? 
    df_list1[6] = attempt_merge




def find_dir(top_level_dir,station,starttime,lower=True):
    year = str(starttime.year)
    day = str('{:03}'.format(starttime.julday))
    if lower:
        network='xa'
        station=station.lower()
    else:
        network='XA'
        station=station.upper()
    return os.path.join(top_level_dir,station,year,day)

def make_dir(top_level_dir,station,starttime,lower=True):
    # makes the directory if one is not found
    #<SDSDIR>/<YEAR>/<NET>/<STA>/<CHAN.TYPE>

# xa/ continuous_waveform/s12/1976/183/
# xa.s12.9.mh1.1976.183.1.mseed


    directory = find_dir(top_level_dir,station,starttime)
    if os.path.exists(directory):
        return directory
    # check that the overall directory exists
    elif not os.path.exists(top_level_dir):
        msg = ("The directory {} doesn't exist".format(top_level_dir))
        raise IOError(msg)
    else:
        year = str(starttime.year)
        day = str('{:03}'.format(starttime.julday))
        if lower:
            network='xa'
            station=station.lower()
        else:
            network='XA'
            station=station.upper()
        directory = os.path.join(top_level_dir, station)
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = os.path.join(directory,year)
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = os.path.join(directory, day)
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory

def filename_from_trace(trace,directory,lower=True):

    station = trace.stats.station
    location = trace.stats.location
    channel = trace.stats.channel
    year = str(trace.stats.starttime.year)
    julday = str('{:03}'.format(trace.stats.starttime.julday))
    starttime = trace.stats.starttime

    if lower:
        station=station.lower()
        channel=channel.lower()
        network='xa'
        ext='mseed'
    else:
        station=station.upper()
        channel=channel.upper()
        network='XA'
        ext='MSEED'

    # the file should have a revision of 0, which indicates that 
    # the timing has not been adjusted yet 
    rev = 0 
    
    # # check if the any exist in the folder
    # test_basename = '%s.%s.%s.%s.%s.%s.*.%s.%s' % (network,station, location, channel,
    #   year, julday, rev, ext)
    # test_filename = os.path.join(directory,test_basename)
    # file_no_lst = []
    # found_file_no = None
    # for filepath in glob.iglob(test_filename):
    #     st = read(filepath)
    #     # if st[0].stats.starttime != starttime:
    #     #     _, _, _, _, _, _, file_no, _, _  = filepath.split('.')
    #     #     file_no_lst.append(int(file_no))
    # 
    #     if st[0].stats.starttime == starttime:
    #         _, _, _, _, _, _, found_file_no, _, _  = filepath.split('.')
    #     else: 
    #         _, _, _, _, _, _, file_no, _, _  = filepath.split('.')
    #         file_no_lst.append(int(file_no))
    # 
    # if found_file_no is not None:
    #     file_no = found_file_no
    # elif len(file_no_lst) > 0:
    #     file_no = max(file_no_lst)
    #     file_no = file_no + 1
    # else:
    #     file_no = 1
    # 
    current_filename = '%s.%s.%s.%s.%s.%s.%s.%s' % (network,station, location, channel,
      year, julday, rev, ext)

    return current_filename

def save_stream(stream,join_dir):
    logging.debug('save_stream()')
    # for station in ('S12', 'S14', 'S15', 'S16'):
    #     for location in ('01','02','03','04','05','06','07','08','09',
    #           '10','11','12','13','14'):
    #         # location_stream = stream.select(location=location)
    #         location_stream = stream.select(station=station, location=location)

    for station in ['S11','S12','S14','S15','S16']:
        station_stream = stream.select(station=station)
        if len(station_stream) > 0:
            locations = set()
            for tr in station_stream:
                locations.add(tr.stats.location)
            for location in locations:
                location_stream = station_stream.select(location=location)
                # need directory with separate station, year and day
                # make the output directory
                station = location_stream[0].stats.station
                starttime = location_stream[0].stats.starttime
                directory = make_dir(join_dir,station,starttime,lower=True)
                # TODO do not overwrite if there is the same 
                # combination later in the day 

                for channel in ('AFR', 'ATT', 'MH1', 'MH2', 'MHZ', 'SHZ',):
                    # note that 'AFR only exists if 'config.combine_ground_stations == False
                    channel_stream = location_stream.select(channel=channel)
                    if len(channel_stream) > 0:
                        filename = filename_from_trace(channel_stream[0],directory)
                        filename = os.path.join(directory,filename)

                        # logging.info('Temporarily not writing {}'.format(filename))
                        # split the streams in order to save them
                        channel_stream = channel_stream.split()
                        
                        # if 'xa.s16..mh1.1976.064.1.0.mseed'  in filename:
                        #     # channel_stream.plot(method='full',size=(1200,600))
                        #     for tr in channel_stream:
                        #         logging.info(tr.stats)
                        #         channel_stream.write(filename, format='MSEED')
                        # save the streams
                        logging.info('Writing file {}'.format(filename))
                        channel_stream.write(filename, format='MSEED')
            print('Wrote files to ', directory)

    
        # if write_gzip:
        #     with open(out_filename, 'rb') as f_in, gzip.open(gzip_out_filename, 'wb') as f_out:
        #         shutil.copyfileobj(f_in, f_out)
        #     os.unlink(out_filename)

def loose_frame_diff(last_frame,new_frame,gap,delta4_1=None,delta4_2=None):
    
    if delta4_1 is not None and delta4_2 is not None:
        delta4 = (delta4_1 + delta4_2)/2
    else:
        if delta4_1 is not None:
            delta4 = delta4_1
        elif delta4_2 is not None:
            delta4 = delta4_1
        else:
            delta4 = DELTA*4

    logging.debug('loose_frame_diff {} {} {}'.format(last_frame,new_frame,gap))
    # tolerance (plus or minus) in seconds
    TOLERANCE = 1

    absolute_error = None
    percent_error = None
    actual_frame_gap = None
    frame_gap_correct = False

    # assuming that the frame diff is always positive
    # gap can be anything above -1 s, but is usually positive 
    if gap > -1:
        # logging.info('condition 3')

        # find difference between the last frame and the new frame 
        # (this number will ignore the full frames)
        diff = minus_frame(old_frame=last_frame,new_frame=new_frame)

        if gap > 0:
            # logging.info('condition 4')

            # get an approximate number for the frame gap 
            est_frame_gap = gap/delta4

            # get an approximate number for the full number of frames
            full_frames_est = int(est_frame_gap // 90)

            # full_frames_est can be wrong if we are very close to 
            # dividing by 90 frames, so we check the esimate minus and 
            # plus one, and then find the value of the single gap closest
            # to the nominal inteval value 
            # actual_frame_gap is the total number of frames in the gap (that 
            # fits with the last_frame and new_frame)
            actual_frame_gap_plus = (full_frames_est+1)*90 + diff

            single_gap_plus = gap/actual_frame_gap_plus
            single_gap = single_gap_plus
            actual_frame_gap = actual_frame_gap_plus

            actual_frame_gap_normal = full_frames_est*90 + diff

            if actual_frame_gap_normal != 0: 
                single_gap_normal = gap/actual_frame_gap_normal
                if abs(single_gap_normal - (delta4)) <  abs(single_gap - (delta4)):
                    single_gap = single_gap_normal
                    actual_frame_gap = actual_frame_gap_normal

            if full_frames_est-1 >= 0:
                actual_frame_gap_minus = (full_frames_est-1)*90 + diff
                if actual_frame_gap_minus != 0: 
                    single_gap_minus = gap/actual_frame_gap_minus
                    if abs(single_gap_minus - (delta4)) <  abs(single_gap - (delta4)):
                        single_gap = single_gap_minus
                        actual_frame_gap = actual_frame_gap_minus

        # if the gap was negative 
        else:
            full_frames_est = 0
            single_gap = gap/diff
            actual_frame_gap = diff
            # logging.info('condition 5')

        est_gap = actual_frame_gap * delta4
        # logging.info('actual_frame_gap {}, gap={}, est_gap={}, diff {},full_frames_est {},single_gap {}'.format(actual_frame_gap,gap,est_gap,diff,full_frames_est,gap/actual_frame_gap))

        # logging.info('gap diff {} single_gap {}'.format(gap-est_gap, single_gap))    
        # logging.info('gap={} est_gap={} actual_frame_gap={}'.format(gap,est_gap,actual_frame_gap))
        # frame_change = actual_frame_gap

        '''
        There are two sources of potential error:
        1. changing the ground station
        2. chnages to the delta which cannot be estimated from our simple method 

        The first is larger, and can be up to 1 s, without any sample gap. 

        The second is much smaller. 12 s in 24 hours = 0.014 %

        The two error types are not connected. 
        '''
        
        absolute_error = gap - est_gap

        if gap != 0:
            percent_error = abs(gap - est_gap)/est_gap
        else:
            percent_error = 0

        max_estimate = gap + (gap*MAX_PERCENT_ERROR/100) + TOLERANCE
        min_estimate = gap - (gap*MAX_PERCENT_ERROR/100) - TOLERANCE
        
        # check if it is within tolerance 
        if (est_gap < max_estimate) and (est_gap > min_estimate):
            # frame_change = actual_frame_gap
            # absolute_error = gap - est_gap
            # not a real exception 
            # if there is a real error, there's not that much we can do - 
            # this is still probably the best estimate
            # logging.info('SEVERE: Suitable gap not found est_gap={} percent_error={}, absolute_error={} '.format(est_gap, percent_error, absolute_error))
            # frame_change = -9999
            frame_gap_correct = True

        # logging.info('{} {}'.format(max_estimate, min_estimate))

    return frame_gap_correct, actual_frame_gap, absolute_error, percent_error

def initial_cleanup(df):

    df['orig_timestamp'] = df['orig_timestamp'].astype('datetime64[ns, UTC]')

    df['orig_station'] = df['orig_station'].astype('string')
    df['bit_synchronizer'] = df['bit_synchronizer'].astype('string')
    df['sync'] = df['sync'].astype('string')

    df['orig_no'] = to_Int64(df['orig_no'])
    df['clock_flag'] = to_Int64(df['clock_flag'])
    df['orig_frame'] = to_Int64(df['orig_frame'])
    df['orig_ground_station'] = to_Int64(df['orig_ground_station'])
    df['orig_mh1_1'] = to_Int64(df['orig_mh1_1'])
    df['orig_mh2_1'] = to_Int64(df['orig_mh2_1'])
    df['orig_mhz_1'] = to_Int64(df['orig_mhz_1'])
    df['orig_mh1_2'] = to_Int64(df['orig_mh1_2'])
    df['orig_mh2_2'] = to_Int64(df['orig_mh2_2'])
    df['orig_mhz_2'] = to_Int64(df['orig_mhz_2'])
    df['orig_mh1_3'] = to_Int64(df['orig_mh1_3'])
    df['orig_mh2_3'] = to_Int64(df['orig_mh2_3'])
    df['orig_mhz_3'] = to_Int64(df['orig_mhz_3'])
    df['orig_mh1_4'] = to_Int64(df['orig_mh1_4'])
    df['orig_mh2_4'] = to_Int64(df['orig_mh2_4'])
    df['orig_mhz_4'] = to_Int64(df['orig_mhz_4'])


    if 'shz_4' in df.columns:
        # also add in the empty columns
        df['shz_2'] = pd.NA

        if 'shz_46' not in df.columns:
            df['shz_46'] = pd.NA

        df['shz_56'] = pd.NA

        if df['orig_station'].iloc[0] == 'S15':
            # a problem with the data on SHZ_24 for S15
            df['shz_24'] = pd.NA

        df['shz_2'] = to_Int64(df['shz_2'])
        df['shz_4'] = to_Int64(df['shz_4'])
        df['shz_6'] = to_Int64(df['shz_6'])
        df['shz_8'] = to_Int64(df['shz_8'])
        # XXXX
        df['shz_10'] = to_Int64(df['shz_10'])
        df['shz_12'] = to_Int64(df['shz_12'])
        df['shz_14'] = to_Int64(df['shz_14'])
        df['shz_16'] = to_Int64(df['shz_16'])
        df['shz_18'] = to_Int64(df['shz_18'])
        df['shz_20'] = to_Int64(df['shz_20'])
        df['shz_22'] = to_Int64(df['shz_22'])
        df['shz_24'] = to_Int64(df['shz_24'])
        df['shz_26'] = to_Int64(df['shz_26'])
        df['shz_28'] = to_Int64(df['shz_28'])
        df['shz_30'] = to_Int64(df['shz_30'])
        df['shz_32'] = to_Int64(df['shz_32'])
        df['shz_34'] = to_Int64(df['shz_34'])
        df['shz_36'] = to_Int64(df['shz_36'])
        df['shz_38'] = to_Int64(df['shz_38'])
        df['shz_40'] = to_Int64(df['shz_40'])
        df['shz_42'] = to_Int64(df['shz_42'])
        df['shz_44'] = to_Int64(df['shz_44'])
        df['shz_46'] = to_Int64(df['shz_46'])
        df['shz_48'] = to_Int64(df['shz_48'])
        df['shz_50'] = to_Int64(df['shz_50'])
        df['shz_52'] = to_Int64(df['shz_52'])
        df['shz_54'] = to_Int64(df['shz_54'])
        df['shz_56'] = to_Int64(df['shz_56'])
        df['shz_58'] = to_Int64(df['shz_58'])
        df['shz_60'] = to_Int64(df['shz_60'])
        df['shz_62'] = to_Int64(df['shz_62'])
        df['shz_64'] = to_Int64(df['shz_64'])

        # since zero is also invalid, replace with pd.NA
        df.loc[df['shz_2'] == 0, 'shz_2'] = pd.NA
        df.loc[df['shz_4'] == 0, 'shz_4'] = pd.NA
        df.loc[df['shz_6'] == 0, 'shz_6'] = pd.NA
        df.loc[df['shz_8'] == 0, 'shz_8'] = pd.NA
        df.loc[df['shz_10'] == 0, 'shz_10'] = pd.NA

        df.loc[df['shz_12'] == 0, 'shz_12'] = pd.NA
        df.loc[df['shz_14'] == 0, 'shz_14'] = pd.NA
        df.loc[df['shz_16'] == 0, 'shz_16'] = pd.NA
        df.loc[df['shz_18'] == 0, 'shz_18'] = pd.NA
        df.loc[df['shz_20'] == 0, 'shz_20'] = pd.NA

        df.loc[df['shz_22'] == 0, 'shz_22'] = pd.NA
        df.loc[df['shz_24'] == 0, 'shz_24'] = pd.NA
        df.loc[df['shz_26'] == 0, 'shz_26'] = pd.NA
        df.loc[df['shz_28'] == 0, 'shz_28'] = pd.NA
        df.loc[df['shz_30'] == 0, 'shz_30'] = pd.NA

        df.loc[df['shz_32'] == 0, 'shz_32'] = pd.NA
        df.loc[df['shz_34'] == 0, 'shz_34'] = pd.NA
        df.loc[df['shz_36'] == 0, 'shz_36'] = pd.NA
        df.loc[df['shz_38'] == 0, 'shz_38'] = pd.NA
        df.loc[df['shz_40'] == 0, 'shz_40'] = pd.NA

        df.loc[df['shz_42'] == 0, 'shz_42'] = pd.NA
        df.loc[df['shz_44'] == 0, 'shz_44'] = pd.NA
        df.loc[df['shz_46'] == 0, 'shz_46'] = pd.NA
        df.loc[df['shz_48'] == 0, 'shz_48'] = pd.NA
        df.loc[df['shz_50'] == 0, 'shz_50'] = pd.NA

        df.loc[df['shz_52'] == 0, 'shz_52'] = pd.NA
        df.loc[df['shz_54'] == 0, 'shz_54'] = pd.NA
        df.loc[df['shz_56'] == 0, 'shz_56'] = pd.NA
        df.loc[df['shz_58'] == 0, 'shz_58'] = pd.NA
        df.loc[df['shz_60'] == 0, 'shz_60'] = pd.NA

        df.loc[df['shz_62'] == 0, 'shz_62'] = pd.NA
        df.loc[df['shz_64'] == 0, 'shz_64'] = pd.NA
    
    # logging.info(df['shz_4'].iloc[0])
    # 
    # logging.info(df.dtypes.to_string())
    # logging.info(df.head(90).to_string())
    


    df['frame'] = to_Int64(df['frame'])

    # replace empty ground station and orig station with values
    df['orig_ground_station'].fillna(-1,inplace=True)
    df['orig_station'].fillna('S0',inplace=True)
    df['orig_frame'].fillna(-99,inplace=True)
    df['frame'].fillna(-99,inplace=True)

    df['corr_timestamp'] = df['corr_timestamp'].astype('datetime64[ns, UTC]')
    df['corr_ground_station'] = to_Int64(df['corr_ground_station'])
    df['corr_gap_count'] = to_Int64(df['corr_gap_count'])
    df['time_index'] = 0
    df['delta4'] = df['delta4'].astype(float)

    return df

def despike3(trace):

    # method masks the trace where spikes are found

    rec_spikes_found = 0 

    TEST_SPIKE_V_SMALL = 5
    TEST_SPIKE_SMALL = 10
    TEST_SPIKE_HIGH = 30

    if len(trace.data) <= 3 or trace.stats.channel not in ('MHZ', 'MH1', 'MH2'):
        # return trace_data, rec_spikes_found
        return rec_spikes_found

    # mean = np.ma.mean(trace.data)

    # find the difference between two adjacent records by shifting the 
    # array and taking away from the original
    shift_plus1 = trace.data - shift(trace.data,1)
    # find the absolute value of the difference
    shift_plus1_abs = abs(shift_plus1)

    # find the mask for any values that were masked
    mask = np.ma.getmask(shift_plus1_abs)
    # find values greater than the TEST_SPIKE value
    idx_tuple = np.ma.where((mask == False) & (shift_plus1_abs > TEST_SPIKE_V_SMALL))
    idx_list = list(idx_tuple[0])
    # logging.info('idx_list {}'.format(idx_list))
        

    spikes_found_list = []

    # return 0

    if len(idx_list) > 0:
        for idx in idx_list:

            if idx > 0 and idx < len(shift_plus1_abs)-1:

                before_spike = trace.data[idx-1]
                current_spike = trace.data[idx]
                after_spike = trace.data[idx+1]
                
                before_shift =  shift_plus1[idx-1]
                current_shift =  shift_plus1[idx]
                after_shift =  shift_plus1[idx+1]

                diff = abs(before_spike - after_spike)

                spike_found = False
                # large spike, not very flat base 
                if spike_found == False:
                    if abs(current_shift) > TEST_SPIKE_HIGH and abs(after_shift) > TEST_SPIKE_HIGH:
                        # look for some symmetry - we need to remove the one that's in the middle, not
                        # the one that's fixing the spike 
                        if (diff < 10) and ((current_shift > 0 and after_shift < 0) or (current_shift < 0 and after_shift >  0)):
                            spike_found = True
                            # print('found large')
            
                # small spike, nearly flat base 
                if spike_found == False:
                    if abs(current_shift) > TEST_SPIKE_SMALL and abs(after_shift) > TEST_SPIKE_SMALL:
                        if (diff < 5) and ((current_shift > 0 and after_shift < 0) or (current_shift < 0 and after_shift >  0)):
                            spike_found = True
                            # print('found small')
            
                # very small spike, very flat base 
                if spike_found == False:
                    if abs(current_shift) > TEST_SPIKE_V_SMALL and abs(after_shift) > TEST_SPIKE_V_SMALL:
                        if (diff < 2) and ((current_shift > 0 and after_shift < 0) or (current_shift < 0 and after_shift >  0)):
                            spike_found = True
                            # print('found very small')


                if spike_found:
                    spikes_found_list.append(idx)


                # also check for NA 
                # a spike will have a null value on one side and a 
                # small value on the other 
                if np.ma.is_masked(before_shift):

                    if abs(after_shift) < TEST_SPIKE_V_SMALL:

                        spikes_found_list.append(idx-1)
                    
                if np.ma.is_masked(after_shift):

                    if abs(before_shift) < TEST_SPIKE_V_SMALL:
                        spikes_found_list.append(idx)

                # or they could both be null
                if  np.ma.is_masked(before_shift) and np.ma.is_masked(after_shift):
                    spikes_found_list.append(idx)

            elif idx == len(shift_plus1_abs)-1:
                # spike on last record of trace
                spikes_found_list.append(idx)



        if len(spikes_found_list) > 0:
            # remove any duplicates 
            spikes_found_list = list(set(spikes_found_list))
            spikes_found_list.sort()
            # print(' Spikes found ', spikes_found_list)
            if np.ma.isMaskedArray(trace.data) == False:
                trace.data = np.ma.masked_array(data=trace.data,
                             mask=False)
            trace.data[spikes_found_list] = np.ma.masked





                # diff < 2 and spike greater than TEST_SPIKE_V_SMALL
                # diff < 5 and spike greater than TEST_SPIKE_SMALL
                # diff < 10 and spike greater than TEST_SPIKE_HIGH







    rec_spikes_found = len(spikes_found_list)

    return rec_spikes_found

# def despike(channel_array):
# 
#     rec_spikes_found = 0 
# 
#     # find the difference between two adjacent records by shifting the 
#     # array and taking away from the original
#     shift_plus1 = channel_array - shift(channel_array,1)
#     # find the absolute value of the difference
#     shift_plus1_abs = abs(shift_plus1)
#     # find the mask for any values that were masked
#     mask = np.ma.getmask(shift_plus1_abs)
#     # find values greater than the TEST_SHIFT value
#     idx_tuple = np.ma.where((mask == False) & (shift_plus1_abs > TEST_SHIFT))
#     idx_list = list(idx_tuple[0])
#     logging.info('idx_list {}'.format(idx_list))
# 
#     if len(idx_list) > 0:
#         for i, idx in enumerate(idx_list):
#             if i > 0:
#                 # look for consecutive entries into the list
#                 prev_idx = idx_list[i-1]
#                 if idx - 1 == prev_idx:
#                     current =  shift_plus1[idx]
#                     previous =  shift_plus1[prev_idx]
#                     # if there is one positive value followed by a negative
#                     # or vice-versa, it's a spike
#                     if ((current > 0 and previous < 0) or
#                       (current < 0 and previous > 0)):
#                         if np.ma.is_masked(channel_array) == False:    
#                             channel_array = np.ma.masked_array(channel_array, mask=False)
#                         channel_array[prev_idx] = ma.masked
#                         rec_spikes_found += 1
# 
#     return channel_array, rec_spikes_found

# preallocate empty array and assign slice 
def shift(arr, num, fill_value=np.nan):
    result = np.ma.empty_like(arr)
    if num > 0:
        result[:num] = ma.masked
        result[num:] = arr[:-num]
    elif num < 0:
        result[num:] = ma.masked
        result[:num] = arr[-num:]
    else:
        result[:] = arr
    return result

    

def stream_import(df,sample_time0,index0,attempt_merge):
    """
    Makes a stream file from the gzipped csv file.

    csv_import() makes streams from the individual csv file.

    :type gzip_filename: str
    :param gzip_filename: gzipped filename of the CSV file to be read.
    :rtype: :class:`~obspy.core.stream.Stream`
    :return: A ObsPy Stream object.

    """

    rec_spikes_found = 0

    stream = Stream()

    mh1 = df[['orig_mh1_1', 'orig_mh1_2','orig_mh1_3','orig_mh1_4']].to_numpy(dtype='int32',na_value=INVALID)
    mh1 = mh1.flatten()
    # if INVALID in mh1:
    #     mh1 = ma.masked_equal(mh1, INVALID)
    # if config.clean_spikes: 
    #     mh1, rec_spikes_found1 = despike(mh1)
    #     rec_spikes_found += rec_spikes_found1

    mh2 = df[['orig_mh2_1', 'orig_mh2_2','orig_mh2_3','orig_mh2_4']].to_numpy(dtype='int32',na_value=INVALID)
    mh2 = mh2.flatten()
    # if INVALID in mh2:
    #     mh2 = ma.masked_equal(mh2, INVALID)
    # if config.clean_spikes: 
    #     mh2, rec_spikes_found1 = despike(mh2)
    #     rec_spikes_found += rec_spikes_found1

    mhz = df[['orig_mhz_1', 'orig_mhz_2','orig_mhz_3','orig_mhz_4']].to_numpy(dtype='int32',na_value=INVALID)
    mhz = mhz.flatten()
    # if INVALID in mhz:
    #     mhz = ma.masked_equal(mhz, INVALID)
    # if config.clean_spikes: 
    #     mhz, rec_spikes_found1 = despike(mhz)
    #     rec_spikes_found += rec_spikes_found1

    if 'shz_4' in df.columns:
    # logging.info('head')
    # logging.info(df.head().to_string())
        shz = df[[
            'shz_2', 'shz_4','shz_6','shz_8', 'shz_10',
            'shz_12', 'shz_14','shz_16','shz_18', 'shz_20',
            'shz_22', 'shz_24','shz_26','shz_28', 'shz_30',
            'shz_32', 'shz_34','shz_36','shz_38', 'shz_40',
            'shz_42', 'shz_44','shz_46','shz_48', 'shz_50',
            'shz_52', 'shz_54','shz_56','shz_58', 'shz_60',
            'shz_62', 'shz_64'
            ]].to_numpy(dtype='int32',na_value=INVALID)
        shz = shz.flatten()

# data_type should be 'float64' for ATT and AFR and
# 'int32' for MH1, MH2 and MHZ

    frame = df[['frame']].to_numpy(dtype='int32',na_value=INVALID)
    # frame = np.repeat(frame,4)
    # if INVALID in frame:
    #     frame = ma.masked_equal(frame, INVALID)
    frame = frame.flatten()

    # get the absolute times in seconds since 1 Jan 1970
    df['abs_times'] = df['corr_timestamp'].values.astype(np.int64) / 10 ** 9

    abs_times = df[['abs_times']].to_numpy(dtype='float64',na_value=None)
    # abs_times = np.repeat(abs_times,4)
    abs_times = abs_times.flatten()

    ground_station = str(df['corr_ground_station'].iloc[0])
    
    station = df['orig_station'].iloc[0]

    # set the starttime based on the estimate using the DELTA, not
    # using the the first time in the timestamp 

    sample_time = sample_time0 + df.time_index.iloc[0] * DELTA * 4
    # adjustment_time = UTCDateTime(abs_times[0]) - starttime_adjust
    # logging.info('clock starttime={} Adjusting time {}s'.format(UTCDateTime(abs_times[0]),adjustment_time))

    # TODO - check that we don't mess up the time differences here!!!!

    tr_MH1 = Trace(data=mh1)
    tr_MH1.stats.delta = DELTA
    tr_MH1.stats.network = NETWORK
    tr_MH1.stats.station = station
    tr_MH1.stats.ground_station = ground_station
    tr_MH1.stats.channel = 'MH1'
    tr_MH1.stats.starttime = sample_time
    tr_MH1.stats.attempt_merge = attempt_merge
    tr_MH1.stats.location = '00'

    tr_MH2 = Trace(data=mh2)
    tr_MH2.stats.delta = DELTA
    tr_MH2.stats.network = NETWORK
    tr_MH2.stats.station = station
    tr_MH2.stats.ground_station = ground_station
    tr_MH2.stats.channel = 'MH2'
    tr_MH2.stats.starttime = sample_time
    tr_MH2.stats.attempt_merge = attempt_merge 
    tr_MH2.stats.location = '00' 

    tr_MHZ = Trace(data=mhz)
    tr_MHZ.stats.delta = DELTA
    tr_MHZ.stats.network = NETWORK
    tr_MHZ.stats.station = station
    tr_MHZ.stats.ground_station = ground_station
    tr_MHZ.stats.channel = 'MHZ'
    tr_MHZ.stats.starttime = sample_time
    tr_MHZ.stats.attempt_merge = attempt_merge
    tr_MHZ.stats.location = '00'

    tr_ATT = Trace(data=abs_times)
    tr_ATT.stats.delta = DELTA*4
    tr_ATT.stats.network = NETWORK
    tr_ATT.stats.station = station
    tr_ATT.stats.ground_station = ground_station
    tr_ATT.stats.channel = 'ATT'
    tr_ATT.stats.starttime = sample_time
    tr_ATT.stats.location = ''

    # append the trace
    stream.append(tr_MH1)
    stream.append(tr_MH2)
    stream.append(tr_MHZ)
    stream.append(tr_ATT)

    if config.combine_ground_stations == False: 
        tr_AFR = Trace(data=frame)
        tr_AFR.stats.delta = DELTA*4
        tr_AFR.stats.network = NETWORK
        tr_AFR.stats.station = station
        tr_AFR.stats.ground_station = ground_station
        tr_AFR.stats.channel = 'AFR'
        tr_AFR.stats.starttime = sample_time
        tr_AFR.stats.location = ''
        stream.append(tr_AFR)


    # only create an SHZ trace if there are any data
    if 'shz_4' in df.columns:
        # if all the values are invalid, don't create an SHZ trace
        if np.where(shz != INVALID, True, False).any():
            tr_SHZ = Trace(data=shz)
            tr_SHZ.stats.delta = DELTA/8
            tr_SHZ.stats.network = NETWORK
            tr_SHZ.stats.station = station
            tr_SHZ.stats.ground_station = ground_station
            tr_SHZ.stats.channel = 'SHZ'
            tr_SHZ.stats.starttime = sample_time
            tr_SHZ.stats.attempt_merge = attempt_merge  
            tr_SHZ.stats.location = ''
            stream.append(tr_SHZ)

    return stream

def to_Int64(df_column):
    # Several columns in the dataframe requre integers that deal with 
    # null values. Int64 does this, but reading the dataframe into Int64 is
    # really slow.
    # Instead, we can transform it like this: 
    df_column.fillna(value=-999999,inplace=True)
    df_column = df_column.astype('int64')
    df_column = df_column.astype('Int64')
    df_column.replace(-999999, pd.NA,inplace=True)
    return df_column

def to_Int64_with_invalid(df_column):
    # Instead, we can transform it like this: 
    df_column.fillna(value=INVALID,inplace=True)
    df_column = df_column.astype('int64')
    df_column = df_column.astype('Int64')
    df_column.replace(INVALID, pd.NA,inplace=True)
    return df_column