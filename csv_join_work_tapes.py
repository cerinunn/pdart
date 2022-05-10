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
from pdart.util import relative_timing_trace
# from pdart.csv_import_work_tapes import find_output_dir, make_output_dir, make_filelist
# import matplotlib.pyplot as plt

# Qt5Agg seems to work best on Mac - try 'TkAgg' if that works for you
# put this after the other imports, otherwise it can be overridden
import matplotlib  
# matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

from pdart.util import maximize_plot

global total
total = 0
global gst_total
gst_total = 0
global df_manual_jump_correction
df_manual_jump_correction = None
global df_manual_clock_correction 
df_manual_clock_correction = None
global df_manual_exclude
df_manual_exclude = None
global manual_only
manual_only = False

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
# MAX_PERCENT_ERROR= 0.014
# try divergence that's a bit more realistic 0.5 seconds in 24 hours 
MAX_PERCENT_ERROR= 0.0005

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


def call_csv_join_work_tapes_manual(
    processed_dir='.',
    join_dir='.',
    log_dir='.',
    wildcard_style='*',
    year_start=None,
    year_end=None,
    day_start=None,
    day_end=None,
    stations=['S11','S12','S14','S15','S16'],
    manual_clock_correction=None,
    manual_jump_correction=None,
    manual_exclude=None,
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
    global manual_only
    manual_only=True

    # note that these call the files again if the same station and day  
    # combination is called more than once
    # if this is used a lot, it should be rewritten

    global df_manual_clock_correction
    if manual_clock_correction is not None:
        df_manual_clock_correction = pd.read_csv(manual_clock_correction, dtype=str, comment='#')
        for f in df_manual_clock_correction.filename:
            # wtn.16.6.S14.3.1977.140.11_30_27_769000.csv.gz
            # pse.a16.4.93.S16.1.1974.001.00_00_00_300000.csv.gz
            lst = f.split('.')
            julday = lst[-4]
            year = lst[-5]
            station = lst[-7]
            call_csv_join_work_tapes(
                processed_dir=processed_dir,
                join_dir=join_dir,
                log_dir=log_dir,
                year_start=int(year),
                year_end=int(year),
                day_start=int(julday),
                day_end=int(julday),
                stations=[station],
                manual_clock_correction=manual_clock_correction,
                manual_jump_correction=manual_jump_correction
            )

    global df_manual_jump_correction
    if manual_jump_correction is not None:
        df_manual_jump_correction = pd.read_csv(manual_jump_correction, dtype=str, comment='#')
        # if len():
        df_manual_jump_correction['correction'] = df_manual_jump_correction['correction'].astype(float)
        for f in df_manual_jump_correction.filename:
            # wtn.16.6.S14.3.1977.140.11_30_27_769000.csv.gz
            # pse.a16.4.93.S16.1.1974.001.00_00_00_300000.csv.gz
            lst = f.split('.')
            julday = lst[-4]
            year = lst[-5]
            station = lst[-7]
            call_csv_join_work_tapes(
                processed_dir=processed_dir,
                join_dir=join_dir,
                log_dir=log_dir,
                year_start=int(year),
                year_end=int(year),
                day_start=int(julday),
                day_end=int(julday),
                stations=[station],
                manual_clock_correction=manual_clock_correction,
                manual_jump_correction=manual_jump_correction
            )

    global df_manual_exclude
    if manual_exclude is not None:
        df_manual_exclude = pd.read_csv(manual_exclude, dtype=str, comment='#')
        # if len():
        for f in df_manual_exclude.filename:
            # wtn.16.6.S14.3.1977.140.11_30_27_769000.csv.gz
            # pse.a16.4.93.S16.1.1974.001.00_00_00_300000.csv.gz
            lst = f.split('.')
            julday = lst[-4]
            year = lst[-5]
            station = lst[-7]
            call_csv_join_work_tapes(
                processed_dir=processed_dir,
                join_dir=join_dir,
                log_dir=log_dir,
                year_start=int(year),
                year_end=int(year),
                day_start=int(julday),
                day_end=int(julday),
                stations=[station],
                manual_clock_correction=manual_clock_correction,
                manual_jump_correction=manual_jump_correction,
                manual_exclude=manual_exclude
            )



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
    manual_clock_correction=None,
    manual_jump_correction=None,
    manual_exclude=None,
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

    global df_manual_clock_correction
    if manual_clock_correction is not None:
        df_manual_clock_correction = pd.read_csv(manual_clock_correction, dtype=str, comment='#')

    global df_manual_jump_correction
    if manual_jump_correction is not None:
        df_manual_jump_correction = pd.read_csv(manual_jump_correction, dtype=str, comment='#')
        # if len():
        df_manual_jump_correction['correction'] = df_manual_jump_correction['correction'].astype(float)

    global df_manual_exclude
    if manual_exclude is not None:
        df_manual_exclude = pd.read_csv(manual_exclude, dtype=str, comment='#')


    # TODO -fix these date ranges - they are not quite right 
    # if no filenames have been passed, then look for them 
    for year in range(year_start,year_end+1):
        for day in range(day_start,day_end+1):

            if config.combine_ground_stations == True:
                if manual_only == False: 
                    log_filename = 'join.{}.{}.log'.format(year,day)
                else:
                    sta = stations[0]
                    log_filename = 'join_manual.{}.{}.{}.log'.format(year,day,sta)
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

            logging.info('config.combine_ground_stations={}'.format(config.combine_ground_stations))  
            logging.info('config.clean_spikes={}'.format(config.clean_spikes))  
            logging.info('config.view_corrected_traces={}'.format(config.view_corrected_traces))  
            logging.info('config.fix_clock_error={}'.format(config.fix_clock_error))  
            logging.info('config.fix_jump_error={}'.format(config.fix_jump_error))
            logging.info('config.exclude_masked_sections={}'.format(config.exclude_masked_sections))
                                
            for station in stations:
                df_list = []
                # wildcard_filename = 'pse.a16.1.83.S16.12.1972.195.00_00_00_159000*XXX.csv.gz'
                # wildcard_filename = 'pse.a15.2.168.S15.12.1972.195.00_00_00_144000.csv.gz'
                wildcard_filename = '{}.*.*.{}.*.{}.{:03}.*.csv.gz'.format(wildcard_style,station,year,day)
                print('wildcard filename ', processed_dir, wildcard_filename)
                # read in each of the csv files for the station
                # for filename in glob.glob('/Users/cnunn/python_packages/pdart/examples/test.csv'):
                for filename in glob.glob(os.path.join(processed_dir,wildcard_filename)):

                    # 
                    # temporary
                    # if os.path.basename(filename) not in ('wtn.1.3.S12.1.1976.061.13_29_58_945000.csv.gz'):
                    #     continue


                    exclude = manual_exclusion(filename)
                    if exclude:
                        logging.info('WARNING: manually excluding file {}'.format(filename))
                        continue



# wtn.15.14.S15.7.1977.120.00_00_00_428000.csv.gz
# wtn.15.14.S15.1.1977.120.01_54_00_296000.csv.gz
# wtn.15.14.S15.11.1977.120.03_00_17_786000.csv.gz
# wtn.15.14.S15.8.1977.120.03_25_41_044000.csv.gz
# wtn.15.15.S15.8.1977.120.06_22_34_958000.csv.gz
# wtn.15.15.S15.4.1977.120.07_50_50_438000.csv.gz
# wtn.15.15.S15.5.1977.120.10_18_09_924000.csv.gz
# wtn.15.15.S15.1.1977.120.17_08_59_828000.csv.gz
# wtn.15.16.S15.1.1977.120.19_00_48_689000.csv.gz
# wtn.15.16.S15.7.1977.120.19_49_59_818000.csv.gz
                    #temp ones
                    # if os.path.basename(filename) in ('wtn.15.15.S15.5.1977.120.10_18_09_924000.csv.gz'):
                    #     continue


        
                    # if os.path.basename(filename) not in ('wtn.15.15.S15.8.1977.120.06_22_34_958000.csv.gz', 
                    #     'wtn.15.15.S15.1.1977.120.17_08_59_828000.csv.gz'):
                    #     continue

                    # # temp ones: 
                    # if os.path.basename(filename) in ('wtn.18.37.S12.3.1977.203.17_35_29_067000.csv.gz', 
                    #    'wtn.18.38.S12.3.1977.203.17_49_36_777000.csv.gz'):
                    #     continue


                    # print(filename)
                    
                    # if filename not in ('/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.6.S12.6.1976.063.00_00_00_199000.csv.gz',
                    #     '/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.6.S12.601.1976.063.00_32_59_285000.csv.gz',
                    #     '/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES/wtn.1.6.S12.9.1976.063.05_35_56_339000.csv.gz'):

                    # 400/800 error
                    # if filename not in ('/Users/cnunn/lunar_data/PDART_PROCESSED/wtn.9.17.S12.4.1976.299.00_00_00_034000.csv.gz'):
                    #     continue

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

                    # if os.path.basename(filename) in ['pse.a16.2.156.S16.3.1973.072.14_48_48_238000.csv.gz']:
                    #     continue

                    # if os.path.basename(filename) not in ('wtn.9.17.S12.4.1976.299.00_00_00_034000.csv.gz'):
                    #     continue

                    # if os.path.basename(filename) not in ['pse.a16.2.156.S16.3.1973.072.16_35_39_642000.csv.gz','pse.a16.2.156.S16.9.1973.072.17_09_12_466000.csv.gz' ]:
                    #     continue

# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a12.6.208.S12.9.1973.073.18_48_36_308000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a12.6.208.S12.9.1973.073.19_22_27_103000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a12.6.208.S12.9.1973.073.20_08_57_714000.csv.gz


#                     if os.path.basename(filename) in ['pse.a12.6.208.S12.9.1973.073.20_08_57_714000.csv.gz','pse.a12.6.208.S12.9.1973.073.19_22_27_103000.csv.gz','pse.a12.6.208.S12.9.1973.073.18_48_36_308000.csv.gz','pse.a12.6.208.S12.9.1973.073.18_35_03_360000.csv.gz','wtn.11.20.S12.700.1976.001.00_04_27_250000.csv.gz',
# # 'pse.a12.10.78.S12.4.1976.001.00_52_41_266000.csv.gz'
#                         ]:
#                         continue

                    # if os.path.basename(filename) in ['pse.a16.5.37.S16.3.1974.123.15_58_30_103000.csv.gz']:
                    #     continue





# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.2.156.S16.9.1973.072.17_09_12_466000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.2.157.S16.9.1973.072.20_20_00_361000.csv.gz

                    if 'dropped' not in filename:
                        try: 
                            gzip_filename = filename

                            # if gzip_filename not in ('/Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.3.1973.184.07_34_51_966000.csv.gz', '/Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.3.1973.184.09_00_45_934000.csv.gz','/Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.9.1973.184.11_52_07_470000.csv.gz'):
                            #     continue

# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.3.1973.184.07_34_51_966000.csv.gz

# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.4.1973.184.00_00_00_011000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.4.1973.184.01_40_49_591000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.4.1973.184.01_54_52_436000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.4.1973.184.02_34_22_154000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.0.1973.184.06_49_03_100000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.3.1973.184.07_34_51_966000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.11.1973.184.15_46_40_917000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.91.S16.11.1973.184.20_20_00_595000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.4.1973.184.00_06_40_902000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.4.1973.184.00_01_05_820000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.91.S16.4.1973.184.23_34_42_581000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.9.1973.184.11_52_07_470000.csv.gz
# /Users/cnunn/lunar_data/PDART_PROCESSED/pse.a16.3.90.S16.3.1973.184.09_00_45_934000.csv.gz


                            logging.info(gzip_filename)
                            df_list = read_file(filename,df_list,logging_level,join_dir,log_dir)
                        except Exception as e:
                            logging.info(traceback.format_exc())
                            logging.info(traceback.print_exc(file=sys.stdout))
                            # reset_config()
                            print('Warning, continuing')
                            # print('Not continuing')

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

def read_file(gzip_filename,df_list,end_logging_level=logging.INFO,join_dir='.',log_dir='.',
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

    # leap second correction 
    df = leap_correction(df,gzip_filename)

    df = manual_correction(df,gzip_filename)

    # df_end_frames = calculate_end_frames(df,df_end_frames=df_end_frames)

    # calculate the delta4 for the segment 
    segment_delta4, segment_est_divergence = calculate_segment_delta4(df)
    
#     if os.path.basename(gzip_filename) in ['pse.a15.4.65.S15.11.1973.074.01_05_34_360000.csv.gz']:
# # 'pse.a16.5.37.S16.3.1974.123.15_58_30_103000.csv.gz']:
#         print('yes')
#         df.clock_flag = 1


    # get some details about the data 
    starttime = UTCDateTime(df.corr_timestamp.iloc[0])
    endtime = UTCDateTime(df.corr_timestamp.iloc[-1])
    start_clock_flag = df.clock_flag.iloc[0]
    end_clock_flag = df.clock_flag.iloc[-1]
    orig_station = df.orig_station.iloc[0]
    orig_ground_station = df.orig_ground_station.iloc[0]
    corr_ground_station = df.corr_ground_station.iloc[0]

    attempt_merge = True

    df_dict = {'starttime' : starttime, 'endtime' : endtime, 
      'orig_station' : orig_station , 'orig_ground_station' : orig_ground_station, 
      'corr_ground_station' : corr_ground_station, 'df' : df, 
      'attempt_merge' : attempt_merge, 'segment_delta4' : segment_delta4,
      'segment_est_divergence' : segment_est_divergence, 'gzip_filename' : gzip_filename}
    # make a list of the dataframes
    df_list.append(df_dict)

    return df_list

def manual_exclusion(gzip_filename):
    exclude = False
    if df_manual_exclude is not None: 
        if len(df_manual_exclude) > 0:
            idx_list = df_manual_exclude[df_manual_exclude['filename'] == os.path.basename(gzip_filename)].index.tolist()
            if len(idx_list) > 0: 
                exclude = True

    return exclude

def manual_correction(df,gzip_filename):
    if df_manual_jump_correction is not None: 
        
        if len(df_manual_jump_correction) > 0: 
            idx_list = df_manual_jump_correction[df_manual_jump_correction['filename'] == os.path.basename(gzip_filename)].index.tolist()
            if len(idx_list) > 0:
                
                idx = idx_list[0]

                correction = df_manual_jump_correction.correction.iloc[idx]

                start_timetamp = df['orig_timestamp'].iloc[0]
                end_timestamp = df['orig_timestamp'].iloc[-1]

                logging.info('WARNING: corrected tapehead jump from original start timestamp {} to original end timestamp {} correction {:.01f} s [manual correction]'.format(start_timetamp,end_timestamp,correction))

                df['corr_timestamp'] = df['corr_timestamp'] - pd.Timedelta(correction, unit='seconds')
                df['orig_timestamp'] = df['orig_timestamp'] - pd.Timedelta(correction, unit='seconds')

                # may need to drop the first couple of records 
                to_drop = []

                julday = UTCDateTime(start_timetamp).julday
                if UTCDateTime(df.corr_timestamp.iloc[0]).julday != julday:
                    to_drop.append(0)
                if UTCDateTime(df.corr_timestamp.iloc[1]).julday != julday:
                    to_drop.append(1)
                if len(to_drop) > 0: 
                    df.drop(df.index[to_drop],inplace=True)
                    df.reset_index(inplace=True,drop=True)
                        
            # logging.info()


        # TODO add the correction thats been made manually to the log 
    if df_manual_clock_correction is not None: 
        if len(df_manual_clock_correction) > 0:
            idx_list = df_manual_clock_correction[df_manual_clock_correction['filename'] == os.path.basename(gzip_filename)].index.tolist()
            if len(idx_list) > 0: 
                df.clock_flag = 1

    return df

def leap_correction(df, filename):
    # make a correction for the leap seconds
    # The UTC added one second at the end of the year from 1972 to 1979.  
    # Please look at the wikipedia web site: https://en.wikipedia.org/wiki/Leap_second.  
    if os.path.basename(filename) in [
      'pse.a12.6.136.S12.5.1973.001.00_00_09_330000.csv.gz',
        'pse.a14.4.166.S14.5.1973.001.00_00_32_799000.csv.gz',
        'pse.a15.3.162.S15.3.1973.001.03_59_06_810000.csv.gz',
        'pse.a15.3.162.S15.3.1973.001.00_00_09_103000.csv.gz',
        'pse.a16.2.85.S16.3.1973.001.00_00_10_037000.csv.gz',

        'pse.a12.7.211.S12.1.1974.001.00_00_06_476000.csv.gz',
        'pse.a14.7.1.S14.1.1974.001.00_00_11_038000.csv.gz',
        'pse.a15.6.1.S15.1.1974.001.00_00_10_969000.csv.gz',
        'pse.a16.4.94.S16.1.1974.001.00_00_11_168000.csv.gz',


        'pse.a12.9.1.S12.2.1975.001.00_00_10_703000.csv.gz',
        'pse.a12.9.1.S12.2.1975.001.03_53_57_272000.csv.gz',
        'pse.a14.9.7.S14.2.1975.001.00_00_23_815000.csv.gz',
        'pse.a15.8.7.S15.2.1975.001.00_00_10_932000.csv.gz',
        'pse.a16.6.100.S16.2.1975.001.00_00_10_238000.csv.gz',

        'pse.a12.10.78.S12.5.1976.001.00_00_09_486000.csv.gz',
        'pse.a14.11.20.S14.5.1976.001.00_00_12_571000.csv.gz',
        'pse.a15.10.23.S15.5.1976.001.00_00_07_307000.csv.gz',
        'pse.a16.8.116.S16.5.1976.001.00_00_12_276000.csv.gz',

        'wtn.11.23.S12.1.1977.001.00_12_50_785000.csv.gz',
        'wtn.11.23.S14.1.1977.001.00_12_50_927000.csv.gz',
        'wtn.11.23.S15.1.1977.001.00_12_50_923000.csv.gz',
        'wtn.11.23.S16.1.1977.001.00_12_50_863000.csv.gz'

]:
        df['corr_timestamp'] = df['corr_timestamp'] + pd.Timedelta(-1, unit='seconds')
        df['orig_timestamp'] = df['orig_timestamp'] + pd.Timedelta(-1, unit='seconds')

    return df 


def process_list(df_list,join_dir):

    # for each station, find the first anchor
    # run through the rest of the dataframes for the station


    orig_station = df_list[0]['orig_station']
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
    # sorted_df_list = sorted(df_list, key=operator.itemgetter(0))
    sorted_df_list =  (sorted(df_list, key = lambda i: i['starttime']))

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
                    starttime = df_list1['starttime']

                    logging.info('WARNING: Deleting two records from trace {}, to continue from trace {}'.format(i, continue_idx))

                    # these starttimes are close, so delete the first couple of 
                    # records in this trace 
                    df = df_list1['df']
                    df = df.iloc[2:]
                    # reindex
                    df.reset_index(inplace=True,drop=True)
                    new_starttime = UTCDateTime(df.corr_timestamp.iloc[0])
                    sorted_df_list[i]['starttime'] = new_starttime
                    sorted_df_list[i]['df'] = df
                    
            # re-sort by starttime
            # sorted_df_list = sorted(df_list, key=operator.itemgetter(0))
            sorted_df_list =  (sorted(df_list, key = lambda i: i['starttime']))
            # now it's been re-sorted, the continuing index will be first 
            logging.info('WARNING: Traces resorted')
            continue_idx = 0

    

    # sometimes, there's more than one record that could be the start record 
    for i, df_list1 in enumerate(sorted_df_list):
        if i == 0:
            starttime_i0 = df_list1['starttime']
        if i > 0:
            starttime = df_list1['starttime']
            if starttime < (starttime_i0+ 1):
                # these starttimes are close, so delete the first couple of 
                # records in this trace 
                df = df_list1['df']
                df = df.iloc[2:]
                logging.info('WARNING: Deleting two records from trace {}, to continue from initial trace {}'.format(i, 0))
                # reindex
                df.reset_index(inplace=True,drop=True)
                sorted_df_list[i]['starttime'] = df.corr_timestamp.iloc[0]
                sorted_df_list[i]['df'] = df
            else:
                break


    # if the first record contains a clock flag, it's not very helpful 
    for i, df_list1 in enumerate(sorted_df_list):
        if pd.isna(df_list1['segment_est_divergence']):
            logging.info('SEVERE: Removing first item from list because it contains clock flags {}'.format(os.path.basename(df_list1['gzip_filename'])))
            sorted_df_list.pop(i)
        else:
            break  
    

    logging.info('#########')

    # calculate if any records might be suspect (use the Interquartile range)
    segment_est_divergence = []
    for x in sorted_df_list:
        logging.info(os.path.basename(x['gzip_filename']))
        if pd.notnull(x['segment_est_divergence']):
            segment_est_divergence.append(x['segment_est_divergence'])

    logging.info('#########')

    segment_est_divergence = np.array(segment_est_divergence)

    sort_data = np.sort(segment_est_divergence)

    Q1 = np.percentile(sort_data, 25, interpolation = 'midpoint') 
    Q2 = np.percentile(sort_data, 50, interpolation = 'midpoint') 
    Q3 = np.percentile(sort_data, 75, interpolation = 'midpoint') 
      
    # print('Q1 25 percentile of the given data is, ', Q1)
    # print('Q1 50 percentile of the given data is, ', Q2)
    # print('Q1 75 percentile of the given data is, ', Q3)
      
    IQR = Q3 - Q1 
    # print('Interquartile range is', IQR)

    low_lim = Q1 - 1.5 * IQR
    up_lim = Q3 + 1.5 * IQR


    config.potentially_suspect = 0
    for i, df_list1 in enumerate(sorted_df_list):
        # print(df_list1['segment_est_divergence'])
        if pd.notnull(df_list1['segment_est_divergence']):
            if (df_list1['segment_est_divergence'] < low_lim) or (df_list1['segment_est_divergence'] > up_lim):
                df_list1['reject'] = 'MAYBE'
                config.potentially_suspect += 1
            else:
                df_list1['reject'] = 'NO'
            if i == 0 and df_list1['reject'] == 'MAYBE':
                orig_station = df_list1['orig_station']
                orig_ground_station = df_list1['orig_ground_station']
                logging.info('SEVERE: Ground station/Station - {} {} First record is potentially suspect. Estimated divergence: {:.02f}'.format(orig_ground_station,orig_station, df_list1['segment_est_divergence']))
        else:
            df_list1['reject'] = 'NO'
    # check last record
    if df_list1['reject'] == 'MAYBE':
        orig_station = df_list1['orig_station']
        orig_ground_station = df_list1['orig_ground_station']
        logging.info('SEVERE: Ground station/Station - {} {} Last record is potentially suspect. Estimated divergence: {:.02f}'.format(orig_ground_station,orig_station, df_list1['segment_est_divergence']))
        

    orig_station = df_list1['orig_station']
    logging.info('INFO: Station - {} Divergence Lower limit={:.01f} Upper limit={:.01f} '.format(orig_station,low_lim,up_lim))
    if config.potentially_suspect > 0: 
        logging.info('WARNING: Station - {} {} record(s) are potentially suspect.'.format(orig_station,config.potentially_suspect))

    config.rejected = 0
    starttime0 = None
    index0 = None
    for i, df_list1 in enumerate(sorted_df_list):
        # df = df_list1['df']

        # logging.info(i)
        # logging.info(df.head().to_string())
        # logging.info(df.tail().to_string())

        # for each station, set up the time_index
        if i==0:
            starttime0, sample_time0 = first_record(index_no=i,sorted_df_list=sorted_df_list)
        elif i > 0:
            later_records(index_no=i,sample_time0=sample_time0,sorted_df_list=sorted_df_list)
                
        df_list1 = sorted_df_list[i]
        df = df_list1['df']
        if len(df) > 0:
            attempt_merge = df_list1['attempt_merge']
            # logging.info('attempt_merge {} {}'.format(i, attempt_merge))

            gaps_7777 = (df['corr_gap_count'] == -7777).sum()
            if gaps_7777 > 0:
        # df_list.append([starttime, endtime, orig_station,orig_ground_station, corr_ground_station, df])
                orig_station = df_list1['orig_station']
                orig_ground_station = df_list1['orig_ground_station']
                logging.info('WARNING: Ground station/Station - {} {} {} frames are reset to zero (-7777 error)'.format(orig_ground_station,orig_station,gaps_7777))

            # # 
            # logging.info('df?')
            # logging.info(df.head().to_string())
            # logging.info(df.dtypes.to_string())

            # import the stream for every record, using the same 
            # time_index0 and index0
            stream1 = stream_import(df,sample_time0,index0,attempt_merge)

            stream += stream1

    # stream.select(channel='CLK').plot()
    # print(stream.select(channel='CLK'))


    # if required, clean the spikes 
    if config.clean_spikes:
        for tr in stream:
            rec_spikes_found = despike3(tr)

    # if required, merge the streams 
    if config.combine_ground_stations: 
        # also includes clock correction
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
                        logging.info('INFO: Unable to merge because {} ({:.1f}%) of the records are different. {} {}'.format(count_overlap, percent_diff, tr.id, tr.stats.starttime))
                    
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
    ''' Make a single stream from the many streams'''

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
        # logging.info('end timestamp') 
        # logging.info(channel_stream[0].data[-1])
        # print('no of traces: ')
        # print(len(channel_stream))
        # channel_stream.plot(method='full')

        # merge the traces so that the timing is correct 


        # channel_stream.merge(method=1,fill_value=-1.0)
        # channel_stream.plot(method='full')


        # print(np.ma.count_masked(trace.data))
        # print(trace[1083])




        # print(trace)
        # fig = trace.plot(handle=True, show=False,method='full')
        # plt.ylim(-2,2)
        # plt.show()
        # print(trace[1083])
        # mask = np.ma.getmask(trace.data)
        # truemask = np.where(mask==1)
        # print(truemask)
        # # print(truemask[0])
        # # 

        # for t in range(0,len(truemask)):
        #     print(t)
        #     print(truemask[t])
        #     print(trace.data[truemask[t]])

# filled(np.nan)  
            
        # print(len(trace.data))
        # trace.plot(method='full')


        # merge and create a masked trace 
        channel_stream.merge(method=1)
        # channel_stream1 = channel_stream.copy()
        # channel_stream1 = channel_stream1.split()
        # channel_stream1.write('/Users/cnunn/Downloads/merged.MSEED', format='MSEED')

        if len(channel_stream) > 1:
            logging.info('Too many values in channel stream')
            raise Exception

        # fig = channel_stream.plot(handle=True, show=False,method='full')
        # # plt.ylim(-2,2)
        # plt.show()
        

        # # mask the trace which has gaps with a fill value
        # trace = channel_stream[0]
        # trace.data = trace.data.filled(fill_value=-1.0)

        # fig = trace.plot(handle=True, show=False,method='full')
        # # plt.ylim(-2,2)
        # plt.show()
        # XXXX
        # print('temp')
        # channel_stream.write('/Users/cnunn/Downloads/filled.MSEED', format='MSEED')
        # exit()

        merged_stream += channel_stream
        

    if config.fix_jump_error:
        for channel in ['DL4']:
            channel_stream = stream.select(channel=channel)
            channel_stream.merge(method=1)
            merged_stream += channel_stream

    if config.fix_clock_error:
 
        for channel in ['CLK']:
            channel_stream = stream.select(channel=channel)
            channel_stream.merge(method=1)
            merged_stream += channel_stream

        # now that the streams are merged, fix the software clock timing if 
        # possible 
        tr_CLK = merged_stream.select(channel='CLK')[0]

        # print('v temp')
        # tr_CLK.data[101] = 1
        # tr_CLK.data[102] = 1  
        # tr_CLK.data[103] = 1  

 
        # find where the clock flag is set 
        clk_set = np.where(tr_CLK.data==1)[0]
        if len(clk_set) > 0:
            date1 = tr_CLK.stats.starttime

            # print(date1)
            
            # mini_stream = Stream()
            # mini_stream.append(tr_CLK)
            # mini_stream.

            # logging.info('WARNING: Updating {} timestamps on the timing trace ATT due to a clock error'.format(len(clk_set)))   

            # get the ATT trace
            tr_ATT = merged_stream.select(channel='ATT')[0]

            # print('v temp')
            # tr_ATT.data[101] = tr_ATT.data[101] + 4
            # tr_ATT.data[102] = tr_ATT.data[102] + 4
            # tr_ATT.data[103] = tr_ATT.data[103] + 4
            # logging.info(UTCDateTime(tr_ATT.data[-1]))
            

            # get the mask (data gaps) for the ATT trace
            mask_ATT = np.ma.getmask(tr_ATT.data)
            # logging.info('ATT')
            # logging.info(mask_ATT[-10:])
            # logging.info(tr_ATT.data[-1])

        
            # # find where the software clock was set 
            # tr_CLK = merged_stream.select(channel='CLK')[0]
            # mask the ATT trace with the software clock 
            tr_CLK.data = np.ma.masked_where(tr_CLK.data==1,tr_CLK.data)                  # mask any values where the software clock has been set 
        
            # get the mask for the CLK trace
            mask_CLK = np.ma.getmask(tr_CLK.data)



            # logging.info('CLK')
            # logging.info(tr_CLK.data[-1])
            # logging.info(tr_CLK[-10:])
            # logging.info(mask_CLK[-10:])
            # 



            if tr_CLK.data[0] is ma.masked:
                logging.info('SEVERE: Initial record contains clock flag - {}.{}'.format(date1.year,date1.julday))   
                logging.info('temp stopping')
                exit()
                # # find first one that is OK
                # C_false = np.where( mask_CLK == False )
                # if len(C_false[0]) > 0: 
                #     first_good = C_false[0][0]
                #     mask_CLK[0:first_good] = False
                #     tr_CLK.data = ma.masked_array(tr_CLK.data, mask=mask_CLK)

            # logging.info(mask_CLK[-10:])
            # logging.info(type(mask_CLK[-1]))
            # logging.info(len(tr_CLK.data))
                
                # XXXX
            C_false = np.where( mask_CLK == False)
            if len(C_false[0]) > 0: 
                last_good = C_false[0][-1]
                first_good = C_false[0][0]
                # logging.info(last_good)
                # logging.info('Last date  {}'.format(UTCDateTime(tr_ATT.data[-1])))
                logging.info('First date without a clock error: {}'.format(UTCDateTime(tr_ATT.data[first_good])))
                logging.info('Last date without a clock error: {}'.format(UTCDateTime(tr_ATT.data[last_good])))

            if tr_CLK.data[-1] is ma.masked:
                logging.info('SEVERE: Last record contains clock flag - {}.{}'.format(date1.year,date1.julday)) 
                # retrieve the valid values 

                C_false = np.where( mask_CLK == False)
                if len(C_false[0]) > 0: 
                    last_good = C_false[0][-1]
                    first_good = C_false[0][0]
                    # logging.info(last_good)
                    # logging.info('Last date  {}'.format(UTCDateTime(tr_ATT.data[-1])))
                    logging.info('First date without a clock error: {}'.format(UTCDateTime(tr_ATT.data[first_good])))
                    logging.info('Last date without a clock error: {}'.format(UTCDateTime(tr_ATT.data[last_good])))
                    mask_CLK[last_good:] = False 
                    tr_CLK.data = ma.masked_array(tr_CLK.data, mask=mask_CLK)

                    # logging.info(mask_CLK.data[-10:])   
                    
            
            logging.info('WARNING: Updating {} timestamps on the timing trace ATT due to a clock error (total={})'.format(len(np.where( mask_CLK == False )[0]), len(tr_ATT) ))




            # [note that the mask doesn't include the empty values, but that doesn't 
            # matter]
            # print(tr_ATT)
            # tr_ATT.data = np.ma.masked_where(np.ma.getmask(mask_CLK), tr_ATT.data)
            tr_ATT.data= ma.masked_array(tr_ATT.data, mask=mask_CLK)
            # print(tr_ATT)

            # logging.info(tr_CLK.data[98:110])
            # logging.info(tr_ATT.data[98:110])
            # 


            # this does the work for interpolating the timing trace - it interpolates the time between good records
            data_series = pd.Series(tr_ATT.data)
            data_series.interpolate(method='linear', axis=0, limit=None, inplace=True, limit_direction=None, limit_area='inside', downcast=None)
            tr_ATT.data=data_series.to_numpy(dtype=np.float64)

            # logging.info(tr_ATT.data[98:110])
        
            # apply the original mask to the trace (to remove the data gaps)
            tr_ATT.data = ma.masked_array(tr_ATT.data, mask=mask_ATT)

            # logging.info(tr_ATT.data[98:110])

            # logging.info('end temp')
            # logging.info(tr_ATT.data[-1])    

        # no need to save the CLK trace 
        merged_stream.remove(tr_CLK)

    if config.fix_jump_error:
        # once the clock error is fixed, also try to fix the 400/800 error 

        tr_DL4 = merged_stream.select(channel='DL4')[0]

        # first make some checks on the overall gradient 

        tr_ATT = merged_stream.select(channel='ATT')[0]
        # tr_DL4 = merged_stream.select(channel='DL4')[0]

        start_timestamp = UTCDateTime(tr_ATT.data[0])
        # logging.info(tr_ATT.data[-1])
        end_timestamp = UTCDateTime(tr_ATT.data[-1])
        # XXXX

        time_diff = end_timestamp - start_timestamp

        overall_samples = len(tr_ATT)
        overall_delta4 = time_diff / (overall_samples - 1)

        overall_divergence = tr_ATT.times()[-1] - time_diff 

        # scale to 24 hours - usually not necessary because it's already 
        # about 24 hours 
        overall_divergence = overall_divergence * 86400/time_diff

        # logging.info('Numbers {} {} {} {} {} {} {} {}'.format(start_timestamp,end_timestamp,tr_ATT.times()[0],tr_ATT.times()[-1], overall_samples, time_diff, overall_delta4, overall_est_divergence))
        # 86400*(overall_delta4/DELTA*4)

        logging.info('Overall: {} Overall divergence: {:.02f}'.format(overall_delta4, overall_divergence))

        # station = tr_ATT.stats.station
        # segment_delta4  = df_end_frames.segment_delta4
        # segment_est_divergence = df_end_frames.segment_est_divergence
        # start_timestamp = df_end_frames.start_timestamp
        # end_timestamp = df_end_frames.end_timestamp
        # divergent_diff_found = False
        # for d4, est_d, t1, t2 in zip(segment_delta4, segment_est_divergence, start_timestamp, end_timestamp):
        #     if d4 is not pd.NA:
        #         divergence_diff = abs(overall_divergence - est_d)
        #         # diff = 100* abs(d4 - overall_delta4)/(DELTA*4)
        #         if divergence_diff > 0.2:
        #             divergent_diff_found = True
        #         logging.info('{} {} {} Segment: {} Divergence Difference: {} Segment estimated divergence: {}'.format(station, t1, t2, d4, divergence_diff, est_d))

        # potential_errors = config.potentially_suspect - config.rejected
        # if potential_errors > 0:
        #     logging.info('SEVERE: Large difference found for at least one of the gradients')
        # else:
        print('no checks before doing the correct shifts')
        correct_shifts(tr_ATT,tr_DL4)

        # no need to save the tr_DL4 trace 
        merged_stream.remove(tr_DL4)

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


def correct_shifts(orig_trace,delta4_trace):

    mask = np.ma.getmask(orig_trace.data)

    divergence_trace = orig_trace.copy() 
    relative_timing_trace(divergence_trace)

    dict_list = []

    # make an empty dataframe to store the results
    # df_results = pd.DataFrame(columns = [ 'station', 'start_time', 'end_time', 'correction'])

    station = orig_trace.stats.station 

    # make a dataframe from the original trace 
    df = pd.DataFrame(data=orig_trace.data, columns=['timestamp_orig'])

    # make a dataframe from the divergent trace (this nicely takes 
    # care of the masked data, if there are any)
    df_rel = pd.DataFrame(data=divergence_trace.data, columns=['divergence'])
    
    # add the second dataframe as a column to the first dataframe 
    df['divergence'] = df_rel.divergence

    # make a dataframe from the delta4 trace (this nicely takes 
    # care of the masked data, if there are any)
    df_dl4 = pd.DataFrame(data=delta4_trace.data, columns=['delta4'])

    # add the third dataframe as a column to the first dataframe 
    df['delta4'] = df_dl4.delta4

    # # get the details for the gradient for the segment
    # if start_clock_flag == 0 and end_clock_flag == 0:
    #     start_timestamp2 = orig_trace.data[-1]  
    #     end_timestamp2 = orig_trace.data[-1]
    #     time_diff = end_timestamp2 - end_timestamp1
    #     overall_delta4 = time_diff / len(df)
    # else:
    #     overall_delta4 = pd.NA

    # ZZZZ


    # make some checks
    
    
#     # print(df.divergence.isna().sum())
# 
    # get the details for the overall gradient
    x1 = divergence_trace.times()[0]
    x2 = divergence_trace.times()[-1]
    y1 = divergence_trace.data[0]
    y2 = divergence_trace.data[-1]

#         y = m * x + b
    m = (y1 - y2) / (x1 - x2)
    b = (x1 * y2 - x2 * y1) / (x1 - x2)

    # make a straight line
    straight_line = m*divergence_trace.times() + b
    trace_straight = divergence_trace.copy()
    trace_straight.data = straight_line
    df['straight_line'] = trace_straight.data

#     # get the details for the overall gradient
#     x1 = divergence_trace.times()[0]
#     x2 = divergence_trace.times()[-1]
#     y1 = divergence_trace.data[0]
#     y2 = divergence_trace.data[-1]
# 
#     # print(y2)

# #         y = m * x + b
#     m = (y1 - y2) / (x1 - x2)
#     b = (x1 * y2 - x2 * y1) / (x1 - x2)
# 
#     # make a straight line
#     straight_line = m*divergence_trace.times() + b
#     trace_straight = divergence_trace.copy()
#     trace_straight.data = straight_line
#     df['straight_line'] = trace_straight.data

    # print('divergence', df['divergence'].iloc[60000])
    # print(df.straight_line.iloc[60000])

    # need to find some way to fill the empty records


    # find first null value in timestamp_orig 
    # while True:
    #     idx_null = df[(df['timestamp_orig'].iloc[0:].isna()].index.tolist()
    #     if len(idx_null) > 0:
    #         idx_null_start = idx_null[0]
    #         print('start ', idx_null_start)
    #         # now find the next not null
    #         idx_null_idx_not_null = df['timestamp_orig'].iloc[idx_null_start:].first_valid_index()
    #         print('end ', idx_not_null)
    #         if idx_not_null is None:
    #             idx_not_null



    # idx_null = df[df['timestamp_orig'].isna()].index.tolist()
    # if len(idx_null) > 0
    # idx = idx_null[0]


    # ZZZZ

    df['last_valid_index'] = np.where((df['timestamp_orig'].isna()), pd.NA, df.index)

    # applying ffill() method to fill the missing values
    df['last_valid_index'].ffill(inplace=True)

    df['last_valid_divergence'] = np.where((df['timestamp_orig'].isna()), pd.NA, df.divergence)
    # applying ffill() method to fill the missing values
    df['last_valid_divergence'].ffill(inplace=True)


    df['last_valid_delta4'] = np.where((df['timestamp_orig'].isna()), pd.NA, df.delta4)
    # applying ffill() method to fill the missing values
    df['last_valid_delta4'].ffill(inplace=True)

# segment_est_divergence1 = time_diff - (len(df)-1) * DELTA *4

    

    # estimate the divergence (for when there is no valid timestamp)
    df['estimated_divergence'] = df['last_valid_divergence'] + ((df['last_valid_index'] - df.index) * (df.last_valid_delta4 - DELTA*4))
    # logging.info(df.dtypes)


    # idx_null = df[df['timestamp_orig'].isna()].index.tolist()
    # print(len(idx_null))
    # 
    # idx_null = df[df['valid_index'].isna()].index.tolist()
    # print(len(idx_null))    




# make a funky little algorithm to find the start null, end null, previous delta, project the gradient, find next null after some good values, project the gradient again. !!!1
            

    
    

    # # if the value is empty, use the estimated divergence, otherwise use the real divergence
    # df['divergence2'] = np.where(df.divergence.isna(), df.estimated_divergence, df.divergence)


    # calculate the jumps (using the estimated divergence column - same as real divergence when it is not an estimate)
    df['time_diff'] = df.estimated_divergence.diff(1)

    df['correction'] = 0.0
    df['possible_correction'] = pd.NA

    # find the difference between the divergence and straight line 
    # (the gradient when the start and ensd times are correct)
    df['divergence_diff'] = df.straight_line - df.divergence 

    idx_list =  df[((df['time_diff'].abs()) > 0.2)].index.tolist()
    if len(idx_list) > 0:
        df['timestamp'] = pd.to_datetime(df['timestamp_orig'], unit='s')
        logging.info('WARNING: Jumps found:')
        logging.info(df.iloc[idx_list].to_string())
    

    
    # jumps are estimated on the time difference from the divergence trace


    # first look for jumps - may be useful when we work out what else is going on 
    # we are looking for jumps around 0.4 s and 0.8 s
    df['jump'] = pd.NA

    df['jump'] = np.where((df.time_diff > 0.3) & (df.time_diff < 0.5), 0.4, df['jump'])
    df['jump'] = np.where((df.time_diff > 0.7) & (df.time_diff < 0.9), 0.8, df['jump'])
    df['jump'] = np.where((df.time_diff > 1.1) & (df.time_diff < 1.3), 1.2, df['jump'])
    df['jump'] = np.where((df.time_diff > 1.5) & (df.time_diff < 1.7), 1.6, df['jump'])
    df['jump'] = np.where((df.time_diff > 1.9) & (df.time_diff < 2.1), 2.0, df['jump'])
    
    df['jump'] = np.where((df.time_diff < -0.3) & (df.time_diff > -0.5), -0.4, df['jump'])
    df['jump'] = np.where((df.time_diff < -0.7) & (df.time_diff > -0.9), -0.8, df['jump'])
    df['jump'] = np.where((df.time_diff < -1.1) & (df.time_diff > -1.5), -1.2, df['jump'])
    df['jump'] = np.where((df.time_diff < -1.5) & (df.time_diff > -1.7), -1.6, df['jump'])
    df['jump'] = np.where((df.time_diff < -1.9) & (df.time_diff > -2.1), -2.0, df['jump'])

    idx_list =  df[(df['jump'].notna())].index.tolist()
    for idx in idx_list:
        jump = df.jump.iloc[idx]
        orig_timestamp = df.timestamp_orig.iloc[idx]
        # only the possible error, because sometimes the gradient is wrong
        possible_correction = round(df.divergence_diff.iloc[idx],1)
        if possible_correction in (-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0):
            df.at[idx,'possible_correction']= possible_correction
        # logging.info('WARNING: corrected tapehead jump from original start timestamp {} to original end timestamp {} correction {:.01f} s [auto]'.format(start_timetamp,end_timestamp,correction))
            logging.info('WARNING: corrected tapehead jump from original start timestamp {} difference {:.01f} s, possible error {:.01f} s [auto]'.format(UTCDateTime(orig_timestamp),jump,possible_correction))



    # df['jump'] = np.where((df.time_diff > 0.3) | (df.time_diff < -0.3), df.time_diff, 0.0)

    # if divergent_diff_found == False:    
    #     df['correction'] = np.where((df.divergence_diff > 0.3) & (df.divergence_diff < 0.5), 0.4, df['correction'])
    #     df['correction'] = np.where((df.divergence_diff > 0.7) & (df.divergence_diff < 0.9), 0.8, df['correction'])
    #     df['correction'] = np.where((df.divergence_diff > 1.1) & (df.divergence_diff < 1.3), 1.2, df['correction'])
    #     df['correction'] = np.where((df.divergence_diff > 1.5) & (df.divergence_diff < 1.7), 1.6, df['correction'])
    #     df['correction'] = np.where((df.divergence_diff > 1.9) & (df.divergence_diff < 2.1), 2.0, df['correction'])
    # 
    #     df['correction'] = np.where((df.divergence_diff < -0.3) & (df.divergence_diff > -0.5), -0.4, df['correction'])
    #     df['correction'] = np.where((df.divergence_diff < -0.7) & (df.divergence_diff > -0.9), -0.8, df['correction'])
    #     df['correction'] = np.where((df.divergence_diff < -1.1) & (df.divergence_diff > -1.5), -1.2, df['correction'])
    #     df['correction'] = np.where((df.divergence_diff < -1.5) & (df.divergence_diff > -1.7), -1.6, df['correction'])
    #     df['correction'] = np.where((df.divergence_diff < -1.9) & (df.divergence_diff > -2.1), -2.0, df['correction'])

    # potentially useless calculations: 


    # df['correction'] = pd.NA
    # df.at[0,'correction']= 0.0
    

    # 

            #     divergence_diff = df.divergence_diff.iloc[idx]
            #     if jump in (-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0):
            #     df.at[idx,'correction']= jump

    # applying ffill() method to fill the missing values
    df.possible_correction.ffill(inplace=True)

    df['divergence_diff'] = df['divergence_diff'].round(1)
    df['correction'] = np.where((df.possible_correction == df.divergence_diff), df['possible_correction'], 0.0)
    
    # 

    

    # idx = 64987
    # logging.info(df.iloc[idx - 10: idx + 11].to_string())
    # logging.info(df.dtypes)
    correction_sum =  ((df['correction'].lt(0)) | (df['correction'].gt(0))).sum()
    if correction_sum > 0:
        logging.info('WARNING: {} timestamps corrected due to 400/800 ms tape head error'.format(correction_sum))

    df['timestamp_corr'] = df.timestamp_orig - df['correction']

    timestamp_corr = df[['timestamp_corr']].to_numpy(dtype='float64',na_value=None)
    # abs_times = np.repeat(abs_times,4)
    timestamp_corr = timestamp_corr.flatten()

    orig_trace.data = timestamp_corr

    # put the original mask back
    orig_trace.data = np.ma.masked_array(orig_trace.data, mask)

def estimate_timestamp(current_index, last_valid_index, last_valid_timestamp, last_valid_delta4):
   return last_valid_timestamp + (current_index-last_valid_index) * last_valid_delta4

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
    starttime0 = df_list1['starttime']
    endtime = df_list1['endtime']
    orig_station = df_list1['orig_station']
    # orig_ground_station = df_list1[3]
    corr_ground_station = df_list1['corr_ground_station']
    df = df_list1['df']
    attempt_merge = df_list1['attempt_merge']
    segment_delta4 = df_list1['segment_delta4']
    segment_est_divergence = df_list1['segment_est_divergence']
    gzip_filename = df_list1['gzip_filename']

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
        
        time_int = round((starttime0 - midnight)/segment_delta4)-1

        est_start_time = starttime0 - time_int *segment_delta4

        # start the index at zero 
        df.time_index = np.arange(len(df))
        
        sample_time0 = est_start_time + time_int *DELTA*4

        # logging.info(sample_time0)
        # logging.info(est_start_time)


    

    # projection = 24*3600 * (end_delta4_final-(DELTA*4))/DELTA*4

    # loc = len(df_end_frames) - 1


    logging.info('INFO: Zero index record. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence={:.02f} Reject={} File={}'.format(index_no, orig_station, corr_ground_station, starttime0, sample_time0, starttime0-sample_time0, segment_est_divergence,df_list1['reject'],gzip_filename))

    calculate_last_valid(df,index_no)

    return starttime0, sample_time0


def later_records(index_no,sample_time0, sorted_df_list):

    # get the information from the sorted_df_list
    df_list1 = sorted_df_list[index_no]
    # starttime = df_list1[0]
    # endtime = df_list1[1]
    orig_station = df_list1['orig_station']
    orig_ground_station = df_list1['orig_ground_station']
    corr_ground_station = df_list1['corr_ground_station']
    attempt_merge = df_list1['attempt_merge']
    reject = df_list1['reject'] 
    segment_est_divergence = df_list1['segment_est_divergence']
    gzip_filename = df_list1['gzip_filename']

    # find the current dataframe
    df = df_list1['df']
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
    start_delta4 = df.delta4.iloc[0]

    # find the previous dataframe
    df_list_prev = sorted_df_list[index_no-1]
    # prev_starttime = df_list_prev[0]
    # prev_endtime = df_list_prev[1]
    
    prev_df = df_list_prev['df']

    # is this a clock merge or a normal one?
    clock_flag = df.clock_flag[0]

    # if this contains the software clock, we need to use the last entry
    # with the software clock - otherwise we will impose gaps 
    # if not, we need to use the correct timestamps

    # if there's a clock flag, use the clock flag record
    if clock_flag == 1 and not (pd.isna(config.clock_end_time_index)):
        end_frame = config.clock_end_frame
        end_timestamp = config.clock_end_timestamp
        end_time_index = config.clock_end_time_index
        end_delta4 = config.clock_end_delta4 
    else: 
        # if not found, use the most recent valid record 
        end_frame = config.valid_end_frame
        end_timestamp = config.valid_end_timestamp
        end_time_index = config.valid_end_time_index
        end_delta4 = config.valid_end_delta4 

    timestamp_time_diff = (start_timestamp - end_timestamp).total_seconds()
    
    starttime= UTCDateTime(start_timestamp)

    found_index = []
    method=''
    
    # logging.info('{} Original starttime (data time)={} 1st frame={} Original endtime={} Last frame={} ground_station={}'.format(index_no, starttime, frame, endtime, df.frame.iloc[-1],corr_ground_station))

    overlap_found = False
    include_trace = False

    # print('changed this temporarily - was 2 s')
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
                if reject == 'MAYBE':
                    reject = 'YES'
                else: 
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
            logging.info('INFO: Overlapping traces found. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence={:.02f} Reject={} File={}'.format(index_no, orig_station, orig_ground_station, UTCDateTime(starttime), sample_time, UTCDateTime(starttime)-sample_time,segment_est_divergence,df_list1['reject'],gzip_filename))

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

        gradient = (start_delta4 + end_delta4)/2

        frame_gap_correct, actual_frame_gap, absolute_error, percent_error= loose_frame_diff(last_frame=end_frame,new_frame=frame,gap=timestamp_time_diff,delta4_1=end_delta4,delta4_2=start_delta4)
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
            logging.info('INFO: Time gap with correct number of missing samples found. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence {:.02f} Reject={} File={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time,segment_est_divergence,df_list1['reject'],gzip_filename))
        else: 
            if reject in ('MAYBE', 'YES'):
                reject = 'YES'
                logging.info('INFO: Frame Gap not correct. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence {:.02f} Reject={} File={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time,segment_est_divergence, df_list1['reject'],gzip_filename))
            else:
                if overlap_found: 
                    logging.info('SEVERE: Overlap found and time gap does not agree with frame gap. Overwriting part of trace. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence {:.02f} Reject={} File={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time,segment_est_divergence,df_list1['reject'],gzip_filename))
                    attempt_merge = False
                else: 
                    logging.info('SEVERE: Time gap does not agree with frame gap. Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence {:.02f} Reject={} File={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time,segment_est_divergence,df_list1['reject'],gzip_filename))

    # TODO how do we sort later - should the timing be updated? 
    df_list1['attempt_merge'] = attempt_merge

    if reject == 'YES':
        config.rejected += 1
        logging.info('INFO: Rejecting Trace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence {:.02f} Reject={} File={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time,segment_est_divergence, df_list1['reject'], gzip_filename))
        df.drop(df.index, inplace=True)
        # logging.info(df.to_string())

    calculate_last_valid(df,index_no)


def calculate_segment_delta4(df):
    
    start_clock_flag = df.clock_flag.iloc[0]
    end_clock_flag = df.clock_flag.iloc[-1]
    
    # get the details for the gradient for the segment
    if start_clock_flag == 0 and end_clock_flag == 0:
        start_timestamp = UTCDateTime(df.corr_timestamp.iloc[0])   
        end_timestamp = UTCDateTime(df.corr_timestamp.iloc[-1])
        time_diff = end_timestamp - start_timestamp
        segment_delta4 = time_diff / (len(df)-1)

        segment_est_divergence1 =  (len(df)-1) * DELTA *4 - time_diff
        # now scale for 24 hours 
        segment_est_divergence = segment_est_divergence1 * 86400 / time_diff
    else:
        segment_delta4 = pd.NA
        segment_est_divergence = pd.NA

    return segment_delta4, segment_est_divergence

def calculate_last_valid(df,index_no):

    if index_no == 0:
            config.valid_end_frame = pd.NA
            config.valid_end_timestamp = pd.NA
            config.valid_end_time_index = pd.NA
            config.valid_end_delta4 = pd.NA

            config.clock_end_frame = pd.NA
            config.clock_end_timestamp = pd.NA
            config.clock_end_time_index = pd.NA
            config.clock_end_delta4 = pd.NA        

    if len(df) > 0: 
        end_clock_flag = df.clock_flag.iloc[-1]
        if end_clock_flag == 0:

            config.valid_end_frame = df.frame.iloc[-1]
            config.valid_end_timestamp = df.corr_timestamp.iloc[-1]
            config.valid_end_time_index = df.time_index.iloc[-1]
            config.valid_end_delta4 = df.delta4.iloc[-1]

            config.clock_end_frame = pd.NA
            config.clock_end_timestamp = pd.NA
            config.clock_end_time_index = pd.NA
            config.clock_end_delta4 = pd.NA

        else: 

            config.clock_end_frame = df.frame.iloc[-1]
            config.clock_end_timestamp = df.corr_timestamp.iloc[-1]
            config.clock_end_time_index = df.time_index.iloc[-1]
            config.clock_end_delta4 = df.delta4.iloc[-1]

    # last_timestamp = df.corr_timestamp.iloc[-1]
    # last_time
    # 
    # 
    # # find the last record in the last timeseries
    # 
    # 
    # 
    # end_clock_flag = df.clock_flag.iloc[-1]
    # 
    # start_timestamp = df.corr_timestamp.iloc[0]
    # 
    # start_clock_flag = df.clock_flag.iloc[0]
    # 
    # # if the clock flag was set, then also get the last good record, 
    # # if there was one
    # if end_clock_flag == 1:
    #     # find the last record where the clock flag wasn't set 
    #     idx_list = df[df.clock_flag == 0].index.tolist()
    #     if len(idx_list) > 0:
    #         last_good_clock_flag  = idx_list[-1]
    #         # find the last record for the good clock flag
    #         end_frame1 = df.frame.iloc[last_good_clock_flag]
    #         end_timestamp1 = df.corr_timestamp.iloc[last_good_clock_flag]
    #         end_time_index1 = df.time_index.iloc[last_good_clock_flag]
    #         end_clock_flag1 = df.clock_flag.iloc[last_good_clock_flag]
    #         end_delta4_1 = df.delta4.iloc[last_good_clock_flag]
    #         df_end_frames1 = {'end_frame' : end_frame1, 'end_timestamp' :end_timestamp1, 'end_time_index' : end_time_index1, 'end_clock_flag' : end_clock_flag1, 'end_delta4' : end_delta4_1}
    #         df_end_frames = df_end_frames.append(df_end_frames1,ignore_index = True)
    # 
    # 
    # df_end_frames2 = {'end_frame' : end_frame, 'end_timestamp' :end_timestamp, 'end_time_index' : end_time_index, 'end_clock_flag' : end_clock_flag, 'end_delta4' : end_delta4,  'start_timestamp' : start_timestamp}
    # df_end_frames = df_end_frames.append(df_end_frames2,ignore_index = True)
    # 
    # return df_end_frames

# def calculate_end_frames(df,df_end_frames=None):
# 
#     # Start!!
#     if df_end_frames is None:
#         # make a dataframe to store the last records 
#         df_end_frames = pd.DataFrame(data=None, columns=['end_frame','end_timestamp','end_time_index','end_clock_flag'])
# 
#     # find the last record in the last timeseries
#     end_frame = df.frame.iloc[-1]
#     end_timestamp = df.corr_timestamp.iloc[-1]
#     end_time_index = df.time_index.iloc[-1]
#     end_clock_flag = df.clock_flag.iloc[-1]
#     end_delta4 = df.delta4.iloc[-1]
#     start_timestamp = df.corr_timestamp.iloc[0]
# 
#     start_clock_flag = df.clock_flag.iloc[0]
# 
#     # if the clock flag was set, then also get the last good record, 
#     # if there was one
#     if end_clock_flag == 1:
#         # find the last record where the clock flag wasn't set 
#         idx_list = df[df.clock_flag == 0].index.tolist()
#         if len(idx_list) > 0:
#             last_good_clock_flag  = idx_list[-1]
#             # find the last record for the good clock flag
#             end_frame1 = df.frame.iloc[last_good_clock_flag]
#             end_timestamp1 = df.corr_timestamp.iloc[last_good_clock_flag]
#             end_time_index1 = df.time_index.iloc[last_good_clock_flag]
#             end_clock_flag1 = df.clock_flag.iloc[last_good_clock_flag]
#             end_delta4_1 = df.delta4.iloc[last_good_clock_flag]
#             df_end_frames1 = {'end_frame' : end_frame1, 'end_timestamp' :end_timestamp1, 'end_time_index' : end_time_index1, 'end_clock_flag' : end_clock_flag1, 'end_delta4' : end_delta4_1}
#             df_end_frames = df_end_frames.append(df_end_frames1,ignore_index = True)
# 
# 
#     df_end_frames2 = {'end_frame' : end_frame, 'end_timestamp' :end_timestamp, 'end_time_index' : end_time_index, 'end_clock_flag' : end_clock_flag, 'end_delta4' : end_delta4,  'start_timestamp' : start_timestamp}
#     df_end_frames = df_end_frames.append(df_end_frames2,ignore_index = True)
# 
#     return df_end_frames


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
                    if station == 'S14' and channel == 'MHZ' and starttime > UTCDateTime('1972-03-12T00:00') and starttime < UTCDateTime('1976-11-17T00:00:00.000000Z'):
                        # we are not saving the MHZ track during this time frame. 
                        continue


                    if len(channel_stream) > 0:
                        filename = filename_from_trace(channel_stream[0],directory,lower=False)
                        filename = os.path.join(directory,filename)

                        if channel == 'ATT':
                            # mask the ATT trace which has gaps with a fill value
                            trace = channel_stream[0]
                            trace.data = trace.data.filled(fill_value=-1.0)
                        
                        # logging.info('Temporarily not writing {}'.format(filename))
                        # split the streams in order to save them
                        # channel_stream = channel_stream.split()
                        
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

        range = (gap*MAX_PERCENT_ERROR/100) + TOLERANCE
        # if range is greater than 20 s, we can't make the estimate 
        if range > 20:
            range = 20

        max_estimate = gap + range
        min_estimate = gap - range

        logging.debug('gap {} max estimate {}, min estimate {}'.format(gap, min_estimate,max_estimate))
        
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

    logging.debug('loose_frame_diff {} {} {} absolute error {} percent error{} frame_gap_correct {}'.format(last_frame,new_frame,gap, absolute_error, percent_error, frame_gap_correct))

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

    # print('v temp test for not displaying when there are gaps')
    # df.at[1000:2000,'orig_no'] = pd.NA


    # logging.info(df.head().to_string())
    # logging.info(df.tail().to_string())

    # XXXX

    # check for any nulls at the beginning and remove    
    if pd.isna(df['orig_no'].iloc[0]) or pd.isna(df['orig_no'].iloc[1]):
        idx_list_nonull =  df[(df['orig_no'].notna())].index.tolist()
        if pd.isna(df['orig_no'].iloc[0]) == True:
            first_good_idx = idx_list_nonull[0]
        else: 
            first_good_idx = idx_list_nonull[1]
        df.drop(df.index[0:first_good_idx], inplace=True)
        df.reset_index(inplace=True,drop=True)
    
    # check for any nulls at the end and remove
    if pd.isna(df['orig_no'].iloc[-1]) or pd.isna(df['orig_no'].iloc[-2]):
        idx_list_nonull =  df[(df['orig_no'].notna())].index.tolist()
        if pd.isna(df['orig_no'].iloc[-1]):
            last_good_idx = idx_list_nonull[-1]
        else:
            last_good_idx = idx_list_nonull[-2]
        df.drop(df.index[last_good_idx+1:len(df)], inplace=True)
    
    # logging.info(df.head().to_string())
    # logging.info(df.tail().to_string())

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

    # make a mask for corr_timestamp 
    mask = df[['orig_no']].to_numpy(dtype='int32',na_value=INVALID)
    mask = mask.flatten()
    # print(mask[0:5])
    mask = ma.masked_equal(mask, INVALID)
    mask = ma.getmask(mask)
    # print(mask[0:5])

    # print(abs_times[0:5])

    # XXXX
    if config.exclude_masked_sections:
        # mask the trace during data gaps 
        abs_times = ma.masked_array(abs_times, mask=mask)


    delta4 = df[['delta4']].to_numpy(dtype='float64',na_value=None)
    # abs_times = np.repeat(abs_times,4)
    delta4 = delta4.flatten()


    ground_station = str(df['corr_ground_station'].iloc[0])
    
    station = df['orig_station'].iloc[0]

    # set the starttime based on the estimate using the DELTA, not
    # using the the first time in the timestamp 

    sample_time = sample_time0 + df.time_index.iloc[0] * DELTA * 4
    # adjustment_time = UTCDateTime(abs_times[0]) - starttime_adjust
    # logging.info('clock starttime={} Adjusting time {}s'.format(UTCDateTime(abs_times[0]),adjustment_time))

    # TODO - check that we don't mess up the time differences here!!!!

    clock_flag = df[['clock_flag']].to_numpy(dtype='int32',na_value=INVALID)
    clock_flag = clock_flag.flatten()


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

    tr_CLK = Trace(data=clock_flag)
    tr_CLK.stats.delta = DELTA*4
    tr_CLK.stats.network = NETWORK
    tr_CLK.stats.station = station
    tr_CLK.stats.ground_station = ground_station
    tr_CLK.stats.channel = 'CLK'
    tr_CLK.stats.starttime = sample_time
    tr_CLK.stats.location = ''

    tr_DL4 = Trace(data=delta4)
    tr_DL4.stats.delta = DELTA*4
    tr_DL4.stats.network = NETWORK
    tr_DL4.stats.station = station
    tr_DL4.stats.ground_station = ground_station
    tr_DL4.stats.channel = 'DL4'
    tr_DL4.stats.starttime = sample_time
    tr_DL4.stats.location = ''

    # append the trace
    stream.append(tr_MH1)
    stream.append(tr_MH2)
    stream.append(tr_MHZ)
    stream.append(tr_ATT)
    stream.append(tr_CLK)
    stream.append(tr_DL4)

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