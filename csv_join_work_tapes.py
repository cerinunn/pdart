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
from copy import deepcopy

from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, Trace, Stats, read

# from pdart.save_24_hours import update_starttimes
import pdart.config as config
from pdart.util import relative_timing_trace
# from pdart.csv_import_work_tapes import find_output_dir, make_output_dir, make_filelist
import matplotlib  
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
global df_manual_grab_before
df_manual_grab_before = None
global df_manual_grab_after
df_manual_grab_after = None

global manual_only
manual_only = False

# global DELTA
DELTA = 0.1509433962
# INTERVAL = '603774us'

EXACT_N = 86400/DELTA/4

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
    manual_grab_before=None,
    manual_grab_after=None,
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

    global df_manual_grab_before
    if manual_grab_before is not None:
        df_manual_grab_before = pd.read_csv(manual_grab_before, dtype=str, comment='#')
        df_manual_grab_before.year = df_manual_grab_before.year.astype('int')
        df_manual_grab_before.julday = df_manual_grab_before.julday.astype('int')
        df_manual_grab_before.station = df_manual_grab_before.station.str.strip()        

    global df_manual_grab_after
    if manual_grab_after is not None:
        df_manual_grab_after = pd.read_csv(manual_grab_after, dtype=str, comment='#')
        df_manual_grab_after.year = df_manual_grab_after.year.astype('int')
        df_manual_grab_after.julday = df_manual_grab_after.julday.astype('int')
        df_manual_grab_after.station = df_manual_grab_after.station.str.strip()


    # TODO -fix these date ranges - they are not quite right 
    # if no filenames have been passed, then look for them 
    for year in range(year_start,year_end+1):
        for day in range(day_start,day_end+1):
            if config.combine_ground_stations == True:
                if manual_only == False: 
                    log_filename = 'join.{}.{:03}.log'.format(year,day)
                else:
                    sta = stations[0]
                    log_filename = 'join_manual.{}.{:03}.{}.log'.format(year,day,sta)
            else:
                log_filename = 'ground_stations.{}.{}.log'.format(year,day)    
            log_filename = os.path.join(log_dir,log_filename)
            logging.basicConfig(filename=log_filename, filemode='w', 
              level=logging_level,format='%(message)s')

            print('Processing {}.{}'.format(year,day))

            logging.info('############################################')
            logging.info('Processing {}.{}'.format(year,day))
            logging.info('############################################')

            logging.info('#########')   

            logging.info('config.combine_ground_stations={}'.format(config.combine_ground_stations))  
            logging.info('config.clean_spikes={}'.format(config.clean_spikes))  
            logging.info('config.view_corrected_traces={}'.format(config.view_corrected_traces))  
            logging.info('config.fix_clock_error={}'.format(config.fix_clock_error))  
            logging.info('config.fix_jump_error={}'.format(config.fix_jump_error))
            logging.info('config.exclude_masked_sections={}'.format(config.exclude_masked_sections))
                                
            for station in stations:
                df_list = []
                days = [day]
                # check whether to get previous or next day
                # make a manual_exclusion type thing !!
                grab_after = manual_after(year,day,station)
                if grab_after:
                    days.append(day+1)
                grab_before = manual_before(year,day,station)
                if grab_before:
                    days.append(day-1)
                
                for day1 in days: 
                    wildcard_filename = '{}.*.*.{}.*.{}.{:03}.*.csv.gz'.format(wildcard_style,station,year,day1)
                    print('wildcard filename ', processed_dir, wildcard_filename)
                    # read in each of the csv files for the station
                    # for filename in glob.glob('/Users/cnunn/python_packages/pdart/examples/test.csv'):
                    for filename in glob.glob(os.path.join(processed_dir,wildcard_filename)):

                        # if os.path.basename(filename) not in ('wtn.1.3.S12.1.1976.061.13_29_58_945000.csv.gz'):
                        #     continue

                        exclude, comment = manual_exclusion(filename)


                        if ('dropped' not in filename) and exclude == False:
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

                        else:
                            # read in the dropped/excluded csv file and make a dataframe
                            df = pd.read_csv(filename, dtype=str, comment='#')

                            if exclude:
                                exclude_type = 'Excluded'
                                start_timestamp = UTCDateTime(df['corr_timestamp'].iloc[0])
                                end_timestamp = UTCDateTime(df['corr_timestamp'].iloc[-1])
                                orig_station = df['orig_station'].iloc[0]
                            else:
                                exclude_type = 'Dropped'
                                comment = ''
                                start_timestamp = UTCDateTime(df['corr_timestamp'].iloc[0])
                                end_timestamp = UTCDateTime(df['corr_timestamp'].iloc[-1])
                                orig_station = df['orig_station'].iloc[0]

                            logging.info('CORRECTION,{},{},{:03},{},,{},{},{},{}'.format(orig_station, start_timestamp.year,start_timestamp.julday, exclude_type, start_timestamp, end_timestamp, os.path.basename(filename),comment))    

                            continue
        

                logging.info('#########')

                if len(df_list) > 0:        
                    try: 
                        # process the list of dataframes
                        process_list(df_list,join_dir,year,day)
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

    if len(df) < 1:
        logging.debug('WARNING: No valid record(s), not importing'.format(len(df)))
        return df_list

    # leap second correction 
    # df = leap_correction(df,gzip_filename)

    # check the clock flags prior to the manual correction
    clock_flag_total = (df['clock_flag'] == 1).sum()
    df_total = (~pd.isna(df['orig_no'])).sum()
    if df_total > 0:
        clock_flag_perc = clock_flag_total/df_total * 100 
    else:
        clock_flag_perc = 0

    df = manual_correction(df,gzip_filename)

    # df_end_frames = calculate_end_frames(df,df_end_frames=df_end_frames)

    # calculate the delta4 for the segment 
    segment_delta4, segment_est_divergence, segment_est_divergence_incl_clock = calculate_segment_delta4(df)

    # calculate the mean for the section
    mh1_mean, mh2_mean, mhz_mean = calculate_segment_mean(df)


    # get some details about the data 
    starttime = UTCDateTime(df.corr_timestamp.iloc[0])
    endtime = UTCDateTime(df.corr_timestamp.iloc[-1])
    start_clock_flag = df.clock_flag.iloc[0]
    end_clock_flag = df.clock_flag.iloc[-1]
    orig_station = df.orig_station.iloc[0]
    orig_ground_station = df.orig_ground_station.iloc[0]
    corr_ground_station = df.corr_ground_station.iloc[0]

    # display a message about the clock flags set originally on the trace
    # ignore just a few because probably not significant
    #
    if clock_flag_perc > 1:
        logging.info('CORRECTION,{},{},{:03},Interpolation (clock flags set on the trace),{:.01f}%,{},{},{},'.format(orig_station, starttime.year,starttime.julday, clock_flag_perc, starttime, endtime, os.path.basename(gzip_filename)))    

    attempt_merge = True

    df_dict = {'starttime' : starttime, 'endtime' : endtime, 
      'orig_station' : orig_station , 'orig_ground_station' : orig_ground_station, 
      'corr_ground_station' : corr_ground_station, 'df' : df, 
      'attempt_merge' : attempt_merge, 'segment_delta4' : segment_delta4,
      'segment_est_divergence' : segment_est_divergence, 'gzip_filename' : gzip_filename,
      'segment_est_divergence_incl_clock' : segment_est_divergence_incl_clock,
      'mh1_mean' : mh1_mean,
      'mh2_mean' : mh2_mean,
      'mhz_mean' : mhz_mean
        }
    # make a list of the dataframes
    df_list.append(df_dict)

    return df_list

def manual_after(year,julday,station):
    grab_after = False
    if df_manual_grab_after is not None: 
        if len(df_manual_grab_after) > 0:
            idx_list = df_manual_grab_after[(df_manual_grab_after['year'] == year) & (df_manual_grab_after['julday'] == julday)  & (df_manual_grab_after['station'] == station) ].index.tolist()
            # idx_list = df_manual_grab_after[(df_manual_grab_after['year'] == year)  & (df_manual_grab_after['julday'] == julday)].index.tolist()
            # idx_list = df_manual_grab_after[(df_manual_grab_after['station'] == station)].index.tolist()
# (df_manual_grab_after['station'] == station) 
            if len(idx_list) > 0: 
                grab_after = True

    return grab_after

def manual_before(year,julday,station):
    grab_before = False
    if df_manual_grab_before is not None: 
        if len(df_manual_grab_before) > 0:
            idx_list = df_manual_grab_before[(df_manual_grab_before['year'] == year) & (df_manual_grab_before['julday'] == julday) & (df_manual_grab_before['station'] == station) ].index.tolist()
            if len(idx_list) > 0: 
                grab_before = True
    return grab_before

def manual_exclusion(gzip_filename):
    exclude = False
    comment = ''
    if df_manual_exclude is not None: 
        if len(df_manual_exclude) > 0:
            idx_list = df_manual_exclude[df_manual_exclude['filename'] == os.path.basename(gzip_filename)].index.tolist()
            if len(idx_list) > 0: 
                exclude = True
                comment = df_manual_exclude['comment'].iloc[idx_list[0]]
    
    return exclude, comment

def manual_correction(df,gzip_filename):
    if df_manual_jump_correction is not None: 
        
        if len(df_manual_jump_correction) > 0: 
            idx_list = df_manual_jump_correction[df_manual_jump_correction['filename'] == os.path.basename(gzip_filename)].index.tolist()
            if len(idx_list) > 0:
                
                idx = idx_list[0]

                correction = df_manual_jump_correction.correction.iloc[idx]

                start_timestamp = UTCDateTime(df['corr_timestamp'].iloc[0])
                end_timestamp = UTCDateTime(df['corr_timestamp'].iloc[-1])
                orig_station = df['orig_station'].iloc[0]
        
                logging.info('CORRECTION,{},{},{:03},Constant,{:.01f} s,{},{},{},'.format(orig_station, start_timestamp.year,start_timestamp.julday, correction, start_timestamp, end_timestamp, os.path.basename(gzip_filename)))    

                df['corr_timestamp'] = df['corr_timestamp'] - pd.Timedelta(correction, unit='seconds')
                df['orig_timestamp'] = df['orig_timestamp'] - pd.Timedelta(correction, unit='seconds')


    if df_manual_clock_correction is not None: 
        if len(df_manual_clock_correction) > 0:
            idx_list = df_manual_clock_correction[df_manual_clock_correction['filename'] == os.path.basename(gzip_filename)].index.tolist()
            if len(idx_list) > 0: 
                logging.info('Clock flag set for {}'.format(gzip_filename))
                df.clock_flag = 1
                start_timestamp = UTCDateTime(df['corr_timestamp'].iloc[0])
                end_timestamp = UTCDateTime(df['corr_timestamp'].iloc[-1])
                orig_station = df['orig_station'].iloc[0]
                logging.info('CORRECTION,{},{},{:03},Interpolation,,{},{},{},'.format(orig_station, start_timestamp.year,start_timestamp.julday, start_timestamp, end_timestamp, os.path.basename(gzip_filename)))    

    return df

# def leap_correction(df, filename):
    # moved to manual correction

    # make a correction for the leap seconds
    # The UTC added one second at the end of the year from 1972 to 1979.  
    # Please look at the wikipedia web site: https://en.wikipedia.org/wiki/Leap_second.
    # We remove one second because the timing was based on the PREVIOUS day's 
    # timing, and is so is recorded 1 second too late 

    # 
    # original = 1,"1975-01-01T00:00:10.703000Z",0,"S12",0,2,532,531,520,533,532,521,534,534,521,534,534,522,7,
    # corrected = 1975-01-01T00:00:09.703000Z

    # In this example, the correct timestamp is 9 s after midnight, but 
    # it was recorded in the previous timeframe as 10 s after midnight, since
    # the leap second is missing. 
#     if os.path.basename(filename) in [
# 
#         'pse.a12.5.241.S12.5.1972.183.00_00_00_118000.csv.gz',
#         'pse.a12.5.241.S12.13.1972.183.05_59_57_593000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.00_00_00_232000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.00_16_50_278000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.00_25_08_961000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.00_32_25_460000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.00_43_18_096000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.00_52_19_645000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.01_02_21_567000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.01_13_14_203000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.01_23_39_066000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.01_39_18_476000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.01_41_22_241000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.01_47_35_952000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.01_55_14_789000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.02_00_22_089000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.02_09_55_636000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.02_14_20_675000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.02_20_28_348000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.02_27_46_659000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.02_32_40_073000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.02_38_44_125000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.02_42_29_317000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.02_51_45_959000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_02_42_217000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_05_05_302000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_09_46_643000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_11_57_048000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_17_40_573000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_24_52_846000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_28_31_397000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_34_01_639000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_37_31_739000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_39_35_504000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_43_11_037000.csv.gz',
#         'pse.a14.3.155.S14.13.1972.183.03_45_42_574000.csv.gz',
#         'pse.a15.2.156.S15.3.1972.183.00_00_00_311000.csv.gz',
#         'pse.a16.1.71.S16.3.1972.183.00_00_00_165000.csv.gz',
# 
#         'pse.a12.6.136.S12.5.1973.001.00_00_09_330000.csv.gz',
#         'pse.a14.4.166.S14.5.1973.001.00_00_32_799000.csv.gz',
#         'pse.a15.3.162.S15.3.1973.001.03_59_06_810000.csv.gz',
#         'pse.a15.3.162.S15.3.1973.001.00_00_09_103000.csv.gz',
#         'pse.a16.2.85.S16.3.1973.001.00_00_10_037000.csv.gz',
# 
#         'pse.a12.7.211.S12.1.1974.001.00_00_06_476000.csv.gz',
#         'pse.a14.7.1.S14.1.1974.001.00_00_11_038000.csv.gz',
#         'pse.a15.6.1.S15.1.1974.001.00_00_10_969000.csv.gz',
#         'pse.a16.4.94.S16.1.1974.001.00_00_11_168000.csv.gz',
# 
#         'pse.a12.9.1.S12.2.1975.001.00_00_10_703000.csv.gz',
#         'pse.a12.9.1.S12.2.1975.001.03_53_57_272000.csv.gz',
#         'pse.a14.9.7.S14.2.1975.001.00_00_23_815000.csv.gz',
#         'pse.a15.8.7.S15.2.1975.001.00_00_10_932000.csv.gz',
#         'pse.a16.6.100.S16.2.1975.001.00_00_10_238000.csv.gz',
# 
#         'pse.a12.10.78.S12.5.1976.001.00_00_09_486000.csv.gz',
#         'pse.a14.11.20.S14.5.1976.001.00_00_12_571000.csv.gz',
#         'pse.a15.10.23.S15.5.1976.001.00_00_07_307000.csv.gz',
#         'pse.a16.8.116.S16.5.1976.001.00_00_12_276000.csv.gz',
# 
#         'wtn.11.23.S12.1.1977.001.00_12_50_785000.csv.gz',
#         'wtn.11.23.S14.1.1977.001.00_12_50_927000.csv.gz',
#         'wtn.11.23.S15.1.1977.001.00_12_50_923000.csv.gz',
#         'wtn.11.23.S16.1.1977.001.00_12_50_863000.csv.gz'
# 
# ]:
#         df['corr_timestamp'] = df['corr_timestamp'] + pd.Timedelta(-1, unit='seconds')
#         df['orig_timestamp'] = df['orig_timestamp'] + pd.Timedelta(-1, unit='seconds')
# 
#     return df 


def process_list(df_list,join_dir,year,julday):

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

    # sort by starttime
    # sorted_df_list = sorted(df_list, key=operator.itemgetter(0))
    sorted_df_list =  (sorted(df_list, key = lambda i: i['starttime']))

    continue_idx = None

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
                sorted_df_list[i]['starttime'] = UTCDateTime(df.corr_timestamp.iloc[0])
                sorted_df_list[i]['df'] = df
            else:
                break

    sorted_df_list_copy = deepcopy(sorted_df_list)
    # if the first record contains a clock flag, it's not very helpful 
    del_count = 0
    for i, df_list1 in enumerate(sorted_df_list_copy):
        if pd.isna(df_list1['segment_est_divergence']):
            logging.info('SEVERE: Removing first item from list because it contains clock flags {}'.format(os.path.basename(df_list1['gzip_filename'])))
            # sorted_df_list.pop(i)
            del sorted_df_list[i-del_count]
            del_count += 1
        else:
            break

    logging.info('#########')

    # calculate if any records might be suspect (use the Interquartile range)
    segment_est_divergence = []
    mean_mh1_all = []
    mean_mh2_all = []
    mean_mhz_all = []
    for x in sorted_df_list:
        logging.info(os.path.basename(x['gzip_filename']))
        if pd.notnull(x['segment_est_divergence']):
            segment_est_divergence.append(x['segment_est_divergence'])
        if pd.notnull(x['mh1_mean']):
            mean_mh1_all.append(x['mh1_mean'])
        if pd.notnull(x['mh2_mean']):
            mean_mh2_all.append(x['mh2_mean'])
        if pd.notnull(x['mhz_mean']):
            mean_mhz_all.append(x['mhz_mean'])

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

    mh1_suspect = 0
    mh2_suspect = 0
    mhz_suspect = 0

    mean_mh1_all = np.array(mean_mh1_all)
    sort_data = np.sort(mean_mh1_all)

    Q1 = np.percentile(sort_data, 25, interpolation = 'midpoint') 
    Q2 = np.percentile(sort_data, 50, interpolation = 'midpoint') 
    Q3 = np.percentile(sort_data, 75, interpolation = 'midpoint') 
      
    IQR = Q3 - Q1 

    low_lim = Q1 - 1.5 * IQR
    up_lim = Q3 + 1.5 * IQR

    for m in sort_data:
        if (m < low_lim) or (m > up_lim):
            mh1_suspect += 1
    if mh1_suspect > 0: 
        logging.info('SEVERE: Station - {} {} MH1 mean is potentially suspect.'.format(orig_station,mh1_suspect))
    

    mean_mh2_all = np.array(mean_mh2_all)
    sort_data = np.sort(mean_mh2_all)

    Q1 = np.percentile(sort_data, 25, interpolation = 'midpoint') 
    Q2 = np.percentile(sort_data, 50, interpolation = 'midpoint') 
    Q3 = np.percentile(sort_data, 75, interpolation = 'midpoint') 
      
    IQR = Q3 - Q1 

    low_lim = Q1 - 1.5 * IQR
    up_lim = Q3 + 1.5 * IQR

    for m in sort_data:
        if (m < low_lim) or (m > up_lim):
            mh2_suspect += 1
    if mh2_suspect > 0: 
        logging.info('SEVERE: Station - {} {} MH2 mean is potentially suspect.'.format(orig_station,mh2_suspect))

    mean_mhz_all = np.array(mean_mhz_all)
    sort_data = np.sort(mean_mhz_all)

    Q1 = np.percentile(sort_data, 25, interpolation = 'midpoint') 
    Q2 = np.percentile(sort_data, 50, interpolation = 'midpoint') 
    Q3 = np.percentile(sort_data, 75, interpolation = 'midpoint') 
      
    IQR = Q3 - Q1 

    low_lim = Q1 - 1.5 * IQR
    up_lim = Q3 + 1.5 * IQR

    for m in sort_data:
        if (m < low_lim) or (m > up_lim):
            mh2_suspect += 1
    if mhz_suspect > 0: 
        logging.info('SEVERE: Station - {} {} MHZ mean is potentially suspect.'.format(orig_station,mhz_suspect))

    config.rejected = 0
    starttime0 = None
    index0 = None
    for i, df_list1 in enumerate(sorted_df_list):
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
                orig_station = df_list1['orig_station']
                orig_ground_station = df_list1['orig_ground_station']
                logging.info('WARNING: Ground station/Station - {} {} {} frames are reset to zero (-7777 error)'.format(orig_ground_station,orig_station,gaps_7777))

            # import the stream for every record, using the same 
            # time_index0 and index0
            stream1 = stream_import(df,sample_time0,index0,attempt_merge)

            stream += stream1

    # if required, clean the spikes 
    if config.clean_spikes:
        for tr in stream:
            rec_spikes_found = despike3(tr)

    # if required, merge the streams 
    if config.combine_ground_stations: 
        # also includes clock correction
        merged_stream = merge_streams(stream)
        stream = merged_stream

    if config.combine_ground_stations:
        # trim_stream_midnight(stream,year=year,julday=julday)
        trim_stream_after_midnight(stream,year=year,julday=julday)
        if len(stream) > 0:
            trim_stream_before_midnight(stream,year=year,julday=julday,sorted_df_list=sorted_df_list)

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
        for tr in channel_stream:
            if INVALID in tr.data:
                tr.data = ma.masked_where(tr.data == INVALID, tr.data)
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
    
    if INVALID in tr_new.data:
        tr_new.data = ma.masked_where(tr_new.data == INVALID, tr_new.data)
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

        # merge the traces so that the timing is correct 

        # merge and create a masked trace 
        channel_stream.merge(method=1)
        # channel_stream1 = channel_stream.copy()
        # channel_stream1 = channel_stream1.split()
        # channel_stream1.write('/Users/cnunn/Downloads/merged.MSEED', format='MSEED')

        if len(channel_stream) > 1:
            logging.info('Too many values in channel stream')
            raise Exception

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

        # find where the clock flag is set 
        clk_set = np.where(tr_CLK.data==1)[0]
        if len(clk_set) > 0:
            date1 = tr_CLK.stats.starttime

            # get the ATT trace
            tr_ATT = merged_stream.select(channel='ATT')[0]

            # get the mask (data gaps) for the ATT trace
            mask_ATT = np.ma.getmask(tr_ATT.data)
        
            # # find where the software clock was set 
            # tr_CLK = merged_stream.select(channel='CLK')[0]
            # mask the ATT trace with the software clock 
            tr_CLK.data = np.ma.masked_where(tr_CLK.data==1,tr_CLK.data)                  # mask any values where the software clock has been set 
        
            # get the mask for the CLK trace
            mask_CLK = np.ma.getmask(tr_CLK.data)

            if tr_CLK.data[0] is ma.masked:
                logging.info('SEVERE: Initial record contains clock flag - {}.{}'.format(date1.year,date1.julday))   
                logging.info('temp stopping')
                exit()
                
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
                    logging.info('First date without a clock error: {}'.format(UTCDateTime(tr_ATT.data[first_good])))
                    logging.info('Last date without a clock error: {}'.format(UTCDateTime(tr_ATT.data[last_good])))
                    mask_CLK[last_good:] = False 
                    tr_CLK.data = ma.masked_array(tr_CLK.data, mask=mask_CLK)
                    
            
            logging.info('WARNING: Updating {} timestamps on the timing trace ATT due to a clock error (total={})'.format(len(np.where( mask_CLK == False )[0]), len(tr_ATT) ))

            tr_ATT.data= ma.masked_array(tr_ATT.data, mask=mask_CLK)

            # this does the work for interpolating the timing trace - it interpolates the time between good records
            data_series = pd.Series(tr_ATT.data)
            data_series.interpolate(method='linear', axis=0, limit=None, inplace=True, limit_direction=None, limit_area='inside', downcast=None)
            tr_ATT.data=data_series.to_numpy(dtype=np.float64)
        
            # apply the original mask to the trace (to remove the data gaps)
            tr_ATT.data = ma.masked_array(tr_ATT.data, mask=mask_ATT)
  
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
        correct_shifts(tr_ATT,tr_DL4)

        # no need to save the tr_DL4 trace 
        merged_stream.remove(tr_DL4)

    if config.combine_ground_stations:
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

    # make some checks
    
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

    df['last_valid_index'] = np.where((df['timestamp_orig'].isna()), pd.NA, df.index)

    # applying ffill() method to fill the missing values
    df['last_valid_index'].ffill(inplace=True)

    df['last_valid_divergence'] = np.where((df['timestamp_orig'].isna()), pd.NA, df.divergence)
    # applying ffill() method to fill the missing values
    df['last_valid_divergence'].ffill(inplace=True)


    df['last_valid_delta4'] = np.where((df['timestamp_orig'].isna()), pd.NA, df.delta4)
    # applying ffill() method to fill the missing values
    df['last_valid_delta4'].ffill(inplace=True)

    # estimate the divergence (for when there is no valid timestamp)
    df['estimated_divergence'] = df['last_valid_divergence'] + ((df['last_valid_index'] - df.index) * (df.last_valid_delta4 - DELTA*4))

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
            logging.info('WARNING: corrected tapehead jump from original start timestamp {} difference {:.01f} s, possible error {:.01f} s [auto]'.format(UTCDateTime(orig_timestamp),jump,possible_correction))


    # applying ffill() method to fill the missing values
    df.possible_correction.ffill(inplace=True)

    df['divergence_diff'] = df['divergence_diff'].round(1)
    df['correction'] = np.where((df.possible_correction == df.divergence_diff), df['possible_correction'], 0.0)

    correction_sum =  ((df['correction'].lt(0)) | (df['correction'].gt(0))).sum()
    if correction_sum > 0: 

        logging.info('WARNING: {} timestamps corrected due to 400/800 ms tape head error'.format(correction_sum))

        # reuse the list we've already created
        for idx in idx_list:
            orig_time = UTCDateTime(df.timestamp_orig.iloc[idx])
            corr = df['correction'].iloc[idx]
            logging.info('CORRECTION,{},{},{:03},Constant [auto],{:.01f} s,{},,,'.format(station, orig_time.year,orig_time.julday, corr, orig_time))


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

    logging.info('INFO: Zero index record. \nTrace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence={:.02f} Est. Divergence={:.02f} Reject={} Mean MH1={:.01f} Mean MH2={:.01f} Mean MHZ={:.01f} File={}'.format(index_no, orig_station, corr_ground_station, starttime0, sample_time0, starttime0-sample_time0, segment_est_divergence,df_list1['segment_est_divergence_incl_clock'],df_list1['reject'],df_list1['mh1_mean'],df_list1['mh2_mean'],df_list1['mhz_mean'], gzip_filename))

    calculate_last_valid(df,index_no,sorted_df_list)

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
    segment_delta4 = df_list1['segment_delta4']
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
    start_delta4 = df_list1['segment_delta4']


    # find the previous dataframe
    df_list_prev = sorted_df_list[index_no-1]
    # prev_starttime = df_list_prev[0]
    # prev_endtime = df_list_prev[1]
    
    prev_df = df_list_prev['df']
    end_delta4 = df_list_prev['segment_delta4']

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

    if pd.isna(start_delta4):
        start_delta4 = clock_flag = df.delta4[0]        
    if pd.isna(end_delta4):
        end_delta4 = config.valid_end_delta4 


    timestamp_time_diff = (start_timestamp - end_timestamp).total_seconds()
    
    starttime= UTCDateTime(start_timestamp)

    found_index = []
    method=''
    
    # logging.info('{} Original starttime (data time)={} 1st frame={} Original endtime={} Last frame={} ground_station={}'.format(index_no, starttime, frame, endtime, df.frame.iloc[-1],corr_ground_station))

    overlap_found = False
    include_trace = False

    if start_timestamp <= (end_timestamp + pd.Timedelta(seconds=25)):
        start_timestamp_match = start_timestamp - pd.Timedelta(seconds=25)
        end_timestamp_match = start_timestamp + pd.Timedelta(seconds=25)
        overlap_found = True

        found_index_all = prev_df[prev_df['frame'] == frame].index

        # look for a match on the frame
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
            # note that time index is every 4 records (because there are 
            # four records per frame number)
            df.time_index = np.arange(len(df)) + prev_time_index
            sample_time = sample_time0 + df.iloc[0].time_index *DELTA*4

            logging.info('INFO: Overlapping traces found. \nTrace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence={:.02f} Est. Divergence={:.02f} Reject={} Mean MH1={:.01f} Mean MH2={:.01f} Mean MHZ={:.01f} File={}'.format(index_no, orig_station, orig_ground_station, UTCDateTime(starttime), sample_time, UTCDateTime(starttime)-sample_time,segment_est_divergence,df_list1['segment_est_divergence_incl_clock'],df_list1['reject'],df_list1['mh1_mean'],df_list1['mh2_mean'],df_list1['mhz_mean'],gzip_filename))

        elif len(found_index) > 1:
            print('EXCEPTION: Too many found.')
            logging.info('EXCEPTION: Too many found.')
            raise Exception


    # if the item has not been found yet, try a different way 
    if len(found_index) == 0:

        gradient = (start_delta4 + end_delta4)/2

        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(last_frame=end_frame,new_frame=frame,gap=timestamp_time_diff,delta4_1=end_delta4,delta4_2=start_delta4)
        # if frame_gap_correct:
        index_diff = actual_frame_gap

        # find the time interval to get the sample time 
        intervals = end_time_index + index_diff
        df.time_index = np.arange(len(df)) + intervals
        sample_time = sample_time0 + df.iloc[0].time_index *DELTA*4

        sample_time_end = sample_time0 + end_time_index *DELTA*4

        if absolute_error > 2:
            logging.info('WARNING: Possible error with the time gap. \nTrace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence={:.02f} Est. Divergence={:.02f} Reject={} Mean MH1={:.01f} Mean MH2={:.01f} Mean MHZ={:.01f} File={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time,segment_est_divergence,df_list1['segment_est_divergence_incl_clock'],df_list1['reject'],df_list1['mh1_mean'],df_list1['mh2_mean'],df_list1['mhz_mean'],gzip_filename))
        else:
            logging.info('INFO: Time gap with correct number of missing samples found. \nTrace index={} Station={} Ground Station={} Original Timestamp={} Sample Time={} Adjustment Time={} Segment Divergence={:.02f} Est. Divergence={:.02f} Reject={} Mean MH1={:.01f} Mean MH2={:.01f} Mean MHZ={:.01f} File={}'.format(index_no, orig_station, orig_ground_station, starttime, sample_time, starttime-sample_time,segment_est_divergence,df_list1['segment_est_divergence_incl_clock'],df_list1['reject'],df_list1['mh1_mean'],df_list1['mh2_mean'],df_list1['mhz_mean'],gzip_filename))

    # TODO how do we sort later - should the timing be updated? 
    df_list1['attempt_merge'] = attempt_merge

    calculate_last_valid(df,index_no,sorted_df_list)

def calculate_divergence(delta4):

    est_divergence =  ((EXACT_N - 1) * DELTA * 4) - ((EXACT_N - 1) * delta4)
    return est_divergence
    
def calculate_segment_mean(df):

    df1 = df[['orig_mh1_1','orig_mh1_2','orig_mh1_3','orig_mh1_4']].stack().reset_index()
    mh1_mean = df1[0].mean(skipna=True)
    df2 = df[['orig_mh2_1','orig_mh2_2','orig_mh2_3','orig_mh2_4']].stack().reset_index()
    mh2_mean = df2[0].mean(skipna=True)
    dfz = df[['orig_mhz_1','orig_mhz_2','orig_mhz_3','orig_mhz_4']].stack().reset_index()
    mhz_mean = dfz[0].mean(skipna=True)
    
    return mh1_mean, mh2_mean, mhz_mean

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

    if pd.isnull(segment_est_divergence):
        start_timestamp = UTCDateTime(df.corr_timestamp.iloc[0])   
        end_timestamp = UTCDateTime(df.corr_timestamp.iloc[-1])
        time_diff = end_timestamp - start_timestamp
        segment_delta4 = time_diff / (len(df)-1)

        segment_est_divergence_incl_clock1 =  (len(df)-1) * DELTA *4 - time_diff
        # now scale for 24 hours 
        segment_est_divergence_incl_clock = segment_est_divergence_incl_clock1 * 86400 / time_diff
    else: 
        segment_est_divergence_incl_clock = segment_est_divergence

    return segment_delta4, segment_est_divergence, segment_est_divergence_incl_clock

def calculate_last_valid(df,index_no,sorted_df_list):

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
            config.valid_end_delta4 = sorted_df_list[index_no]['segment_delta4']

            config.clock_end_frame = pd.NA
            config.clock_end_timestamp = pd.NA
            config.clock_end_time_index = pd.NA
            config.clock_end_delta4 = pd.NA

        else: 

            config.clock_end_frame = df.frame.iloc[-1]
            config.clock_end_timestamp = df.corr_timestamp.iloc[-1]
            config.clock_end_time_index = df.time_index.iloc[-1]
            config.clock_end_delta4 = sorted_df_list[index_no]['segment_delta4']


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
    
    current_filename = '%s.%s.%s.%s.%s.%s.%s.%s' % (network,station, location, channel,
      year, julday, rev, ext)

    return current_filename


def trim_stream_after_midnight(stream,year,julday):
    # trim stream to the end of the day if necessary 
    stream_ATT_orig = stream.select(channel='ATT')
    start_time = (UTCDateTime(year=year,julday=julday)).timestamp
    end_time = (UTCDateTime(year=year,julday=julday) +3600*24).timestamp

    if len(stream_ATT_orig) > 1:
        print('Exception - too many traces in stream.')
        raise Exception
    
    if stream_ATT_orig[0].data[-1] < end_time:
        return   

    index_end = np.where(stream_ATT_orig[0].data < end_time)[0]
    if len(index_end) > 0: 
        index_end = index_end[-1]
        trim_end_time = stream_ATT_orig[0].times(type='utcdatetime')[index_end]
    else:
        trim_end_time = UTCDateTime(end_time)

    # add a small margin to deal with issue that all the other traces come
    # after the ATT trace
    trim_end_time = trim_end_time + DELTA*4 - 1/106
    stream = stream.trim(endtime=trim_end_time)

    if len(stream) > 0:
        logging.info('Trimming: Trim_end_time {} Stream End Time {}'.format(trim_end_time, stream.select(channel='ATT')[0].stats.endtime))
    else:
        logging.info('Removed stream: Trim_end_time {}'.format(trim_end_time))





def trim_stream_before_midnight(stream,year,julday,sorted_df_list):
    stream_ATT_orig = stream.select(channel='ATT')
    start_time = UTCDateTime(year=year,julday=julday).timestamp
    end_time = (UTCDateTime(year=year,julday=julday) +3600*24).timestamp


    if len(stream_ATT_orig) > 1:
        print('Exception - too many traces in stream.')
        raise Exception

    trim_start_time = None

    if stream_ATT_orig[0].data[0] > start_time:
        return   

    index_start= np.where(stream_ATT_orig[0].data > start_time)[0]

    if len(index_start) > 0: 
        index_start = index_start[0]
        # print(UTCDateTime(stream_ATT_orig[0].data[index_start]))
        trim_start_time = stream_ATT_orig[0].times(type='utcdatetime')[index_start]
    else:
        trim_start_time = UTCDateTime(start_time)

    

    stream = stream.trim(starttime=trim_start_time)

    if len(stream) > 0:
        logging.info('Trimming: Trim_start_time {} Stream Start Time {}'.format(trim_start_time, stream.select(channel='ATT')[0].stats.starttime))
        
        starttime0 = UTCDateTime(stream_ATT_orig[0].data[0])
        # first check that the starttime0 is near midnight
        midnight = UTCDateTime(year=year,julday=julday)
        test_time = midnight + 0.6091
        if starttime0 < test_time:
            sample_time0 = starttime0

        else:
            logging.info('WARNING: Recalculated starttime not immediately after midnight')

            segment_delta4 = DELTA*4
            for x in sorted_df_list:
                if pd.notnull(x['segment_delta4']):
                    segment_delta4 = x['segment_delta4']
                    break

            time_int = round((starttime0 - midnight)/segment_delta4)-1
            
            est_start_time = starttime0 - time_int *segment_delta4
            
            sample_time0 = est_start_time + time_int *DELTA*4

        if len(stream) > 5:
            print('There should not be more than five traces in the stream.')
            raise Exception
        for tr in stream:
            tr.stats.starttime = sample_time0
                
            


    else:
        logging.info('Removed stream: Trim_start_time {}'.format(trim_start_time))



def save_stream(stream,join_dir):
    logging.debug('save_stream()')
    logging.debug(stream)
    # for station in ('S12', 'S14', 'S15', 'S16'):
    #     for location in ('01','02','03','04','05','06','07','08','09',
    #           '10','11','12','13','14'):
    #         # location_stream = stream.select(location=location)
    #         location_stream = stream.select(station=station, location=location)

    for station in ['S11','S12','S14','S15','S16']:
        station_stream = stream.select(station=station)

        if len(station_stream) > 0:
            locations = set()
            if config.combine_ground_stations:
                for tr in station_stream:
                    locations.add(tr.stats.location)
            else: 
                for tr in station_stream:
                    locations.add(tr.stats.ground_station)
                    tr.stats.location = tr.stats.ground_station
            logging.debug(station_stream.select(channel='ATT'))
            for location in locations:
                location_stream = station_stream.select(location=location)
                # need directory with separate station, year and day
                # make the output directory
                station = location_stream[0].stats.station
                starttime = location_stream[0].stats.starttime
                directory = make_dir(join_dir,station,starttime,lower=True)

                for channel in ('AFR', 'ATT', 'MH1', 'MH2', 'MHZ', 'SHZ',):
                    # note that 'AFR only exists if 'config.combine_ground_stations == False
                    channel_stream = location_stream.select(channel=channel)   
                    if station == 'S14' and channel == 'MHZ' and starttime > UTCDateTime('1972-03-12T00:00') and starttime < UTCDateTime('1976-11-17T00:00:00.000000Z'):
                        # we are not saving the MHZ track during this time frame. 
                        continue

                    if station == 'S12' and channel == 'SHZ':
                        # we are not saving the SHZ track for S12
                        continue


                    if len(channel_stream) > 0:
                        filename = filename_from_trace(channel_stream[0],directory,lower=False)
                        filename = os.path.join(directory,filename)

                        if channel == 'ATT':
                            for trace in channel_stream:
                                if np.ma.isMaskedArray(trace.data):
                                    # mask the ATT trace which has gaps with a fill value
                                    trace.data = trace.data.filled(fill_value=-1.0)

                        if channel == 'AFR':
                            for trace in channel_stream:
                                if np.ma.isMaskedArray(trace.data):
                                    # mask the ATT trace which has gaps with a fill value
                                    trace.data = trace.data.filled(fill_value=-1.0)

                        if channel in ('MH1','MH2','MHZ','SHZ'):
                            for trace in channel_stream:
                                if np.ma.isMaskedArray(trace.data):
                                    # mask the ATT trace which has gaps with a fill value
                                    trace.data = trace.data.filled(fill_value=INVALID)
                
                        logging.info('Writing file {}'.format(filename))

                        channel_stream.write(filename, format='MSEED')

            print('Wrote files to ', directory)

def loose_frame_diff(last_frame,new_frame,gap,delta4_1=None,delta4_2=None):
    if delta4_1 is not None and delta4_2 is not None:
        delta4 = (delta4_1 + delta4_2)/2
    else:
        if delta4_1 is not None:
            delta4 = delta4_1
        elif delta4_2 is not None:
            delta4 = delta4_2
        else:
            delta4 = DELTA*4

    logging.debug('delta4_1 {} {} delta4_2 {} {}'.format(delta4_1, calculate_divergence(delta4_1), delta4_2, calculate_divergence(delta4_2) ))

    # tolerance (plus or minus) in seconds
    TOLERANCE = 27.171

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
            # full_frames_est = int(est_frame_gap // 90)
            full_frames_est = est_frame_gap // 90
            full_frames_est = int(full_frames_est)



            # full_frames_est can be wrong if we are very close to 
            # dividing by 90 frames, so we check the esimate minus and 
            # plus one, and then find the value of the single gap closest
            # to the nominal inteval value 
            # actual_frame_gap is the total number of frames in the gap (that 
            # fits with the last_frame and new_frame)
            actual_frame_gap_plus = (full_frames_est+1)*90 + diff
            # single_gap_plus = gap/actual_frame_gap_plus
            # est_divergence_plus = calculate_divergence(single_gap_plus)
            est_gap_plus = actual_frame_gap_plus * delta4
            absolute_error_plus = gap - est_gap_plus
            # logging.debug('est_divergence_plus {}'.format(est_divergence_plus))
            # logging.debug('absolute_error_plus {}'.format(absolute_error_plus))

            actual_frame_gap_normal = (full_frames_est)*90 + diff
            # single_gap_normal = gap/actual_frame_gap_normal
            # est_divergence_normal = calculate_divergence(single_gap_normal)
            est_gap_normal = actual_frame_gap_normal * delta4
            absolute_error_normal = gap - est_gap_normal
            # logging.debug('est_divergence_normal {}'.format(est_divergence_normal))
            # logging.debug('absolute_error_normal {}'.format(absolute_error_normal))

            actual_frame_gap_minus = (full_frames_est-1)*90 + diff
            # single_gap_minus = gap/actual_frame_gap_minus
            # est_divergence_minus = calculate_divergence(single_gap_minus)
            est_gap_minus = actual_frame_gap_minus * delta4
            absolute_error_minus = gap - est_gap_minus
            # logging.debug('est_divergence_minus {}'.format(est_divergence_minus))
            # logging.debug('absolute_error_minus {}'.format(absolute_error_minus))

            # find the best solution using the absolute error 
            if abs(absolute_error_normal) < abs(absolute_error_plus):
                absolute_error = absolute_error_normal
                actual_frame_gap = actual_frame_gap_normal
                # est_divergence = est_divergence_normal
            else: 
                absolute_error = absolute_error_plus
                actual_frame_gap = actual_frame_gap_plus
                # est_divergence = est_divergence_plus

            if abs(absolute_error_minus) < abs(absolute_error):
                absolute_error = absolute_error_minus
                actual_frame_gap = actual_frame_gap_minus
                # est_divergence = est_divergence_minus



        # if the gap was negative 
        else:
            full_frames_est = 0
            # single_gap = gap/diff
            actual_frame_gap = diff

            # actual_frame_gap_normal = (full_frames_est)*90 + diff
            single_gap = gap/actual_frame_gap
            # est_divergence = calculate_divergence(single_gap)


        est_gap = actual_frame_gap * delta4

        logging.info('actual_frame_gap {}, gap={}, est_gap={}, diff {},full_frames_est {}'.format(actual_frame_gap,gap,est_gap,diff,full_frames_est))


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

    # check for any nulls at the beginning and remove    
    if pd.isna(df['orig_no'].iloc[0]) or pd.isna(df['orig_no'].iloc[1]):
        first_good_idx = df['orig_no'][1:].first_valid_index()
        if first_good_idx is None:
            df = df.iloc[0:0]
            return df
        else:
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

    # start by masking the INVALID data
    if INVALID in trace.data:
        trace.data = ma.masked_where(trace.data == INVALID, trace.data)

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

    # put the INVALID values back again 
    if np.ma.isMaskedArray(trace.data):
        trace.data = trace.data.filled(fill_value=INVALID)



                # diff < 2 and spike greater than TEST_SPIKE_V_SMALL
                # diff < 5 and spike greater than TEST_SPIKE_SMALL
                # diff < 10 and spike greater than TEST_SPIKE_HIGH



    return rec_spikes_found


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

    mh2 = df[['orig_mh2_1', 'orig_mh2_2','orig_mh2_3','orig_mh2_4']].to_numpy(dtype='int32',na_value=INVALID)
    mh2 = mh2.flatten()

    mhz = df[['orig_mhz_1', 'orig_mhz_2','orig_mhz_3','orig_mhz_4']].to_numpy(dtype='int32',na_value=INVALID)
    mhz = mhz.flatten()

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
    frame = frame.flatten()

    # get the absolute times in seconds since 1 Jan 1970
    df['abs_times'] = df['corr_timestamp'].values.astype(np.int64) / 10 ** 9

    abs_times = df[['abs_times']].to_numpy(dtype='float64',na_value=None)
    abs_times = abs_times.flatten()

    # make a mask for corr_timestamp 
    mask = df[['orig_no']].to_numpy(dtype='int32',na_value=INVALID)
    mask = mask.flatten()

    mask = ma.masked_equal(mask, INVALID)
    mask = ma.getmask(mask)

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