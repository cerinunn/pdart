#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Check the traces for duplicates (when one station is copied to another)

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
# slightly shorter delta, so that the timing doesn't run over 24 hours 
# SHORT_DELTA = 0.1509000000
INTERVAL = '603774us'

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



def call_check_duplicates(
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
  
    ):

    '''
    Make files which are joined into the right days
    Calls csv_import_work_tapes()
    '''


    # TODO -fix these date ranges - they are not quite right 
    # if no filenames have been passed, then look for them 
    for year in range(year_start,year_end+1):
        log_filename = 'duplicate.{}.log'.format(year)
        log_filename = os.path.join(log_dir,log_filename)

        logging.basicConfig(filename=log_filename, filemode='w', 
          level=logging_level,format='%(message)s')

        logging.info('############################################')
        logging.info('Processing {}'.format(year))
        logging.info('############################################')



        for day in range(day_start,day_end+1):
            logging.info('#########')  
            print('Processing {}.{}'.format(year,day))
            logging.info('Processing {}.{}'.format(year,day))

            df_list = []
                                

            wildcard_filename = '{}.*.*.*.*.{}.{:03}.*.csv.gz'.format(wildcard_style,year,day)
            print('wildcard filename ', processed_dir, wildcard_filename)
                    # read in each of the csv files for the station
                    # for filename in glob.glob('/Users/cnunn/python_packages/pdart/examples/test.csv'):
            for filename in glob.glob(os.path.join(processed_dir,wildcard_filename)):

                # if os.path.basename(filename) not in ('wtn.1.3.S12.1.1976.061.13_29_58_945000.csv.gz'):
                #     continue

                if 'dropped' not in filename:
                    try: 
                        gzip_filename = filename

                        # logging.info(gzip_filename)
                        df_list = read_file(filename,df_list,logging_level,join_dir,log_dir)
                    except Exception as e:
                        logging.info(traceback.format_exc())
                        logging.info(traceback.print_exc(file=sys.stdout))
                        # reset_config()
                        print('Warning, continuing')
                        # print('Not continuing')

                # logging.info('#########')

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

    # create a hash on the combined columns
    df['hash'] = pd.util.hash_pandas_object(df[[
        'orig_mh1_1','orig_mh2_1','orig_mhz_1',
        'orig_mh1_2','orig_mh2_2','orig_mhz_2',
        'orig_mh1_3','orig_mh2_3','orig_mhz_3',
        'orig_mh1_4','orig_mh2_4','orig_mhz_4',
        'frame']], 
        index=False)

    # ignore the frame, to make it easier to check for nulls
    # create a hash on the combined columns
    df['hash2'] = pd.util.hash_pandas_object(df[[
        'orig_mh1_1','orig_mh2_1','orig_mhz_1',
        'orig_mh1_2','orig_mh2_2','orig_mhz_2',
        'orig_mh1_3','orig_mh2_3','orig_mhz_3',
        'orig_mh1_4','orig_mh2_4','orig_mhz_4',]], 
        index=False)

        

    # df.set_index('hash', inplace=True)


    if len(df) < 1:
        logging.debug('WARNING: No valid record(s), not importing'.format(len(df)))
        return df_list

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
      'gzip_filename' : gzip_filename,  'gzip_duplicate' : None, 'duplicate_count' : 0 }
    # make a list of the dataframes
    df_list.append(df_dict)

    return df_list


def process_list(df_list,join_dir,year,julday):

    # for each station, find the first anchor
    # run through the rest of the dataframes for the station

    sorted_df_list =  (sorted(df_list, key = lambda i: i['orig_station']))
    # print(type(sorted_df_list))
    
    file_list = []

    for x in sorted_df_list:
        file_list.append(x['gzip_filename'])

    from itertools import combinations

    for combo in combinations(file_list,2):
        combo1 = file_list.index(combo[0])
        combo2 = file_list.index(combo[1])
        
        station1 = sorted_df_list[combo1]['orig_station']
        station2 = sorted_df_list[combo2]['orig_station']

        # no need to check for duplicates for the same station
        if station1 == station2:
            continue

        # print(df_list[combo1])
        df1 = sorted_df_list[combo1]['df']
        # logging.info(df1.index.head().to_string())
        df2 = sorted_df_list[combo2]['df']
        # logging.info(df2['hash'].head().to_string())

        df1['found'] = df1.hash.isin(df2.hash)
        combo_sum = df1['found'].sum()
        if combo_sum > 0:
            df1['found'] = np.where((df1['hash2']== 14531169447873018637), False,  df1['found'])
            combo_sum = df1['found'].sum()
            percentage = 100 * combo_sum/len(df1)
            if percentage > 5:
                logging.info('Duplicate found: {} {} duplicates={} left={} right={}'.format(combo[0],combo[1],combo_sum,len(df1),len(df2)))
                if logging.DEBUG >= logging.root.level:
                    df_display = df1[(df1['found']==True)]            
                    logging.debug(df_display.head().to_string())
        # else: 
        #     logging.info('Not found    : {} {} 0'.format(combo[0],combo[1]))
        



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

    # logging.info(df.head().to_string())
    # logging.info(df.tail().to_string())

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