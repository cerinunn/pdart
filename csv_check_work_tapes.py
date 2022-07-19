#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run initial checks on csv files

Make checked csv files (will drop some damaged data)
out of the raw csv files extracted from binary. 

This uses the following workflow:

# split the dataframe into station and ground station
df_st_list = all_split(df,gzip_filename)

config.cumsum_test = cumsum_test
# This is a key bit of code to test
# check the timeseries data and remove any that don't fit
df_gst, df_dropped, df_orig, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps, rec_repeated_frame, compliance_passed = check_compliance(df_gst,gzip_filename)

    # Fix the frame number if it was incorrectly recorded 
    df_gst, rec_fixed_frame, rec_repeated_frame = station_fix_frames(df_gst)

    # Where possible, do a really simple fix on the timestamp
    df_gst, rec_fixed_simple_timestamp = station_simple_fix_timestamp(df_gst)

    # find the cumulative sum of the good records 
    df_gst = calculate_cumsum(df_gst)

    # estimate the delta as it changes over time
    df_gst = calculate_delta4(df_gst)

    # fix timestamps if possible (called from check_compliance)
    df_gst, df_dropped, df_orig, rec_adjusted_timestamps, compliance_passed = fix_timestamps(df_gst, df_dropped)

        # fix any gaps (called from fix_timestamps, which is called from check_compliance)
        Search for gaps in the record. Blank records are inserted if the 
        following conditions are met:
        Can either be called by the bad record, or for the whole 
        dataframe.
        df_gst = station_fix_missing_timestamps2(df_gst)

    split_good_records
    finally, save the good files and the dropped files (note that the dropped files 
    are concatenated - so bear the timestamp of the first record)
    if there are a lot of missed records, then reducing config.cumsum_final_test
    to a smaller number (e.g. 10), may reduced the missed records
    (but make require some manual corrections at the join stage)

    See also the notes in config.py - there are some adjustable parameters.

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
import math
import numpy as np
import numpy.ma as ma
import logging, logging.handlers
import csv
import fnmatch
import shutil
import pandas as pd
from random import choice, seed
from itertools import groupby
import sys, traceback
from collections import OrderedDict

from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, Trace, Stats, read

# from pdart.save_24_hours import update_starttimes
import pdart.config as config
# from pdart.csv_import_work_tapes import find_output_dir, make_output_dir, make_filelist
import matplotlib.pyplot as plt





# global DELTA
DELTA = 0.1509433962
INTERVAL = '603774us'

INVALID = -99999
# X should point north, Y east, but this is not always the case
# so we rename LPX to MH1, and LPY to MH2
NETWORK='XA'


# global first_record
# global last_record

# consecutive frames
# FRAMES = list(range(0,90))

def call_csv_check_work_tapes(
    checked_dir='.',
    processed_dir='.',
    log_dir='.',
    filenames=None,
    # logging_level=logging.DEBUG
    logging_level=logging.INFO,
    single_station=None,
    single_ground_station=None
    ):

    '''
    Makes initial stream files from gzipped csv files which have been checked.
    Calls csv_import_work_tapes()
    '''

    # if no filenames have been passed, then look for them 
    if filenames is None:
        filenames = []
        for filename in glob.glob(os.path.join(checked_dir,'wtn*.csv.gz')):
            filenames.append(os.path.basename(filename))

    for filename in filenames: 
        print('Searching for :', os.path.join(checked_dir,filename))
        try:
            for in_file in glob.glob(os.path.join(checked_dir,filename)):    
                orig_timestamps = csv_check_work_tapes(in_file,single_station,
                  single_ground_station,logging_level,processed_dir,log_dir)
        except Exception as e:
            logging.info(traceback.format_exc())
            logging.info(traceback.print_exc(file=sys.stdout))
            logging.info('FINAL: Statistics - Exited before end')
            print('Warning, continuing')

        # close the log file so that we can open a new one
        logger = logging.getLogger()
        handlers = logger.handlers[:]
        for handler in handlers:
            handler.close()
            logger.removeHandler(handler)
            

def csv_check_work_tapes(gzip_filename,single_station=None,
  single_ground_station=None,logging_level=logging.INFO,processed_dir='.',log_dir='.',delete_previous=True):
    """
    Makes a stream file from the gzipped csv file. The timing is based on
    the timing in the csv file, and so the seismograms are not yet
    continuous.

    csv_import() makes streams from the individual csv file.

    :type gzip_filename: str
    :param gzip_filename: gzipped filename of the CSV file to be read.
    :type single_station: str
    :param single_station: The name of a single station can be specified (the 
    entire file is read in first)
    :type single_ground_station: int
    :param single_station: The name of a single ground station can be specified 
    (the entire file is read in first)
    :rtype: :class:`~obspy.core.stream.Stream`
    :return: A ObsPy Stream object.
    """

    print('Processing file: ', gzip_filename)
    print('Writing to ', processed_dir)
    print('Log dir ', log_dir)
    

    log_filename = os.path.basename(gzip_filename)
    log_filename = log_filename + '.log'
    log_filename = os.path.join(log_dir,log_filename)
    print('log filename ', log_filename)
    
    logging.basicConfig(filename=log_filename, filemode='w', 
      level=logging_level,format='%(message)s')

    # reset the global variables
    global first_record
    first_record = None
    global last_record 
    last_record = None
    global rec_initial_total
    rec_initial_total = None

    global gst_total
    gst_total = 0
    global total
    total = 0

    logging.info('############################################')
    logging.info(gzip_filename)
    logging.info('############################################')
    print(gzip_filename)

    logging.info('config.low_tolerance={}'.format(config.low_tolerance))
    logging.info('config.high_tolerance={}'.format(config.high_tolerance))
    logging.info('config.lower_tolerance={}'.format(config.lower_tolerance))
    logging.info('config.higher_tolerance={}'.format(config.higher_tolerance))
    
    df = None
    df = pd.read_csv(gzip_filename, dtype=str)


    if len(df) == 0:
        delete_previous = False

    if delete_previous:
        out_gzip_filename = os.path.basename(gzip_filename)
        out_gzip_filename = os.path.join(processed_dir,out_gzip_filename)
        delete_wildcard = out_gzip_filename.replace('.csv.gz', '.*')

        delete_list = glob.glob(delete_wildcard, recursive=False)
        for file1 in delete_list:
            if '.log' not in file1:
                os.remove(file1)
                print('Removing file ', file1)



    default_config()

    df = initial_cleanup(df)

    df = drop_last_frame(df)

    initial_report(df)

    all_station_order(df,gzip_filename)

    df, rec_station_duplicates = all_station_duplicates_v2(df,gzip_filename)
    
    all_flat_seismogram(df,gzip_filename)

    if config.initial == False:
        df, rec_drop_damaged_total = all_known_errors(df,gzip_filename)
    else:
        rec_drop_damaged_total = 0
    
    df, rec_damaged_sync_total = all_delete_damaged(df)

    # station
    # df = all_station_check(df,station_order)

    # duplicates are not currently removed (checking at the end instead)
    rec_duplicates = 0
    # df, rec_duplicates = all_drop_duplicates(df)

    detailed_report(df)

    # check for non-consecutive stations  
    # all_check_stations(df)

    # split into station and ground station

    logging.debug('df total {}'.format(len(df)))

    # split the dataframe into station and ground station
    df_st_list = all_split(df,gzip_filename)

    # rec_missing_total = 0
    # rec_final_total = 0 
    # rec_negative_total = 0
    # rec_fixed_frame_total = 0
    # rec_adjusted_timestamps_total = 0 
    # rec_deleted_timestamps_total = 0
    # rec_repeated_frame_total = 0
    # rec_fixed_simple_timestamp_total = 0
    # # timestamps may be adjusted more than once - this gives the final number
    # rec_adjusted_timestamps_final_total = 0 

    rec_final_total = 0 

    rec_negative_total = 0
    rec_clock_flag_total = 0
    rec_adjusted_timestamps_final_total = 0
    rec_deleted_timestamps_total = 0
    rec_missing_total = 0

    rec_fixed_frame_total = 0
    rec_fixed_simple_timestamp_total = 0
    rec_adjusted_timestamps_total = 0
    rec_deleted_non_uniq_total = 0
    rec_repeated_frame_total = 0

    rec_dropped_total = 0 

    overall_start_times = []


    # go through each station and ground station combination
    for dict1 in df_st_list:
        
        start_times = []
        end_times = []
        lengths = []
        final_station = []
        final_ground_station = []

        df_gst = dict1['df_gst']
        corr_ground_station = df_gst['corr_ground_station'].iloc[0]
        orig_station = df_gst['orig_station'].iloc[0]
        orig_timestamp = df_gst['orig_timestamp'].iloc[0]   

        # logging.info('Temp 8- Length of this section: {}'.format(len(df_gst)))

        # the station and timestamp of the last records for the whole 
        # ground station
        config.last_station = dict1['last_station']
        config.last_timestamp = dict1['last_timestamp']
        config.last_orig_idx = df_gst['orig_idx'].iloc[-1]

        # If required, process a single station or ground station (or both)
        if single_station is not None:
            if orig_station != single_station:
                continue
        if single_ground_station is not None:
            if corr_ground_station != single_ground_station:
                continue

        # if orig_timestamp != pd.Timestamp('1977-04-20 02:28:58.870000+00:00'):
        #     print(orig_timestamp)
        #     continue

        logging.info('############ Ground station: {} Station: {} {}'.format(corr_ground_station,orig_station,orig_timestamp))

        if len(df_gst) == 0:
            continue

        global df_orig_shown
        df_orig_shown = False

        if config.initial:
            # use this setting to view large breaks and gaps 
            station_test_missing_timestamps(df_gst,gzip_filename)
            return df_gst, None, 0, 0
        else: 
            station_test_missing_timestamps(df_gst,gzip_filename)

        df_original = df_gst.copy()

        cumsum_test_list = [10, 20, 90, 180]

        # loop through, increasing the number of records in the cumulative sum test 
        # [the idea is that a big section of records that work together are probably right]
        for cumsum_test in cumsum_test_list:

            logging.debug('config.cumsum_test={}'.format(config.cumsum_test))

            config.cumsum_test = cumsum_test
            # each time, start from the original data
            df_gst = df_original.copy()
            
            # check the timeseries data and remove any that don't fit
            # note that both df_gst and df_dropped are taken from the last one (that hopefully works)
            df_gst, df_dropped, df_orig, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps, rec_repeated_frame, compliance_passed = check_compliance(df_gst,gzip_filename)
            if cumsum_test > 10:
                logging.info('WARNING: Ground station/Station - {} {} Using cumsum_test={}'.format(
                  corr_ground_station, orig_station,cumsum_test))

            if cumsum_test == 180:
                logging.debug('------------------')
                logging.debug(df_orig.tail().to_string())
                logging.debug(df_dropped.tail().to_string())
                logging.debug(df_gst.tail().to_string())
                logging.debug('------------------')
            
            # # there's no reason to carry on if all records have been deleted 
            # if len(df_gst) == 0:
            #     break

            # break out if bad records have been deleted (also inlcudes only bad reocrds)
            if compliance_passed:
                # logging.info('Compliance passed at cumsum_test={}'.format(cumsum_test))
                # logging.info(len(df_dropped))
                break
            # else:
            #     logging.info('Compliance not passed at cumsum_test={}'.format(cumsum_test))
            #     logging.info(len(df_dropped))

                # logging.debug('--------------------------')
                # logging.debug(df_gst.to_string())
                # df_orig_shown = True
                # logging.debug('--------------------------')


        rec_fixed_frame_total += rec_fixed_frame
        rec_fixed_simple_timestamp_total += rec_fixed_simple_timestamp
        rec_adjusted_timestamps_total += rec_adjusted_timestamps
        rec_deleted_non_uniq_total += rec_deleted_non_uniq
        rec_deleted_timestamps_total += rec_deleted_timestamps
        rec_repeated_frame_total += rec_repeated_frame


        
        if config.initial:
            continue

        if len(df_gst) > 0:

            # find the cumulative sum of the good records 
            # this looks at any that have already been fixed 
            df_gst = calculate_cumsum2(df_gst)

            logging.info('still there?')
            logging.debug(df_gst.to_string())
            logging.info('still there?')
            

            df_list1, df_list_bad = split_good_records(df_gst)

            

            for df1 in df_list1:

                # make some final checks to make sure it works 
                df1 = station_final_check(df1, df_orig)

                start_time = UTCDateTime(df1['corr_timestamp'].iloc[0])
                start_times.append(start_time)
                end_time = UTCDateTime(df1['corr_timestamp'].iloc[-1])
                end_times.append(end_time)
                lengths.append(len(df1))
                final_station.append(orig_station)
                final_ground_station.append(corr_ground_station)

                record_diff_time = end_time - start_time
                if record_diff_time > 3600*12 or record_diff_time < 0:
                    logging.info('SEVERE: Ground station/Station - {} {} Start to end of record longer than 12 hours {} {} ({} s)'.format(
                      corr_ground_station, orig_station,start_time,end_time,record_diff_time))
                    print('SEVERE: Ground station/Station - {} {} Start to end of record longer than 12 hours {} {} ({} s)'.format(
                      corr_ground_station, orig_station,start_time,end_time,record_diff_time))
                else:
                    logging.info('INFO: Ground station/Station -  {} {} Start to end of record {} {} ({} s)'.format(
                      corr_ground_station, orig_station,start_time,end_time,record_diff_time))
                
                df_list2 = split_midnight(df1)

                # write the files
                for df2 in df_list2:
                    # write the files
                    write_files(df2,gzip_filename,processed_dir)

                # add up the totals
                gst_total += len(df1)
                rec_final_total += len(df1)

                # make a warning? if there's a warning, make a single severe
                # tests on the start times 
                # +or- 24 hours of each other 
                # print out all the start - end times - number of records
                # not more than 4 hours long
                # SEVERE if they are 


            if len(df_list_bad) > 0:
                for df_bad in df_list_bad:
                # vertical_stack = pd.concat([survey_sub, survey_sub_last10], axis=0)
                    df_dropped = pd.concat([df_dropped, df_bad], ignore_index = False)

                rec_drop_local = len(df_dropped)
                rec_clock_flag_local = (df_dropped['clock_flag'] == 1).sum()
                if rec_dropped_total > 0:
                    if rec_dropped_total < 180:
                        logging.info('WARNING: Ground station/Station - {} {} - Dropped {} record(s)'.format(
                          corr_ground_station, orig_station,rec_dropped_total))
                    else: 
                        logging.info('SEVERE: Ground station/Station - {} {} - Dropped {} record(s)'.format(
                          corr_ground_station, orig_station,rec_dropped_total))
                        print('SEVERE: Ground station/Station - {} {} - Dropped {} record(s)'.format(
                          corr_ground_station, orig_station,rec_dropped_total))
                    if rec_clock_flag_local > 179:
                        logging.info('SEVERE: Ground station/Station - {} {} - {} clock flag record(s) have been dropped'.format(
                          corr_ground_station, orig_station,rec_dropped_total))
                        print('SEVERE: Ground station/Station - {} {} - {} clock flag record(s) have been dropped'.format(
                          corr_ground_station, orig_station,rec_dropped_total))

        for st, en, l, gr, sta in zip(start_times, end_times, lengths, final_ground_station, final_station):
            logging.info('START_TIMES: Ground station/Station - {} {} - {} {} {}'.format(gr, sta, str(st), str(en), l))

        logging.info('############################################')
        earlier = False
        for i, st in enumerate(start_times):
            if i > 0:
                if st < start_times[i-1]:
                    earlier = True
        if earlier == True:
            logging.info('SEVERE: File possibly contains records with the wrong timing')
        
        overall_start_times.extend(start_times)

        rec_dropped_total += len(df_dropped)
        write_dropped_files(df_dropped,gzip_filename,processed_dir)

        rec_clock_flag_total += (df_gst['clock_flag'] == 1).sum()

        rec_missing1 = df_gst['orig_no'].isna().sum()
        rec_missing_total += rec_missing1
        # rec_final_total += len(df_gst)
        if rec_missing1 > 0: 
            logging.info('WARNING: Ground station/Station - {} {} Inserted {} blank record(s)'.format(
              corr_ground_station, orig_station,rec_missing1))

        # logging.info('##INFORMATION: Ground station/Station - 7 S15########## Ground station: {} Station: {} {}'.format(corr_ground_station,orig_station,orig_timestamp))
        # 
        # logging.info('this is the end ')
        # logging.info(df_gst[(df_gst['orig_timestamp'] != df_gst['corr_timestamp']) & (~df_gst['orig_no'].isna())].to_string())
        # exit()



        rec_adjusted_timestamps_final = ((df_gst['orig_timestamp'] != df_gst['corr_timestamp']) & (~df_gst['orig_no'].isna())).sum()
        if rec_adjusted_timestamps_final > 0:
            logging.info('WARNING: Ground station/Station - {} {} Adjusted {} timestamp(s) (simple/complex)'.format(
              corr_ground_station, orig_station, rec_adjusted_timestamps_final))

        rec_adjusted_timestamps_final_total += rec_adjusted_timestamps_final



        # df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

    if config.initial:
        logging.info('############################################')
        logging.info('WARNING: Used initial=True')
        
    else:

        logging.info('############################################')

        earlier = False
        for i, st in enumerate(overall_start_times):
            if i > 0:
                if st < overall_start_times[i-1]:
                    earlier = True
        if earlier == True:
            logging.info('SEVERE: File (overall) possibly contains records with the wrong timing')

        if len(overall_start_times) >= 2:
            overall_start_times.sort()
            total_diff_time = overall_start_times[-1] - overall_start_times[0]
            if abs(total_diff_time) > 3600*24:
                logging.info('SEVERE: Ground station/Station - {} {}: This file is more than 24 hours long: {} s ({:01f}) {} {}'.format(
                  corr_ground_station, orig_station,total_diff_time,total_diff_time/3600,overall_start_times[0],overall_start_times[-1]))

        rec_damaged_total = (rec_deleted_timestamps_total + rec_duplicates + 
          rec_station_duplicates + rec_drop_damaged_total + rec_dropped_total)
        rec_valid_total = rec_final_total - rec_missing_total

        if rec_valid_total != 0:
            rec_damaged_percentage = 100*(rec_damaged_total/rec_valid_total)
        else: 
            rec_damaged_percentage = 100
        if rec_valid_total != 0:
            rec_missing_percentage = 100*(rec_missing_total/rec_valid_total)
        else:
            rec_missing_percentage = 100
        rec_valid_of_intial_percentage = 100*(rec_valid_total/rec_initial_total)

        logging.info('############################################')

        config.extra_ground_stations.sort()

        if len(config.extra_ground_stations) > 0:
            logging.info('INFO: Final Check - Extra ground stations added {}'.format(config.extra_ground_stations))

        if rec_clock_flag_total > 0:
            logging.info('WARNING: Final Check - {} clock flag(s)'.format(rec_clock_flag_total))

        if rec_negative_total > 0:
            logging.info('Negative gaps: {}'.format(
                 rec_negative_total))

        if rec_missing_total > 0: 
            logging.info('WARNING: Final Check - Inserted {} blank record(s) ({:0.02f}%)'.format(
              rec_missing_total,rec_missing_percentage))


        if rec_dropped_total > 0: 
            logging.info('WARNING: Final Check - Dropped {} record(s)'.format(
              rec_dropped_total))

        if rec_damaged_sync_total > 0:
            logging.info('WARNING: Final Check -  Removed {} damaged sync code(s)'.format(
              rec_damaged_sync_total))

        if rec_damaged_total > 0:
            logging.info('WARNING: Final Check -  Removed {} damaged record(s) ({:0.02f}%) known damaged {} duplicates {} station duplicates {} deleted timestamp(s) {} dropped record(s) {}'.format(
              rec_damaged_total,rec_damaged_percentage,rec_drop_damaged_total,rec_duplicates,rec_station_duplicates,
              rec_deleted_timestamps_total, rec_dropped_total))

        if rec_fixed_frame_total > 0:
            logging.info('WARNING: Final Check -  Simple adjustments to {} frame number(s)'.format(
                rec_fixed_frame_total))

        # if rec_fixed_simple_timestamp_total > 0:
        #     logging.info('WARNING: Final Check -  Simple adjustments to {} timestamp(s)'.format(
        #         rec_fixed_simple_timestamp_total))

        # using the final total, which is calculated from the difference 
        # between corr_timestamp and orig_timestamp
        if rec_adjusted_timestamps_final_total > 0:
            logging.info('WARNING: Final Check -  Adjusted {} timestamp(s) (simple/complex)'.format(
                rec_adjusted_timestamps_final_total))

        if rec_deleted_non_uniq_total > 0:
            logging.info('WARNING: Final Check -  {} Non unique record(s) removed'.format(
                rec_deleted_non_uniq_total))

        if single_station is not None or single_ground_station is not None:
            logging.info('FINAL: Processed only {} and {} '.format(single_station, single_ground_station))

        logging.info('FINAL: Number of records: {}'.format(rec_final_total))
        if rec_valid_of_intial_percentage < 95:
            logging.info('FINAL: WARNING Number of valid records: {} ({:0.02f}% of initial)'.format(rec_valid_total, rec_valid_of_intial_percentage))
        else: 
            logging.info('FINAL: Number of valid records: {} ({:0.02f}% of initial)'.format(rec_valid_total, rec_valid_of_intial_percentage))

        rec_percent_sync = 100*(rec_damaged_sync_total/rec_initial_total)
        rec_percent_damaged = 100*(rec_damaged_total/rec_initial_total)


        # making a distinction between sync damaged and other damaged - because there might be some potential to include the other damaged
        
        logging.info('FINAL: Statistics Initial={} Final={} Sync errors={} Percent sync errors={:0.02f}% Dropped due to timing={} Percent dropped={:0.02f}% Inserted={} Final Valid={} Final Valid Percent={:0.02f}%'.format(
            rec_initial_total, 
            rec_final_total, 
            rec_damaged_sync_total, 
            rec_percent_sync,
            rec_damaged_total, 
            rec_percent_damaged, 
            rec_missing_total, 
            rec_valid_total, 
            rec_valid_of_intial_percentage))

# FINAL: Number of records: 110581
# WARNING: Final Check - Inserted 2386 blank record(s) (2.21%)
# WARNING: Final Check -  Removed 767 damaged record(s) (0.71%) - sync code 767 known damaged 0 duplicates 0 station duplicates 0 deleted timestamp(s) 0
# WARNING: Final Check -  Simple adjustments to 1 frame number(s)
# WARNING: Final Check -  Adjusted 3 timestamp(s) (simple/complex)

        # logging.info('total {}'.format(total))
        # logging.info('gst total {}'.format(gst_total))

def write_files(df,gzip_filename,processed_dir):
    # wtn.17.19.csv.gz
        
    if len(df) > 0: 
        corr_ground_station = df['corr_ground_station'].iloc[0]
        orig_station = df['orig_station'].iloc[0]
        starttime = df['corr_timestamp'].iloc[0]
        julday = str('{:03}'.format(starttime.dayofyear))
        year = starttime.year
        str_starttime = starttime.strftime('%H_%M_%S_%f')
        # xa.s12.01.afr.1976.066.1.0.mseed
        out_gzip_filename = gzip_filename.replace('.csv.gz', '')
        out_gzip_filename = '{}.{}.{}.{}.{}.{}.csv.gz'.format(out_gzip_filename,orig_station,corr_ground_station,year,julday,str_starttime)
        out_gzip_filename = os.path.basename(out_gzip_filename)
        out_gzip_filename = os.path.join(processed_dir,out_gzip_filename)
        
        print('Writing File ', out_gzip_filename)
        df.to_csv(out_gzip_filename,index=False,date_format='%Y-%m-%dT%H:%M:%S.%fZ',quoting=csv.QUOTE_NONNUMERIC)

def split_midnight(df_gst):
    # if the datafile crosses midnight, split it  

    #TODO Note that if there are no valid values around midnight, but the 
    #dataframe does have values at midnight, the records can start or end with 
    #the empty ones, which isn't ideal 
    

    df_gst['days'] = df_gst.corr_timestamp.dt.normalize()
    gb = df_gst.groupby('days')    
    df_list = [gb.get_group(x) for x in gb.groups]
    return df_list

def default_config():
    config.extra_ground_stations = []
    config.station_order = []
    config.last_station=None
    config.last_timestamp=None

    # the other parameters can be set in the calling file, 
    # and shouldn't be reset here




def write_dropped_files(df_dropped,gzip_filename,processed_dir):
    # wtn.17.19.csv.gz

    if (df_dropped is not None) and (len(df_dropped) > 0): 
        df_dropped.sort_values(by=['orig_idx'], inplace=True)
        corr_ground_station = df_dropped['corr_ground_station'].iloc[0]
        orig_station = df_dropped['orig_station'].iloc[0]
        starttime = df_dropped['corr_timestamp'].iloc[0]
        julday = str('{:03}'.format(starttime.dayofyear))
        year = starttime.year
        str_starttime = starttime.strftime('%H_%M_%S_%f')
        out_gzip_filename = gzip_filename.replace('.csv.gz', '')
        # out_gzip_filename = '{}.{}.{}.dropped.csv.gz'.format(out_gzip_filename,orig_station,corr_ground_station)
        out_gzip_filename = '{}.{}.{}.{}.{}.{}.dropped.csv.gz'.format(out_gzip_filename,orig_station,corr_ground_station,year,julday,str_starttime)
        out_gzip_filename = os.path.basename(out_gzip_filename)
        out_gzip_filename = os.path.join(processed_dir,out_gzip_filename)
        df_dropped.to_csv(out_gzip_filename,index=False,date_format='%Y-%m-%dT%H:%M:%S.%fZ',quoting=csv.QUOTE_NONNUMERIC)
        logging.info('WARNING: Ground station/Station - {} {} - Dropped {} records(s)'.format(corr_ground_station,orig_station,len(df_dropped)))
        # logging.debug('Dropped files:\n{}'.format(df_dropped.to_string()))


def initial_cleanup(df):

    df['orig_timestamp'] = df['orig_timestamp'].astype('datetime64[ns, UTC]')


    # main tapes extracted without the S 
    df['orig_station'].replace('11', 'S11',inplace=True)
    df['orig_station'].replace('12', 'S12',inplace=True)
    df['orig_station'].replace('14', 'S14',inplace=True)
    df['orig_station'].replace('15', 'S15',inplace=True)
    df['orig_station'].replace('16', 'S16',inplace=True)
    df['orig_station'] = df['orig_station'].astype('string')

    if 'bit_synchronizer' in df.columns:
        df['bit_synchronizer'] = df['bit_synchronizer'].astype('string')
    else:
        df['bit_synchronizer'] = pd.NA
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
    df['frame'] = to_Int64(df['frame'])

    if 'shz_4' in df.columns:
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
        # shz 46 doesn't work for all stations, so replace with null
        # if necessary
        if 'shz_46' not in df.columns:
            df['shz_46'] = pd.NA
        df['shz_46'] = to_Int64(df['shz_46'])        
        df['shz_48'] = to_Int64(df['shz_48'])
        df['shz_50'] = to_Int64(df['shz_50'])
        df['shz_52'] = to_Int64(df['shz_52'])
        df['shz_54'] = to_Int64(df['shz_54'])
        df['shz_58'] = to_Int64(df['shz_58'])
        df['shz_60'] = to_Int64(df['shz_60'])
        df['shz_62'] = to_Int64(df['shz_62'])
        df['shz_64'] = to_Int64(df['shz_64'])
    else: 
        logging.info('WARNING: No SHZ values found.')    
    

    # replace empty ground station and orig station with values
    df['orig_ground_station'].fillna(-1,inplace=True)
    df['orig_station'].fillna('S0',inplace=True)
    df['orig_frame'].fillna(-99,inplace=True)
    df['frame'].fillna(-99,inplace=True)


    # add new columns for timestamps
    df['corr_timestamp'] = df['orig_timestamp']   
    df['corr_gap'] = np.NaN
    df['frame_change'] = 1
    df['corr_gap_count'] = 1
    df['corr_ground_station'] = df['orig_ground_station']
    df['cumsum'] = 0

    # df['begin'] = False
    # df['end'] = False

    df['begin_end'] = 'None'
    df['delta4'] = pd.NA
    # df['guess_count'] = pd.NA    

    df['orig_idx'] = df.index

    return df

def drop_last_frame(df):
    # drop the last frame if it is empty
    # if df['orig_no'].iloc[-60] == 4369 and pd.isna(df['orig_station'].iloc[-60]):
    if df['orig_no'].iloc[-60] == 4369 and (df['orig_station'].iloc[-60] == 'S0'):
        df = df[:-60]
        logging.info('INFO: Last frame removed because it was empty')
    else:
        logging.info('WARNING: Last frame was not empty (record will probably continue to next csv file)')

    return df

def initial_report(df):

    global rec_initial_total
    logging.info('Total records: {}'.format(len(df)))
    rec_initial_total = (df['orig_station'] != 'S17').sum()
    logging.info('Total records (excluding S17): {}'.format(rec_initial_total))

    df_filter = df[df['orig_station'] != 'S17']
    df_filter2 = df_filter[df_filter['sync'] != '1110001001000011101101']
    
    damaged_records = len(df_filter2)
    if damaged_records > 0:
        logging.info('WARNING: Data Check - {} ({:0.02f}%) damaged record(s) - sync code is incorrect'.format(damaged_records, (100*damaged_records/rec_initial_total)))
    # if damaged_records > 200:
    # 
    #     logging.info('WARNING: Data Check - {} ({:0.02f}%) damaged record(s) - sync code is incorrect'.format(damaged_records, (100*damaged_records/rec_initial_total)))
    #     raise Exception



def detailed_report(df):
    # logging.info(df['corr_ground_station'].dropna().values.astype('int').tolist())
    groups = [ list(group) for key, group in groupby(df['corr_ground_station'].dropna().values.astype('int').tolist())]

    logging.info('Ground station, no_of_records')
    for group in groups:
        logging.info('{}, {}'.format(int(group[0]), len(group)))

    missing_timestamps = (df['orig_timestamp'].isna()).sum()
    if missing_timestamps > 0:
        logging.info('WARNING: Data Check - {} missing timestamp(s)'.format(missing_timestamps))
    else: 
        logging.info('PASSED: Data Check - {} missing timestamp(s)'.format(missing_timestamps))

    # missing_frames = (df['frame'] == -99).sum()
    # if missing_frames > 0:
    #     logging.info('WARNING: Data Check - {} missing frame(s)'.format(missing_frames))

    missing_clock_flags = df['clock_flag'].isna().sum()
    if missing_clock_flags > 0:
        logging.info('WARNING: Data Check - {} missing clock flag(s)'.format(missing_clock_flags))

    missing_stations = (df['orig_station'] == 'S0').sum()
    if missing_stations > 0:
        logging.info('WARNING: Data Check - {} missing station(s)'.format(missing_stations))
    else: 
        logging.info('PASSED: Data Check - {} missing station(s)'.format(missing_stations))

    missing_ground_stations = (df['corr_ground_station'] == -1).sum()
    if missing_ground_stations > 0:
        logging.info('WARNING: Data Check - {} missing ground station(s)'.format(missing_ground_stations))
    else: 
        logging.info('PASSED: Data Check - {} missing ground station(s)'.format(missing_ground_stations))

    set_clock_flag = (df['clock_flag'] == 1).sum()
    if set_clock_flag > 0:
        logging.info('WARNING: Data Check - {} clock flag(s)'.format(set_clock_flag))
    else: 
        logging.info('PASSED: Data Check - {} clock flag(s)'.format(set_clock_flag))

    global first_record
    first_record = df['orig_timestamp'].iloc[0:5].min()
    global last_record 
    last_record = df['orig_timestamp'].iloc[-5:].max()

    logging.info('First record: {} Last record: {}'.format(first_record, last_record))
    logging.info('First record: {} Last record: {}'.format(UTCDateTime(first_record), last_record))

    elapsed = last_record-first_record
    if elapsed > pd.Timedelta(hours=12):
        logging.info('WARNING: Time Elapsed: {}'.format(elapsed)) 
    else:
        logging.info('Time Elapsed: {}'.format(elapsed))

def all_known_errors(df,gzip_filename):

    rec_drop_damaged_total = 0

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.1.3.csv.gz', 
              orig_idx_start=337105, orig_idx_end=337250,single_station='S12')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.1.41.csv.gz', 
              start_timestamp ='1976-03-18T21:56:43.634000Z', 
              end_timestamp='1976-01-06T08:37:33.856000Z',
              start_station='S12',
              end_station=pd.NA,
              corr_ground_station1=1,
              corr_ground_station2=pd.NA)
    rec_drop_damaged_total += rec_drop_damaged1


    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.1.43.csv.gz', 
              start_timestamp ='1976-06-09T13:08:47.744000Z', 
              end_timestamp='1976-03-20T09:15:41.240000Z',
              start_station=pd.NA,
              end_station='S17',
              corr_ground_station1=4,
              corr_ground_station2=4,
              start_orig_no=16,
              end_orig_no=16)

    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.10.19.csv.gz', 
              start_timestamp ='1976-08-01T14:55:50.865000Z', 
              end_timestamp='1976-11-30T00:00:51.674000Z',
              start_station='S15',
              end_station='S14',
              corr_ground_station1=10,
              corr_ground_station2=1,
              start_orig_no=3395,
              end_orig_no=57344)

    rec_drop_damaged_total += rec_drop_damaged1

# 70      3395 1976-08-01 14:55:50.865000+00:00          11          S15           1                   10 

    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.10.24.csv.gz', 
              start_timestamp ='1976-12-02T11:58:04.049000Z', 
              end_timestamp='1976-01-06T09:12:31.008000Z',
              start_station='S12',
              end_station=pd.NA,
              corr_ground_station1=6,
              corr_ground_station2=pd.NA,
              start_orig_no=846,
              end_orig_no=1792

# begin: 
# 846,"1976-12-02T11:58:04.049000Z",1,"S12",0,6,53
# end: 
# 1792,"1976-01-06T09:12:31.008000Z",48,"",0,"",8,136,968,8

)
    rec_drop_damaged_total += rec_drop_damaged1
    # new



    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.10.33.csv.gz', 
              start_timestamp ='1976-01-01T18:01:01.795000Z', 
              end_timestamp='1976-12-07T18:59:16.887000Z',
              start_station=pd.NA,
              end_station='S14',
              corr_ground_station1=pd.NA,
              corr_ground_station2=9,
              start_orig_no=3787,
              end_orig_no=48
)
# beginning of bad:
# 3787,"1976-01-01T18:01:01.795000Z",23,"",1,"",584,518,448,584,582,512,520,518,0,520,582,576,"00000",3,"1
# end of bad:
# 48,"1976-12-07T18:59:16.887000Z",60,"S14",0,9,499,452,519,500,452,519,500,452,519,500,452,519,"00011",59,

    rec_drop_damaged_total += rec_drop_damaged1


    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.10.39.csv.gz', 
              start_timestamp ='1976-12-11T09:14:44.010000Z', 
              end_timestamp='1976-03-09T21:56:27.395000Z',
              start_station='S12',
              end_station=pd.NA,
              corr_ground_station1=6,
              corr_ground_station2=pd.NA)
    rec_drop_damaged_total += rec_drop_damaged1


    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.10.55.csv.gz', 
              start_timestamp ='1976-12-21T03:53:04.950000Z', 
              end_timestamp='1976-12-21T03:52:56.107000Z',
              start_station='S15',
              end_station='S14',
              corr_ground_station1=9,
              corr_ground_station2=9,
              start_orig_no=339,
              end_orig_no=16)
    rec_drop_damaged_total += rec_drop_damaged1

    
    # 10.6
    # a lot of problems with 10.15 (mainly deleted)

# 339,"1976-12-21T03:53:04.950000Z",58,"S15",0,9,505,500,491,505,500,491,504,500,491,519,263,711,"00011",47,

# beginning of bad:
# 16,"1976-06-01T00:22:59.584000Z",1,"S12",0,4,576,576,576,576,587,583,512,576,576,576,7,490,"00100",35,
# 
# end of bad:
# 16,"1976-12-21T03:52:56.107000Z",60,"S14",0,9,511,490,509,511,489,509,511,489,508,511,489,509,"00011",17

    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.10.55.csv.gz', 
              start_timestamp ='1976-01-06T07:45:08.128000Z', 
              end_timestamp='1976-12-21T05:26:21.885000',
              start_station=pd.NA,
              end_station='S14',
              corr_ground_station1=pd.NA,
              corr_ground_station2=9,
              start_orig_no=1152,
              end_orig_no=1152)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.11.18.csv.gz', 
              orig_idx_start=101796, orig_idx_end=101838,single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.11.24.csv.gz', 
              orig_idx_start=211596, orig_idx_end=211630,single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1

    # orig_no=1507
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.11.25.csv.gz', 
              orig_idx_start=60, orig_idx_end=120)
    rec_drop_damaged_total += rec_drop_damaged1
    # 

    # orig_no=49152
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.11.25.csv.gz', 
              orig_idx_start=120, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1


    # orig_no=40960
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.11.41.csv.gz', 
              orig_idx_start=120, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.11.44.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.11.48.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.11.55.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.11.8.csv.gz', 
              orig_idx_start=101796, orig_idx_end=101838,single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.12.34.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.12.40.csv.gz', 
              orig_idx_start=93406, orig_idx_end=93419)
    rec_drop_damaged_total += rec_drop_damaged1

    # POSSIBLE COPIED STATION ERROR
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.13.1.csv.gz', 
              orig_idx_start=37740, orig_idx_end=40639)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.13.12.csv.gz',
              orig_idx_start=186708, orig_idx_end=186768,single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1

    # COPIED STATION ERROR
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.13.3.csv.gz', 
              orig_idx_start=0, orig_idx_end=21179,single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.14.16.csv.gz', 
              orig_idx_start=167335, orig_idx_end=167381)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.14.44.csv.gz', 
              orig_idx_start=70, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.14.55.csv.gz', 
              orig_idx_start=60, orig_idx_end=239)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.17.csv.gz', 
              orig_idx_start=237224, orig_idx_end=237257,single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.17.csv.gz', 
              orig_idx_start=237849, orig_idx_end=237882,single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.18.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.21.csv.gz', 
              orig_idx_start=68, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.3.csv.gz', 
              orig_idx_start=37560, orig_idx_end=37619)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.3.csv.gz', 
              orig_idx_start=125400, orig_idx_end=125459)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.30.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.30.csv.gz', 
              orig_idx_start=246600, orig_idx_end=246599)
    rec_drop_damaged_total += rec_drop_damaged1

    # a lot of damaged records 
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.37.csv.gz', 
              orig_idx_start=214837, orig_idx_end=239955)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.15.50.csv.gz', 
              orig_idx_start=352320, orig_idx_end=352434)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.16.13.csv.gz', 
              orig_idx_start=34440, orig_idx_end=34514)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.16.15.csv.gz', 
              orig_idx_start=129128, orig_idx_end=129161,single_station='S16')
    rec_drop_damaged_total += rec_drop_damaged1




    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.16.26.csv.gz', 
              orig_idx_start=68, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.16.26.csv.gz', 
              orig_idx_start=257100, orig_idx_end=257159)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.16.40.csv.gz', 
              orig_idx_start=250184, orig_idx_end=250259)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.16.44.csv.gz', 
              orig_idx_start=82, orig_idx_end=239)
    rec_drop_damaged_total += rec_drop_damaged1

    # df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.17.47.csv.gz', 
    #           orig_idx_start=51614, orig_idx_end=51726)
    # rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz', 
              orig_idx_start=83134, orig_idx_end=83219)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz', 
              orig_idx_start=83347, orig_idx_end=83519)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz', 
              orig_idx_start=83624, orig_idx_end=83700)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz', 
              orig_idx_start=83831, orig_idx_end=83940)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.17.47.csv.gz', 
              orig_idx_start=32236, orig_idx_end=32249, single_station='S12')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.17.47.csv.gz', 
              orig_idx_start=51611, orig_idx_end=51760,single_station='S12')
    rec_drop_damaged_total += rec_drop_damaged1

    # the frame counts are not always working
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.17.47.csv.gz', 
              orig_idx_start=42569, orig_idx_end=79022,single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1

    # df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.18.48.csv.gz', 
    #           orig_idx_start=39860, orig_idx_end=39940)
    # rec_drop_damaged_total += rec_drop_damaged1
    # 
    # df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.18.48.csv.gz', 
    #           orig_idx_start=37095, orig_idx_end=52125,single_station='S12')
    # rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.18.48.csv.gz', 
              orig_idx_start=37095, orig_idx_end=52125)
    rec_drop_damaged_total += rec_drop_damaged1


    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.17.9.csv.gz', 
              orig_idx_start=60, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.18.16.csv.gz', 
              orig_idx_start=44640, orig_idx_end=44759)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.18.19.csv.gz', 
              orig_idx_start=239677, orig_idx_end=239879)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.18.25.csv.gz', 
              orig_idx_start=173003, orig_idx_end=173073)
    rec_drop_damaged_total += rec_drop_damaged1


    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.18.32.csv.gz', 
              orig_idx_start=271680, orig_idx_end=271799)
    rec_drop_damaged_total += rec_drop_damaged1


    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.18.34.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.18.46.csv.gz', 
              orig_idx_start=191655, orig_idx_end=191759)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.19.18.csv.gz', 
              orig_idx_start=60, orig_idx_end=239)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.19.35.csv.gz', 
              orig_idx_start=60, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.19.46.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.19.50.csv.gz', 
              orig_idx_start=280507, orig_idx_end=280679)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.19.50.csv.gz', 
              orig_idx_start=143160, orig_idx_end=143279)
    rec_drop_damaged_total += rec_drop_damaged1

    # XFAIL
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.55.csv.gz',start_timestamp ='1977-08-26T02:37:26.628000Z',end_timestamp='1977-08-26T03:44:54.689000Z',start_station='S12',end_station='S17',start_orig_no=607,end_orig_no=1165,orig_ground_station1=2)
# 607,"1977-08-26T02:37:26.628000Z",1,"S12",0,2,514,513,513,514,513,513,514,513,513,514,513,513,"00011",20,"1110001001000011101101"
# 1165,"1977-08-26T03:44:54.689000Z",60,"S17",0,2,0,0,0,514,508,0,0,0,0,0,0,0,"00011",61,"1110001001000011101101"

    # XFAIL
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.19.6.csv.gz', 
              orig_idx_start=106265, orig_idx_end=151017,single_station='S15')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.2.29.csv.gz', 
              orig_idx_start=70, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.2.29.csv.gz', 
              orig_idx_start=16860, orig_idx_end=16919)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.2.29.csv.gz', 
              orig_idx_start=186632, orig_idx_end=186751)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.2.51.csv.gz', 
              orig_idx_start=149368, orig_idx_end=149387)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.2.52.csv.gz', 
              orig_idx_start=204788, orig_idx_end=204899)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.20.15.csv.gz', 
              orig_idx_start=59, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.20.27.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.20.29.csv', 
              orig_idx_start=263110, orig_idx_end=263155)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.20.30.csv', 
              orig_idx_start=97167, orig_idx_end=97541,single_station='S15')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.20.6.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.21.37.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.3.15.csv', 
              orig_idx_start=161649, orig_idx_end=161707,single_station='S12')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.3.5.csv.gz', 
              orig_idx_start=17209, orig_idx_end=17282)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.3.54.csv.gz', 
              orig_idx_start=17099, orig_idx_end=17159)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.3.54.csv.gz', 
              orig_idx_start=92736, orig_idx_end=92855)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.3.55.csv.gz', 
              orig_idx_start=68, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.3.55.csv.gz', 
              orig_idx_start=15596, orig_idx_end=15660)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.4.10.csv', 
              orig_idx_start=278412, orig_idx_end=278520,single_station='S15')
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.4.14.csv.gz', 
              orig_idx_start=68, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.4.50.csv.gz', 
              orig_idx_start=240531, orig_idx_end=240655)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.4.52.csv.gz', 
              orig_idx_start=83, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.6.20.csv.gz', 
              orig_idx_start=364950, orig_idx_end=369719)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.5.12.csv.gz', 
              orig_idx_start=67, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.5.26.csv.gz', 
              orig_idx_start=73, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.5.35.csv.gz', 
              orig_idx_start=61427, orig_idx_end=61559)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.5.50.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.6.24.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.6.52.csv.gz', 
              orig_idx_start=60, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.7.31.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.7.40.csv.gz', 
              orig_idx_start=69, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.8.17.csv.gz', 
              orig_idx_start=67, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.8.46.csv.gz', 
              orig_idx_start=70500, orig_idx_end=70508)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.8.50.csv.gz', 
              orig_idx_start=85659, orig_idx_end=85859)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.8.53.csv.gz', 
              orig_idx_start=12179, orig_idx_end=29459)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.8.54.csv.gz', 
              orig_idx_start=23520, orig_idx_end=23579)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.9.13.csv.gz', 
              orig_idx_start=208339, orig_idx_end=208460)
    rec_drop_damaged_total += rec_drop_damaged1

    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.9.29.csv.gz', 
              orig_idx_start=59, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    # XFAIL
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.9.50.csv.gz', 
              orig_idx_start=66, orig_idx_end=179)
    rec_drop_damaged_total += rec_drop_damaged1

    old_timestamp ='1976-12-28 12:48:53.375000Z'
    # add two gaps 
    new_timestamp=pd.Timestamp(old_timestamp) + pd.Timedelta(seconds=(DELTA*4*2))
    df = timestamp_adjust(df, gzip_filename, problem_gzip_filename='wtn.11.13.csv.gz', 
              old_timestamp=old_timestamp, 
              # add two gaps 
              new_timestamp=new_timestamp,
              station1='S15',
              corr_ground_station1=9)

    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.11.7.csv.gz', 
              start_timestamp ='1976-12-24 16:21:59.664000', 
              end_timestamp='1976-12-24 16:21:59.664000',
              start_station='S12',
              end_station='S12',
              corr_ground_station1=9,
              corr_ground_station2=9)
    rec_drop_damaged_total += rec_drop_damaged1

# 2201,"1976-12-24T16:21:59.664000Z",59,"S12",0,9,513,510,472,513,510,472,513,510,471,514,510,471,"00011",39,"1110001001000011101101"

    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.18.18.csv.gz', 
              start_timestamp ='1977-07-13T11:24:03.539000Z', 
              end_timestamp='1977-07-13T11:24:08.368000Z',
              start_station='S14',
              end_station='S14',
              corr_ground_station1=8,
              single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1

    old_timestamp ='1977-07-13T13:11:13.310000Z'
    # use the estimated gap 
    new_timestamp =pd.Timestamp('1977-07-13T13:11:14.343092Z')
    df = timestamp_adjust(df, gzip_filename, problem_gzip_filename='wtn.18.18.csv.gz', 
              old_timestamp=old_timestamp, 
              # add two gaps 
              new_timestamp=new_timestamp,
              station1='S14',
              corr_ground_station1=8)

    df, rec_drop_damaged1 = all_drop_data(df, gzip_filename, problem_gzip_filename='wtn.12.40.csv.gz', 
              start_timestamp ='1977-02-10 18:05:48.726000+00:00', 
              end_timestamp='1977-02-10 18:06:56.947000+00:00',
              start_station='S14',
              end_station='S14',
              corr_ground_station1=5,
              corr_ground_station2=5,
              single_station='S14')
    rec_drop_damaged_total += rec_drop_damaged1


    # New, Improved!
    # XFAIL
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.10.csv.gz',start_timestamp ='1976-03-05 01:09:47.019000+00:00',end_timestamp='1976-03-05 03:39:54.226000+00:00',start_station='S12',end_station='S17',start_orig_no=58358,end_orig_no=3970,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.24.csv.gz',start_timestamp ='1976-03-10 20:00:22.091000+00:00',end_timestamp='1976-03-10 20:00:24.217000+00:00',start_station='S12',end_station='S17',start_orig_no=61,end_orig_no=61,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.25.csv.gz',start_timestamp ='1976-03-11 07:36:01.552000+00:00',end_timestamp='1976-03-11 07:36:03.033000+00:00',start_station='S12',end_station='S17',start_orig_no=62,end_orig_no=62,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.27.csv.gz',start_timestamp ='1976-03-12 12:53:02.688000+00:00',end_timestamp='1976-03-12 13:26:53.010000+00:00',start_station='S12',end_station='S17',start_orig_no=2715,end_orig_no=2993,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.32.csv.gz',start_timestamp ='1976-03-14 14:19:25.645000+00:00',end_timestamp='1976-03-14 17:27:54.277000+00:00',start_station='S12',end_station='S17',start_orig_no=2787,end_orig_no=4347,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.34.csv.gz',start_timestamp ='1976-03-15 18:07:08.333000+00:00',end_timestamp='1976-03-15 19:47:57.090000+00:00',start_station='S12',end_station='S17',start_orig_no=3045,end_orig_no=3879,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.35.csv.gz',start_timestamp ='1976-03-16 06:27:15.047000+00:00',end_timestamp='1976-03-16 07:38:31.934000+00:00',start_station='S12',end_station='S17',start_orig_no=1848,end_orig_no=2436,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.4.csv.gz',start_timestamp ='1976-03-02 03:02:33.424000+00:00',end_timestamp='1976-03-02 04:24:53.801000+00:00',start_station='S12',end_station='S17',start_orig_no=689,end_orig_no=1370,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.4.csv.gz',start_timestamp ='1976-03-02 03:21:59.866000+00:00',end_timestamp='1976-03-02 04:24:53.801000+00:00',start_station='S12',end_station='S17',start_orig_no=850,end_orig_no=1370,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.4.csv.gz',start_timestamp ='1976-03-02 03:30:55.996000+00:00',end_timestamp='1976-03-02 04:24:53.801000+00:00',start_station='S12',end_station='S17',start_orig_no=924,end_orig_no=1370,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.4.csv.gz',start_timestamp ='1976-03-02 04:03:10.407000+00:00',end_timestamp='1976-03-02 04:24:53.801000+00:00',start_station='S12',end_station='S17',start_orig_no=1191,end_orig_no=1370,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.4.csv.gz',start_timestamp ='1976-03-02 04:23:59.556000+00:00',end_timestamp='1976-03-02 08:13:52.998000+00:00',start_station='S12',end_station='S17',start_orig_no=1371,end_orig_no=3301,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.4.csv.gz',start_timestamp ='1976-03-02 05:05:10.094000+00:00',end_timestamp='1976-03-02 08:13:52.998000+00:00',start_station='S12',end_station='S17',start_orig_no=1712,end_orig_no=3301,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.42.csv.gz',start_timestamp ='1976-03-19 16:45:42.623000+00:00',end_timestamp='1976-03-19 18:27:53.546000+00:00',start_station='S12',end_station='S17',start_orig_no=3189,end_orig_no=3865,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.5.csv.gz',start_timestamp ='1976-03-02 15:15:44.324000+00:00',end_timestamp='1976-03-02 19:10:52.552000+00:00',start_station='S12',end_station='S17',start_orig_no=2134,end_orig_no=4079,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.5.csv.gz',start_timestamp ='1976-03-02 20:07:57.172000+00:00',end_timestamp='1976-03-02 23:59:05.468000+00:00',start_station='S12',end_station='S17',start_orig_no=3,end_orig_no=1915,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.5.csv.gz',start_timestamp ='1976-03-02 22:34:34.379000+00:00',end_timestamp='1976-03-02 23:59:05.468000+00:00',start_station='S12',end_station='S17',start_orig_no=1216,end_orig_no=1915,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.5.csv.gz',start_timestamp ='1976-03-02 23:35:47.581000+00:00',end_timestamp='1976-03-02 23:59:05.468000+00:00',start_station='S12',end_station='S17',start_orig_no=1723,end_orig_no=1915,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.50.csv.gz',start_timestamp ='1976-03-23 21:35:06.298000+00:00',end_timestamp='1976-03-23 22:41:52.084000+00:00',start_station='S12',end_station='S17',start_orig_no=2864,end_orig_no=1734,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.52.csv.gz',start_timestamp ='1976-03-24 16:04:35.097000+00:00',end_timestamp='1976-03-24 17:23:58.242000+00:00',start_station='S12',end_station='S17',start_orig_no=2490,end_orig_no=3015,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.53.csv.gz',start_timestamp ='1976-03-24 22:58:01.280000+00:00',end_timestamp='1976-03-25 02:40:56.458000+00:00',start_station='S12',end_station='S17',start_orig_no=2247,end_orig_no=3723,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.55.csv.gz',start_timestamp ='1976-03-25 23:57:59.087000+00:00',end_timestamp='1976-03-26 03:57:58.306000+00:00',start_station='S12',end_station='S17',start_orig_no=2439,end_orig_no=4028,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.6.csv.gz',start_timestamp ='1976-03-03 04:11:44.780000+00:00',end_timestamp='1976-03-03 04:40:57.425000+00:00',start_station='S12',end_station='S17',start_orig_no=1814,end_orig_no=2055,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.6.csv.gz',start_timestamp ='1976-03-03 12:15:28.760000+00:00',end_timestamp='1976-03-03 14:16:04.142000+00:00',start_station='S12',end_station='S17',start_orig_no=1,end_orig_no=998,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.6.csv.gz',start_timestamp ='1976-03-03 12:38:08.400000+00:00',end_timestamp='1976-03-03 14:16:04.142000+00:00',start_station='S12',end_station='S17',start_orig_no=190,end_orig_no=998,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.7.csv.gz',start_timestamp ='1976-03-03 18:10:47.101000+00:00',end_timestamp='1976-03-03 20:19:58.553000+00:00',start_station='S12',end_station='S17',start_orig_no=2974,end_orig_no=4043,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.7.csv.gz',start_timestamp ='1976-03-03 18:11:01.590000+00:00',end_timestamp='1976-03-03 20:19:58.553000+00:00',start_station='S12',end_station='S17',start_orig_no=2976,end_orig_no=4043,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.7.csv.gz',start_timestamp ='1976-03-03 22:42:47.610000+00:00',end_timestamp='1976-03-04 02:11:41.268000+00:00',start_station='S12',end_station='S17',start_orig_no=1221,end_orig_no=2950,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.8.csv.gz',start_timestamp ='1976-03-04 02:46:35.204000+00:00',end_timestamp='1976-03-04 03:41:52.798000+00:00',start_station='S12',end_station='S12',start_orig_no=3240,end_orig_no=3697,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.8.csv.gz',start_timestamp ='1976-03-04 02:53:28.167000+00:00',end_timestamp='1976-03-04 03:41:52.798000+00:00',start_station='S12',end_station='S12',start_orig_no=3297,end_orig_no=3697,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.8.csv.gz',start_timestamp ='1976-03-04 03:16:24.712000+00:00',end_timestamp='1976-03-04 03:41:52.798000+00:00',start_station='S12',end_station='S12',start_orig_no=3487,end_orig_no=3697,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.8.csv.gz',start_timestamp ='1976-03-04 08:11:12.711000+00:00',end_timestamp='1976-03-04 09:24:59.350000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=2811,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.14.csv.gz',start_timestamp ='1976-11-27 15:45:38.970000+00:00',end_timestamp='1976-11-27 16:19:55.598000+00:00',start_station='S12',end_station='S14',start_orig_no=2176,end_orig_no=1487,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.15.csv.gz',start_timestamp ='1976-11-28 09:20:40.256000+00:00',end_timestamp='1976-11-28 12:17:51.469000+00:00',start_station='S12',end_station='S14',start_orig_no=2134,end_orig_no=3307,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.15.csv.gz',start_timestamp ='1976-11-28 14:04:52.328000+00:00',end_timestamp='1976-11-28 14:04:54.716000+00:00',start_station='S12',end_station='S14',start_orig_no=806,end_orig_no=806,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.20.csv.gz',start_timestamp ='1976-11-30 14:20:30.518000+00:00',end_timestamp='1976-11-30 14:34:52.984000+00:00',start_station='S12',end_station='S14',start_orig_no=125,end_orig_no=218,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.23.csv.gz',start_timestamp ='1976-11-26 20:37:28.464000+00:00',end_timestamp='1976-11-26 20:37:30.716000+00:00',start_station='S12',end_station='S14',start_orig_no=1022,end_orig_no=1022,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.24.csv.gz',start_timestamp ='1976-12-02 22:34:39.517000+00:00',end_timestamp='1976-12-03 01:38:55.956000+00:00',start_station='S12',end_station='S14',start_orig_no=2516,end_orig_no=3736,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.25.csv.gz',start_timestamp ='1976-12-03 08:38:15.812000+00:00',end_timestamp='1976-12-03 11:00:55.926000+00:00',start_station='S12',end_station='S14',start_orig_no=2862,end_orig_no=3805,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.25.csv.gz',start_timestamp ='1976-12-03 08:38:16.416000+00:00',end_timestamp='1976-12-03 11:00:55.926000+00:00',start_station='S12',end_station='S14',start_orig_no=33654,end_orig_no=3805,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.32.csv.gz',start_timestamp ='1976-12-07 07:07:14.942000+00:00',end_timestamp='1976-12-07 09:39:52.893000+00:00',start_station='S12',end_station='S14',start_orig_no=2936,end_orig_no=3945,orig_ground_station1=6)
    # doesn't help
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.33.csv.gz',start_timestamp ='1976-12-07 18:59:06.968000+00:00',end_timestamp='1976-12-07 18:59:53.111000+00:00',start_station='S12',end_station='S14',start_orig_no=3787,end_orig_no=3790,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.37.csv.gz',start_timestamp ='1976-12-10 08:15:05.217000+00:00',end_timestamp='1976-12-10 08:45:52.308000+00:00',start_station='S12',end_station='S14',start_orig_no=2317,end_orig_no=2520,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.38.csv.gz',start_timestamp ='1976-12-11 06:03:35.675000+00:00',end_timestamp='1976-12-11 07:38:54.421000+00:00',start_station='S12',end_station='S14',start_orig_no=2421,end_orig_no=3059,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.39.csv.gz',start_timestamp ='1976-12-11 20:30:53.329000+00:00',end_timestamp='1976-12-11 21:52:04.661000+00:00',start_station='S12',end_station='S14',start_orig_no=2899,end_orig_no=3436,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.42.csv.gz',start_timestamp ='1976-12-13 06:04:46.012000+00:00',end_timestamp='1976-12-13 09:20:57.238000+00:00',start_station='S12',end_station='S14',start_orig_no=2246,end_orig_no=3549,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.43.csv.gz',start_timestamp ='1976-12-14 00:23:38.813000+00:00',end_timestamp='1976-12-14 00:23:49.720000+00:00',start_station='S12',end_station='S14',start_orig_no=3212,end_orig_no=17430,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.46.csv.gz',start_timestamp ='1976-12-15 17:15:46.831000+00:00',end_timestamp='1976-12-15 17:15:48.832000+00:00',start_station='S12',end_station='S14',start_orig_no=2869,end_orig_no=2869,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.46.csv.gz',start_timestamp ='1976-12-16 00:00:51.187000+00:00',end_timestamp='1976-12-16 02:24:50.814000+00:00',start_station='S12',end_station='S14',start_orig_no=2724,end_orig_no=3677,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.5.csv.gz',start_timestamp ='1976-11-23 06:30:37.083000+00:00',end_timestamp='1976-11-23 07:01:24.527000+00:00',start_station='S12',end_station='S16',start_orig_no=1997,end_orig_no=2200,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.50.csv.gz',start_timestamp ='1976-12-18 12:02:29.092000+00:00',end_timestamp='1976-12-18 13:05:52.097000+00:00',start_station='S12',end_station='S14',start_orig_no=3040,end_orig_no=3459,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.54.csv.gz',start_timestamp ='1976-12-20 17:35:55.614000+00:00',end_timestamp='1976-12-20 19:59:54.463000+00:00',start_station='S12',end_station='S14',start_orig_no=2399,end_orig_no=3352,orig_ground_station1=3)
    # not helping
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.55.csv.gz',start_timestamp ='1976-12-21 03:52:57.411000+00:00',end_timestamp='1976-12-21 07:29:58.902000+00:00',start_station='S12',end_station='S14',start_orig_no=16,end_orig_no=1776,orig_ground_station1=9)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.55.csv.gz',start_timestamp ='1976-12-21 05:26:16.302000+00:00',end_timestamp='1976-12-21 07:29:58.902000+00:00',start_station='S12',end_station='S14',start_orig_no=958,end_orig_no=1776,orig_ground_station1=9)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.55.csv.gz',start_timestamp ='1976-12-21 05:26:19.924000+00:00',end_timestamp='1976-12-21 07:29:58.902000+00:00',start_station='S12',end_station='S14',start_orig_no=1152,end_orig_no=1776,orig_ground_station1=9)
    
    # Added manually
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.55.csv.gz',start_timestamp ='1976-12-21T05:26:15.847000Z',end_timestamp='1976-12-21 07:29:58.902000+00:00',start_station='S14',end_station='S14',start_orig_no=958,end_orig_no=1776,orig_ground_station1=9)
    # 958,"1976-12-21T05:26:15.847000Z",20,"S14",0,9,511,490,509,511,490,508,511,489,508,511,489,508,"00011",22,"1

    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.2.csv.gz',start_timestamp ='1976-12-22 02:13:32.893000+00:00',end_timestamp='1976-12-22 04:00:51.428000+00:00',start_station='S12',end_station='S14',start_orig_no=2400,end_orig_no=3110,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.24.csv.gz',start_timestamp ='1977-01-01 16:54:59.129000+00:00',end_timestamp='1977-01-01 17:04:56.422000+00:00',start_station='S12',end_station='S14',start_orig_no=3359,end_orig_no=3424,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.25.csv.gz',start_timestamp ='1977-01-02 16:07:59.373000+00:00',end_timestamp='1977-01-02 17:59:57.120000+00:00',start_station='S12',end_station='S14',start_orig_no=2485,end_orig_no=3226,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.28.csv.gz',start_timestamp ='1977-01-04 00:31:30.972000+00:00',end_timestamp='1977-01-04 03:58:38.948000+00:00',start_station='S12',end_station='S14',start_orig_no=2498,end_orig_no=3905,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.28.csv.gz',start_timestamp ='1977-01-04 00:49:59.466000+00:00',end_timestamp='1977-01-04 03:58:38.948000+00:00',start_station='S12',end_station='S14',start_orig_no=2656,end_orig_no=3905,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.29.csv.gz',start_timestamp ='1977-01-04 09:10:52.560000+00:00',end_timestamp='1977-01-04 11:00:53.359000+00:00',start_station='S12',end_station='S14',start_orig_no=2154,end_orig_no=2882,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.31.csv.gz',start_timestamp ='1977-01-05 17:08:59.657000+00:00',end_timestamp='1977-01-05 19:59:50.189000+00:00',start_station='S12',end_station='S14',start_orig_no=955,end_orig_no=2086,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.32.csv.gz',start_timestamp ='1977-01-06 03:47:59.285000+00:00',end_timestamp='1977-01-06 03:59:53.720000+00:00',start_station='S12',end_station='S14',start_orig_no=3181,end_orig_no=3259,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.35.csv.gz',start_timestamp ='1977-01-06 19:38:58.933000+00:00',end_timestamp='1977-01-06 20:59:52.441000+00:00',start_station='S12',end_station='S14',start_orig_no=2293,end_orig_no=2828,orig_ground_station1=5)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.41.csv.gz',start_timestamp ='1977-01-09 22:38:00.189000+00:00',end_timestamp='1977-01-09 22:45:51.485000+00:00',start_station='S12',end_station='S14',start_orig_no=40960,end_orig_no=4090,orig_ground_station1=9)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.44.csv.gz',start_timestamp ='1977-01-11 18:19:34.269000+00:00',end_timestamp='1977-01-11 19:19:57.215000+00:00',start_station='S12',end_station='S14',start_orig_no=1058,end_orig_no=1969,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.46.csv.gz',start_timestamp ='1977-01-12 22:00:28.857000+00:00',end_timestamp='1977-01-13 01:36:54.749000+00:00',start_station='S12',end_station='S14',start_orig_no=2803,end_orig_no=4236,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.50.csv.gz',start_timestamp ='1977-01-15 01:57:27.278000+00:00',end_timestamp='1977-01-15 02:39:51.653000+00:00',start_station='S12',end_station='S14',start_orig_no=2977,end_orig_no=3257,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.54.csv.gz',start_timestamp ='1977-01-17 17:25:10.371000+00:00',end_timestamp='1977-01-17 19:29:50.254000+00:00',start_station='S12',end_station='S14',start_orig_no=2172,end_orig_no=2997,orig_ground_station1=3)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.55.csv.gz',start_timestamp ='1977-01-17 21:55:45.018000+00:00',end_timestamp='1977-01-17 22:29:58.091000+00:00',start_station='S12',end_station='S14',start_orig_no=8192,end_orig_no=1226,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.11.9.csv.gz',start_timestamp ='1976-12-26 18:55:49.897000+00:00',end_timestamp='1976-12-26 22:59:53.131000+00:00',start_station='S12',end_station='S14',start_orig_no=2038,end_orig_no=3654,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.12.csv.gz',start_timestamp ='1977-01-25 17:08:18.053000+00:00',end_timestamp='1977-01-25 18:39:54.559000+00:00',start_station='S12',end_station='S14',start_orig_no=2743,end_orig_no=3349,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.15.csv.gz',start_timestamp ='1977-01-27 12:05:59.097000+00:00',end_timestamp='1977-01-27 12:14:52.711000+00:00',start_station='S12',end_station='S14',start_orig_no=3068,end_orig_no=3126,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.17.csv.gz',start_timestamp ='1977-01-28 19:03:16.027000+00:00',end_timestamp='1977-01-28 19:19:54.554000+00:00',start_station='S12',end_station='S14',start_orig_no=2502,end_orig_no=2610,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.2.csv.gz',start_timestamp ='1977-01-19 06:28:20.175000+00:00',end_timestamp='1977-01-19 07:19:53.925000+00:00',start_station='S12',end_station='S14',start_orig_no=2975,end_orig_no=3316,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.24.csv.gz',start_timestamp ='1977-01-31 23:11:59.012000+00:00',end_timestamp='1977-02-01 01:04:52.295000+00:00',start_station='S12',end_station='S14',start_orig_no=849,end_orig_no=1596,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.29.csv.gz',start_timestamp ='1977-02-03 18:45:59.166000+00:00',end_timestamp='1977-02-03 19:59:56.075000+00:00',start_station='S12',end_station='S14',start_orig_no=2242,end_orig_no=2731,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.29.csv.gz',start_timestamp ='1977-02-04 01:20:58.943000+00:00',end_timestamp='1977-02-04 02:59:53.478000+00:00',start_station='S12',end_station='S14',start_orig_no=2209,end_orig_no=2863,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.35.csv.gz',start_timestamp ='1977-02-07 19:08:58.944000+00:00',end_timestamp='1977-02-07 22:29:49.794000+00:00',start_station='S12',end_station='S14',start_orig_no=2266,end_orig_no=3615,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.35.csv.gz',start_timestamp ='1977-02-07 20:26:59.258000+00:00',end_timestamp='1977-02-07 22:29:49.794000+00:00',start_station='S12',end_station='S14',start_orig_no=2802,end_orig_no=3615,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.37.csv.gz',start_timestamp ='1977-02-08 18:37:59.172000+00:00',end_timestamp='1977-02-08 20:48:50.002000+00:00',start_station='S12',end_station='S14',start_orig_no=2386,end_orig_no=3252,orig_ground_station1=5)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.40.csv.gz',start_timestamp ='1977-02-10 18:49:49.744000+00:00',end_timestamp='1977-02-10 21:08:51.130000+00:00',start_station='S12',end_station='S14',start_orig_no=1557,end_orig_no=2477,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.48.csv.gz',start_timestamp ='1977-02-14 20:35:14.935000+00:00',end_timestamp='1977-02-14 23:14:58.670000+00:00',start_station='S12',end_station='S14',start_orig_no=2047,end_orig_no=1485,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.48.csv.gz',start_timestamp ='1977-02-14 20:35:14.935000+00:00',end_timestamp='1977-02-14 23:14:58.670000+00:00',start_station='S12',end_station='S14',start_orig_no=50305,end_orig_no=1485,orig_ground_station1=9)
    
    # modified
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.49.csv.gz',start_timestamp ='1977-02-15T14:02:59.394000Z',end_timestamp='1977-02-15 14:55:57.946000+00:00',start_station='S12',end_station='S14',start_orig_no=3114,end_orig_no=3464,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.50.csv.gz',start_timestamp ='1977-02-16 04:53:59.374000+00:00',end_timestamp='1977-02-16 06:15:55.989000+00:00',start_station='S12',end_station='S14',start_orig_no=1789,end_orig_no=2331,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.55.csv.gz',start_timestamp ='1977-02-18 17:55:31.696000+00:00',end_timestamp='1977-02-18 20:29:51.346000+00:00',start_station='S12',end_station='S14',start_orig_no=1289,end_orig_no=2310,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.7.csv.gz',start_timestamp ='1977-01-22 19:37:21.569000+00:00',end_timestamp='1977-01-22 20:59:58.005000+00:00',start_station='S12',end_station='S14',start_orig_no=311,end_orig_no=856,orig_ground_station1=2)
    # ???
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.12.8.csv.gz',start_timestamp ='1977-01-23 09:00:41.431000+00:00',end_timestamp='1977-01-23 10:30:55.690000+00:00',start_station='S12',end_station='S14',start_orig_no=2925,end_orig_no=3521,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.13.csv.gz',start_timestamp ='1977-02-25 18:00:32.481000+00:00',end_timestamp='1977-02-25 20:29:57.538000+00:00',start_station='S12',end_station='S14',start_orig_no=2379,end_orig_no=3368,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.16.csv.gz',start_timestamp ='1977-02-28 00:28:23.105000+00:00',end_timestamp='1977-02-28 02:04:57.907000+00:00',start_station='S12',end_station='S14',start_orig_no=2187,end_orig_no=2826,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.22.csv.gz',start_timestamp ='1977-03-03 12:42:47.098000+00:00',end_timestamp='1977-03-03 12:42:58.460000+00:00',start_station='S12',end_station='S14',start_orig_no=2163,end_orig_no=2849,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.34.csv.gz',start_timestamp ='1977-03-09 15:43:34.622000+00:00',end_timestamp='1977-03-09 22:59:50.920000+00:00',start_station='S12',end_station='S14',start_orig_no=127,end_orig_no=3654,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.34.csv.gz',start_timestamp ='1977-03-09 15:43:34.622000+00:00',end_timestamp='1977-03-09 22:59:50.920000+00:00',start_station='S12',end_station='S14',start_orig_no=58228,end_orig_no=3654,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.34.csv.gz',start_timestamp ='1977-03-09 15:43:17.717000+00:00',end_timestamp='1977-03-09 22:59:50.920000+00:00',start_station='S12',end_station='S14',start_orig_no=755,end_orig_no=3654,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.34.csv.gz',start_timestamp ='1977-03-09 19:43:47.472000+00:00',end_timestamp='1977-03-09 22:59:50.920000+00:00',start_station='S12',end_station='S14',start_orig_no=2356,end_orig_no=3654,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.4.csv.gz',start_timestamp ='1977-02-20 12:27:04.412000+00:00',end_timestamp='1977-02-20 14:34:57.353000+00:00',start_station='S12',end_station='S14',start_orig_no=2212,end_orig_no=3057,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.42.csv.gz',start_timestamp ='1977-03-14 17:03:38.297000+00:00',end_timestamp='1977-03-14 17:19:56.187000+00:00',start_station='S12',end_station='S14',start_orig_no=2201,end_orig_no=2308,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.44.csv.gz',start_timestamp ='1977-05-02 21:45:42.946000+00:00',end_timestamp='1977-05-02 21:45:44.678000+00:00',start_station='S12',end_station='S17',start_orig_no=1374,end_orig_no=1374,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.48.csv.gz',start_timestamp ='1977-03-16 17:00:47.988000+00:00',end_timestamp='1977-03-16 19:19:57.080000+00:00',start_station='S12',end_station='S14',start_orig_no=2342,end_orig_no=3263,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.51.csv.gz',start_timestamp ='1977-03-18 14:04:46.222000+00:00',end_timestamp='1977-03-18 17:34:58.357000+00:00',start_station='S12',end_station='S14',start_orig_no=80,end_orig_no=1472,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.57.csv.gz',start_timestamp ='1977-03-22 00:41:59.222000+00:00',end_timestamp='1977-03-22 02:29:51.467000+00:00',start_station='S12',end_station='S14',start_orig_no=828,end_orig_no=1485,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.13.8.csv.gz',start_timestamp ='1977-02-23 06:29:45.954000+00:00',end_timestamp='1977-02-23 06:29:57.036000+00:00',start_station='S12',end_station='S14',start_orig_no=3615,end_orig_no=2857,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.1.csv.gz',start_timestamp ='1977-03-22 15:34:50.513000+00:00',end_timestamp='1977-03-22 18:59:57.094000+00:00',start_station='S12',end_station='S14',start_orig_no=2937,end_orig_no=4295,orig_ground_station1=1)
    
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.16.csv.gz',start_timestamp ='1977-04-02 03:51:20.962000+00:00',end_timestamp='1977-04-02 05:39:54.490000+00:00',start_station='S12',end_station='S14',start_orig_no=2876,end_orig_no=3507,orig_ground_station1=4)
    
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.22.csv.gz',start_timestamp ='1977-04-05 10:19:53.812000+00:00',end_timestamp='1977-04-05 13:54:51.987000+00:00',start_station='S12',end_station='S14',start_orig_no=1164,end_orig_no=3624,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.24.csv.gz',start_timestamp ='1977-04-06 13:54:01.559000+00:00',end_timestamp='1977-04-06 14:49:51.525000+00:00',start_station='S12',end_station='S14',start_orig_no=2918,end_orig_no=3287,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.27.csv.gz',start_timestamp ='1977-04-07 23:49:40.691000+00:00',end_timestamp='1977-04-08 05:31:53.918000+00:00',start_station='S12',end_station='S14',start_orig_no=551,end_orig_no=2816,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.35.csv.gz',start_timestamp ='1977-04-13 02:11:03.948000+00:00',end_timestamp='1977-04-13 03:22:53.593000+00:00',start_station='S12',end_station='S14',start_orig_no=3201,end_orig_no=3676,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.38.csv.gz',start_timestamp ='1977-04-14 17:00:34.569000+00:00',end_timestamp='1977-04-14 18:44:51.746000+00:00',start_station='S12',end_station='S14',start_orig_no=3176,end_orig_no=3866,orig_ground_station1=8)
    # new one
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.39.csv.gz',start_timestamp ='1977-04-15 09:41:35.159000+00:00',end_timestamp='1977-04-15 12:04:57.351000+00:00',start_station='S15',end_station='S14',start_orig_no=2088,end_orig_no=3037,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.40.csv.gz',start_timestamp ='1977-04-15 16:22:22.073000+00:00',end_timestamp='1977-04-15 19:49:53.458000+00:00',start_station='S12',end_station='S14',start_orig_no=1738,end_orig_no=3139,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.40.csv.gz',start_timestamp ='1977-04-15 18:00:37.412000+00:00',end_timestamp='1977-04-15 19:49:53.458000+00:00',start_station='S12',end_station='S14',start_orig_no=2416,end_orig_no=3139,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.42.csv.gz',start_timestamp ='1977-04-17 04:08:29.670000+00:00',end_timestamp='1977-04-17 04:19:57.076000+00:00',start_station='S12',end_station='S14',start_orig_no=3145,end_orig_no=3220,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.43.csv.gz',start_timestamp ='1977-04-17 16:41:45.862000+00:00',end_timestamp='1977-04-17 19:44:50.430000+00:00',start_station='S12',end_station='S14',start_orig_no=2467,end_orig_no=3679,orig_ground_station1=8)
    # not required 
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.44.csv.gz',start_timestamp ='1977-04-17 22:04:00.216000+00:00',end_timestamp='1977-04-17 23:59:56.013000+00:00',start_station='S12',end_station='S14',start_orig_no=8192,end_orig_no=1750,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.44.csv.gz',start_timestamp ='1977-04-18 15:30:50.994000+00:00',end_timestamp='1977-04-18 15:32:47.930000+00:00',start_station='S12',end_station='S14',start_orig_no=1780,end_orig_no=1792,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.49.csv.gz',start_timestamp ='1977-04-21 16:00:38.313000+00:00',end_timestamp='1977-04-21 16:40:02.055000+00:00',start_station='S12',end_station='S14',start_orig_no=1389,end_orig_no=1649,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.51.csv.gz',start_timestamp ='1977-04-21 19:00:44.228000+00:00',end_timestamp='1977-04-21 20:29:55.854000+00:00',start_station='S12',end_station='S14',start_orig_no=2610,end_orig_no=3200,orig_ground_station1=3)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.55.csv.gz',start_timestamp ='1977-04-23 15:58:44.349000+00:00',end_timestamp='1977-04-23 22:28:51.916000+00:00',start_station='S12',end_station='S14',start_orig_no=32768,end_orig_no=3246,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.14.8.csv.gz',start_timestamp ='1977-03-28 09:22:35.492000+00:00',end_timestamp='1977-03-28 10:59:50.255000+00:00',start_station='S12',end_station='S14',start_orig_no=325,end_orig_no=963,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.14.csv.gz',start_timestamp ='1977-04-29 19:27:48.065000+00:00',end_timestamp='1977-04-29 22:39:53.643000+00:00',start_station='S12',end_station='S17',start_orig_no=1515,end_orig_no=3105,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.21.csv.gz',start_timestamp ='1977-05-03 13:20:31.827000+00:00',end_timestamp='1977-05-03 16:39:52.725000+00:00',start_station='S12',end_station='S17',start_orig_no=1652,end_orig_no=3302,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.23.csv.gz',start_timestamp ='1977-05-04 08:50:00.305000+00:00',end_timestamp='1977-05-04 12:29:52.261000+00:00',start_station='S12',end_station='S17',start_orig_no=2195,end_orig_no=4015,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.25.csv.gz',start_timestamp ='1977-05-05 04:58:59.139000+00:00',end_timestamp='1977-05-05 05:44:51.864000+00:00',start_station='S12',end_station='S17',start_orig_no=4006,end_orig_no=4385,orig_ground_station1=7)
    # not requured 
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.3.csv.gz',start_timestamp ='1977-04-25 07:42:44.433000+00:00',end_timestamp='1977-04-25 08:19:56.448000+00:00',start_station='S12',end_station='S17',start_orig_no=24576,end_orig_no=1044,orig_ground_station1=6)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.3.csv.gz',start_timestamp ='1977-04-25 10:34:32.235000+00:00',end_timestamp='1977-04-25 10:34:33.686000+00:00',start_station='S12',end_station='S17',start_orig_no=1219,end_orig_no=1219,orig_ground_station1=5)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.30.csv.gz',start_timestamp ='1977-05-07 09:46:13.627000+00:00',end_timestamp='1977-05-07 09:51:19.762000+00:00',start_station='S12',end_station='S17',start_orig_no=62098,end_orig_no=2200,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.38.csv.gz',start_timestamp ='1977-05-10 19:26:40.167000+00:00',end_timestamp='1977-05-10 20:21:57.843000+00:00',start_station='S12',end_station='S17',start_orig_no=2198,end_orig_no=2655,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.43.csv.gz',start_timestamp ='1977-05-12 23:35:35.513000+00:00',end_timestamp='1977-05-13 04:29:56.357000+00:00',start_station='S12',end_station='S17',start_orig_no=2873,end_orig_no=3900,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.46.csv.gz',start_timestamp ='1977-05-14 09:49:18.424000+00:00',end_timestamp='1977-05-14 10:32:56.167000+00:00',start_station='S12',end_station='S17',start_orig_no=15,end_orig_no=4121,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.15.51.csv.gz',start_timestamp ='1977-05-16 11:40:33.192000+00:00',end_timestamp='1977-05-16 13:24:53.380000+00:00',start_station='S12',end_station='S17',start_orig_no=64640,end_orig_no=1614,orig_ground_station1=11)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.11.csv.gz',start_timestamp ='1977-05-23 10:17:23.821000+00:00',end_timestamp='1977-05-23 10:50:03.372000+00:00',start_station='S12',end_station='S17',start_orig_no=3043,end_orig_no=3322,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.11.csv.gz',start_timestamp ='1977-05-23 10:41:43.775000+00:00',end_timestamp='1977-05-23 10:50:03.372000+00:00',start_station='S12',end_station='S17',start_orig_no=3254,end_orig_no=3322,orig_ground_station1=5)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.13.csv.gz',start_timestamp ='1977-05-23 23:12:16.484000+00:00',end_timestamp='1977-05-23 23:13:04.809000+00:00',start_station='S12',end_station='S17',start_orig_no=57344,end_orig_no=2200,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.14.csv.gz',start_timestamp ='1977-05-24 18:08:07.725000+00:00',end_timestamp='1977-05-24 19:29:58.757000+00:00',start_station='S12',end_station='S17',start_orig_no=2064,end_orig_no=2741,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.18.csv.gz',start_timestamp ='1977-05-26 05:38:06.539000+00:00',end_timestamp='1977-05-26 06:19:55.241000+00:00',start_station='S12',end_station='S17',start_orig_no=51198,end_orig_no=1526,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.18.csv.gz',start_timestamp ='1977-05-26 07:07:39.616000+00:00',end_timestamp='1977-05-26 09:32:54.327000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=3403,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.19.csv.gz',start_timestamp ='1977-05-26 22:56:44.431000+00:00',end_timestamp='1977-05-26 22:56:53.738000+00:00',start_station='S12',end_station='S17',start_orig_no=1867,end_orig_no=2858,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.21.csv.gz',start_timestamp ='1977-05-27 23:37:30.429000+00:00',end_timestamp='1977-05-28 01:49:57.546000+00:00',start_station='S12',end_station='S17',start_orig_no=3057,end_orig_no=4153,orig_ground_station1=1)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.26.csv.gz',start_timestamp ='1977-05-30 10:09:47.434000+00:00',end_timestamp='1977-05-30 10:09:49.254000+00:00',start_station='S12',end_station='S17',start_orig_no=0,end_orig_no=0,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.31.csv.gz',start_timestamp ='1977-06-01 02:19:12.026000+00:00',end_timestamp='1977-06-01 03:59:55.976000+00:00',start_station='S12',end_station='S17',start_orig_no=824,end_orig_no=1656,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.31.csv.gz',start_timestamp ='1977-06-01 02:19:12.026000+00:00',end_timestamp='1977-06-01 03:59:55.976000+00:00',start_station='S12',end_station='S17',start_orig_no=33044,end_orig_no=1656,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.33.csv.gz',start_timestamp ='1977-06-02 08:15:27.808000+00:00',end_timestamp='1977-06-02 08:26:42.174000+00:00',start_station='S12',end_station='S17',start_orig_no=2858,end_orig_no=1848,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.33.csv.gz',start_timestamp ='1977-06-02 08:22:11.185000+00:00',end_timestamp='1977-06-02 08:22:12.964000+00:00',start_station='S12',end_station='S17',start_orig_no=4,end_orig_no=4,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.35.csv.gz',start_timestamp ='1977-06-03 02:37:09.490000+00:00',end_timestamp='1977-06-03 02:59:52.594000+00:00',start_station='S12',end_station='S17',start_orig_no=656,end_orig_no=1200,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.39.csv.gz',start_timestamp ='1977-06-04 22:44:47.322000+00:00',end_timestamp='1977-06-04 22:44:56.407000+00:00',start_station='S12',end_station='S17',start_orig_no=2154,end_orig_no=2870,orig_ground_station1=5)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.40.csv.gz',start_timestamp ='1977-06-05 02:36:49.306000+00:00',end_timestamp='1977-06-05 06:17:23.124000+00:00',start_station='S12',end_station='S17',start_orig_no=8,end_orig_no=1834,orig_ground_station1=1)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.40.csv.gz',start_timestamp ='1977-06-05 02:36:49.306000+00:00',end_timestamp='1977-06-05 06:17:23.124000+00:00',start_station='S12',end_station='S17',start_orig_no=9222,end_orig_no=1834,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.42.csv.gz',start_timestamp ='1977-06-05 14:18:27.983000+00:00',end_timestamp='1977-06-05 16:30:56.873000+00:00',start_station='S12',end_station='S17',start_orig_no=2416,end_orig_no=3511,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.43.csv.gz',start_timestamp ='1977-06-05 22:32:59.566000+00:00',end_timestamp='1977-06-05 23:34:55.570000+00:00',start_station='S12',end_station='S17',start_orig_no=3065,end_orig_no=3577,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.46.csv.gz',start_timestamp ='1977-06-07 07:45:06.939000+00:00',end_timestamp='1977-06-07 08:49:56.648000+00:00',start_station='S12',end_station='S17',start_orig_no=1657,end_orig_no=2193,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.47.csv.gz',start_timestamp ='1977-06-07 22:42:58.935000+00:00',end_timestamp='1977-06-07 23:39:58.037000+00:00',start_station='S12',end_station='S17',start_orig_no=876,end_orig_no=1347,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.5.csv.gz',start_timestamp ='1977-05-20 04:39:58.847000+00:00',end_timestamp='1977-05-20 09:14:53.218000+00:00',start_station='S12',end_station='S17',start_orig_no=1839,end_orig_no=3769,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.16.53.csv.gz',start_timestamp ='1977-06-10 12:11:38.985000+00:00',end_timestamp='1977-06-10 13:49:57.828000+00:00',start_station='S12',end_station='S17',start_orig_no=1631,end_orig_no=2443,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.11.csv.gz',start_timestamp ='1977-06-16 17:39:59.177000+00:00',end_timestamp='1977-06-16 19:05:56.994000+00:00',start_station='S12',end_station='S14',start_orig_no=3371,end_orig_no=4082,orig_ground_station1=7)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz',start_timestamp ='1977-06-16 23:13:48.560000+00:00',end_timestamp='1977-06-17 02:43:55.670000+00:00',start_station='S12',end_station='S17',start_orig_no=2205,end_orig_no=3947,orig_ground_station1=4)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz',start_timestamp ='1977-06-16 23:13:48.560000+00:00',end_timestamp='1977-06-17 02:43:55.670000+00:00',start_station='S12',end_station='S17',start_orig_no=62976,end_orig_no=3947,orig_ground_station1=4)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz',start_timestamp ='1977-06-16 23:14:06.070000+00:00',end_timestamp='1977-06-17 02:43:55.670000+00:00',start_station='S12',end_station='S17',start_orig_no=2208,end_orig_no=3947,orig_ground_station1=4)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz',start_timestamp ='1977-06-16 23:14:06.070000+00:00',end_timestamp='1977-06-17 02:43:55.670000+00:00',start_station='S12',end_station='S17',start_orig_no=40960,end_orig_no=3947,orig_ground_station1=4)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz',start_timestamp ='1977-06-16 23:14:10.900000+00:00',end_timestamp='1977-06-17 02:43:55.670000+00:00',start_station='S12',end_station='S17',start_orig_no=0,end_orig_no=3947,orig_ground_station1=4)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz',start_timestamp ='1977-06-16 23:14:24.787000+00:00',end_timestamp='1977-06-17 02:43:55.670000+00:00',start_station='S12',end_station='S17',start_orig_no=2210,end_orig_no=3947,orig_ground_station1=4)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz',start_timestamp ='1977-06-16 23:14:24.787000+00:00',end_timestamp='1977-06-17 02:43:55.670000+00:00',start_station='S12',end_station='S17',start_orig_no=0,end_orig_no=3947,orig_ground_station1=4)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.12.csv.gz',start_timestamp ='1977-06-16 23:14:42.901000+00:00',end_timestamp='1977-06-17 02:43:55.670000+00:00',start_station='S12',end_station='S17',start_orig_no=2213,end_orig_no=3947,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.15.csv.gz',start_timestamp ='1977-06-18 06:20:53.823000+00:00',end_timestamp='1977-06-18 07:53:59.004000+00:00',start_station='S12',end_station='S17',start_orig_no=3479,end_orig_no=4249,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.21.csv.gz',start_timestamp ='1977-06-21 00:10:32.113000+00:00',end_timestamp='1977-06-21 02:49:56.720000+00:00',start_station='S12',end_station='S17',start_orig_no=0,end_orig_no=3144,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.21.csv.gz',start_timestamp ='1977-06-21 00:10:32.113000+00:00',end_timestamp='1977-06-21 02:49:56.720000+00:00',start_station='S12',end_station='S17',start_orig_no=64671,end_orig_no=3144,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.22.csv.gz',start_timestamp ='1977-06-22 02:19:19.444000+00:00',end_timestamp='1977-06-22 02:49:58.900000+00:00',start_station='S12',end_station='S17',start_orig_no=3148,end_orig_no=3401,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.27.csv.gz',start_timestamp ='1977-06-24 02:48:54.149000+00:00',end_timestamp='1977-06-24 03:40:58.357000+00:00',start_station='S12',end_station='S17',start_orig_no=2866,end_orig_no=4025,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.30.csv.gz',start_timestamp ='1977-06-25 03:00:51.564000+00:00',end_timestamp='1977-06-25 05:09:55.834000+00:00',start_station='S12',end_station='S17',start_orig_no=3144,end_orig_no=4212,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.32.csv.gz',start_timestamp ='1977-06-25 23:44:23.924000+00:00',end_timestamp='1977-06-25 23:49:59.079000+00:00',start_station='S12',end_station='S17',start_orig_no=2158,end_orig_no=2203,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.33.csv.gz',start_timestamp ='1977-06-26 08:53:08.854000+00:00',end_timestamp='1977-06-26 10:16:58.222000+00:00',start_station='S12',end_station='S17',start_orig_no=3691,end_orig_no=4383,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.34.csv.gz',start_timestamp ='1977-06-27 00:54:23.501000+00:00',end_timestamp='1977-06-27 02:49:55.925000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=3157,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.36.csv.gz',start_timestamp ='1977-06-28 01:16:43.780000+00:00',end_timestamp='1977-06-28 01:49:55.915000+00:00',start_station='S12',end_station='S17',start_orig_no=1823,end_orig_no=2097,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.38.csv.gz',start_timestamp ='1977-06-29 02:28:39.598000+00:00',end_timestamp='1977-06-29 02:39:52.771000+00:00',start_station='S12',end_station='S17',start_orig_no=3752,end_orig_no=3844,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.39.csv.gz',start_timestamp ='1977-06-29 05:50:48.460000+00:00',end_timestamp='1977-06-29 05:50:57.781000+00:00',start_station='S12',end_station='S17',start_orig_no=1736,end_orig_no=2868,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.40.csv.gz',start_timestamp ='1977-06-29 16:09:27.424000+00:00',end_timestamp='1977-06-29 18:27:56.529000+00:00',start_station='S12',end_station='S17',start_orig_no=787,end_orig_no=1933,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.41.csv.gz',start_timestamp ='1977-06-29 22:57:43.036000+00:00',end_timestamp='1977-06-30 02:43:28.903000+00:00',start_station='S12',end_station='S14',start_orig_no=2200,end_orig_no=4067,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.5.csv.gz',start_timestamp ='1977-06-14 04:03:02.239000+00:00',end_timestamp='1977-06-14 04:31:33.687000+00:00',start_station='S12',end_station='S17',start_orig_no=1966,end_orig_no=2200,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.17.7.csv.gz',start_timestamp ='1977-06-15 01:06:41.980000+00:00',end_timestamp='1977-06-15 03:02:50.682000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=3162,orig_ground_station1=12)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.19.csv.gz',start_timestamp ='1977-07-13 19:02:22.356000+00:00',end_timestamp='1977-07-13 19:04:53.856000+00:00',start_station='S12',end_station='S17',start_orig_no=7,end_orig_no=2193,orig_ground_station1=8)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.19.csv.gz',start_timestamp ='1977-07-13 19:02:27.186000+00:00',end_timestamp='1977-07-13 19:04:53.856000+00:00',start_station='S12',end_station='S17',start_orig_no=61462,end_orig_no=2193,orig_ground_station1=8)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.19.csv.gz',start_timestamp ='1977-07-14 01:15:45.244000+00:00',end_timestamp='1977-07-14 01:31:56.825000+00:00',start_station='S12',end_station='S17',start_orig_no=4116,end_orig_no=2327,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.21.csv.gz',start_timestamp ='1977-07-15 01:51:03.995000+00:00',end_timestamp='1977-07-15 03:51:57.816000+00:00',start_station='S12',end_station='S17',start_orig_no=0,end_orig_no=1239,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.22.csv.gz',start_timestamp ='1977-07-15 13:20:20.379000+00:00',end_timestamp='1977-07-15 13:37:43.213000+00:00',start_station='S12',end_station='S17',start_orig_no=2057,end_orig_no=2200,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.25.csv.gz',start_timestamp ='1977-07-17 05:03:49.167000+00:00',end_timestamp='1977-07-17 05:17:05.568000+00:00',start_station='S12',end_station='S17',start_orig_no=2091,end_orig_no=2200,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.25.csv.gz',start_timestamp ='1977-07-17 05:08:17.249000+00:00',end_timestamp='1977-07-17 05:17:05.568000+00:00',start_station='S12',end_station='S17',start_orig_no=2128,end_orig_no=2200,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.25.csv.gz',start_timestamp ='1977-07-17 05:08:17.249000+00:00',end_timestamp='1977-07-17 05:17:05.568000+00:00',start_station='S12',end_station='S17',start_orig_no=2128,end_orig_no=2200,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.25.csv.gz',start_timestamp ='1977-07-17 05:08:17.249000+00:00',end_timestamp='1977-07-17 05:17:05.568000+00:00',start_station='S12',end_station='S17',start_orig_no=2128,end_orig_no=2200,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.25.csv.gz',start_timestamp ='1977-07-17 05:08:17.249000+00:00',end_timestamp='1977-07-17 05:17:05.568000+00:00',start_station='S12',end_station='S17',start_orig_no=2128,end_orig_no=2200,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.28.csv.gz',start_timestamp ='1977-07-18 07:28:57.146000+00:00',end_timestamp='1977-07-18 08:19:53.720000+00:00',start_station='S12',end_station='S17',start_orig_no=2194,end_orig_no=2615,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.29.csv.gz',start_timestamp ='1977-07-18 19:20:11.032000+00:00',end_timestamp='1977-07-18 21:54:58.751000+00:00',start_station='S12',end_station='S17',start_orig_no=2047,end_orig_no=3328,orig_ground_station1=6)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.32.csv.gz',start_timestamp ='1977-07-20 04:04:25.524000+00:00',end_timestamp='1977-07-20 06:19:02.476000+00:00',start_station='S12',end_station='S17',start_orig_no=1111,end_orig_no=2225,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.33.csv.gz',start_timestamp ='1977-07-20 06:59:07.884000+00:00',end_timestamp='1977-07-20 09:59:53.175000+00:00',start_station='S12',end_station='S17',start_orig_no=2560,end_orig_no=4056,orig_ground_station1=9)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.34.csv.gz',start_timestamp ='1977-07-20 20:15:26.204000+00:00',end_timestamp='1977-07-20 22:29:58.524000+00:00',start_station='S12',end_station='S17',start_orig_no=0,end_orig_no=1118,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.37.csv.gz',start_timestamp ='1977-07-22 09:17:27.393000+00:00',end_timestamp='1977-07-22 12:19:55.831000+00:00',start_station='S12',end_station='S17',start_orig_no=1968,end_orig_no=3478,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.39.csv.gz',start_timestamp ='1977-07-23 08:29:01.611000+00:00',end_timestamp='1977-07-23 10:59:57.244000+00:00',start_station='S12',end_station='S17',start_orig_no=2929,end_orig_no=4178,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.39.csv.gz',start_timestamp ='1977-07-23 14:51:30.395000+00:00',end_timestamp='1977-07-23 16:04:54.495000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=2808,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.4.csv.gz',start_timestamp ='1977-07-07 19:02:14.513000+00:00',end_timestamp='1977-07-07 21:59:52.450000+00:00',start_station='S12',end_station='S17',start_orig_no=2848,end_orig_no=3519,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.4.csv.gz',start_timestamp ='1977-07-08 01:00:59.080000+00:00',end_timestamp='1977-07-08 02:00:55.125000+00:00',start_station='S12',end_station='S14',start_orig_no=1558,end_orig_no=2039,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.41.csv.gz',start_timestamp ='1977-07-24 09:53:49.397000+00:00',end_timestamp='1977-07-24 10:44:53.520000+00:00',start_station='S12',end_station='S17',start_orig_no=2973,end_orig_no=3395,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.44.csv.gz',start_timestamp ='1977-07-25 20:42:53.384000+00:00',end_timestamp='1977-07-25 22:49:39.632000+00:00',start_station='S12',end_station='S17',start_orig_no=1809,end_orig_no=2858,orig_ground_station1=1)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.46.csv.gz',start_timestamp ='1977-07-26 16:28:04.831000+00:00',end_timestamp='1977-07-26 20:14:58.722000+00:00',start_station='S12',end_station='S17',start_orig_no=1408,end_orig_no=1976,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.18.52.csv.gz',start_timestamp ='1977-07-29 04:44:04.705000+00:00',end_timestamp='1977-07-29 08:59:56.984000+00:00',start_station='S12',end_station='S17',start_orig_no=34639,end_orig_no=2727,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.16.csv.gz',start_timestamp ='1977-08-06 01:14:04.782000+00:00',end_timestamp='1977-08-06 01:15:32.173000+00:00',start_station='S12',end_station='S16',start_orig_no=42326,end_orig_no=692,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.17.csv.gz',start_timestamp ='1977-08-06 20:49:44.187000+00:00',end_timestamp='1977-08-06 20:49:45.962000+00:00',start_station='S12',end_station='S17',start_orig_no=2860,end_orig_no=2860,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.2.csv.gz',start_timestamp ='1977-07-31 19:19:59.086000+00:00',end_timestamp='1977-07-31 21:00:57.215000+00:00',start_station='S12',end_station='S17',start_orig_no=3194,end_orig_no=3862,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.20.csv.gz',start_timestamp ='1977-08-07 23:49:01.832000+00:00',end_timestamp='1977-08-08 00:34:54.195000+00:00',start_station='S12',end_station='S17',start_orig_no=2974,end_orig_no=3353,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.22.csv.gz',start_timestamp ='1977-08-09 06:32:45.236000+00:00',end_timestamp='1977-08-09 06:32:54.382000+00:00',start_station='S12',end_station='S14',start_orig_no=1504,end_orig_no=2855,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.28.csv.gz',start_timestamp ='1977-08-11 07:12:41.514000+00:00',end_timestamp='1977-08-11 08:34:56.753000+00:00',start_station='S12',end_station='S17',start_orig_no=8402,end_orig_no=1493,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.28.csv.gz',start_timestamp ='1977-08-11 14:30:36.519000+00:00',end_timestamp='1977-08-11 17:04:55.237000+00:00',start_station='S12',end_station='S17',start_orig_no=829,end_orig_no=2106,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.31.csv.gz',start_timestamp ='1977-08-13 03:37:58.990000+00:00',end_timestamp='1977-08-13 04:35:56.205000+00:00',start_station='S12',end_station='S17',start_orig_no=3081,end_orig_no=3560,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.31.csv.gz',start_timestamp ='1977-08-13 10:39:19.449000+00:00',end_timestamp='1977-08-13 10:39:57.178000+00:00',start_station='S12',end_station='S17',start_orig_no=15,end_orig_no=1170,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.39.csv.gz',start_timestamp ='1977-08-17 08:29:34.797000+00:00',end_timestamp='1977-08-17 09:10:58.976000+00:00',start_station='S12',end_station='S17',start_orig_no=3780,end_orig_no=4122,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.40.csv.gz',start_timestamp ='1977-08-18 07:33:45.002000+00:00',end_timestamp='1977-08-18 08:51:15.742000+00:00',start_station='S12',end_station='S17',start_orig_no=2975,end_orig_no=3616,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.46.csv.gz',start_timestamp ='1977-08-21 09:23:43.708000+00:00',end_timestamp='1977-08-21 10:19:58.829000+00:00',start_station='S12',end_station='S17',start_orig_no=2982,end_orig_no=3447,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.49.csv.gz',start_timestamp ='1977-08-23 09:22:44.044000+00:00',end_timestamp='1977-08-23 10:09:54.153000+00:00',start_station='S12',end_station='S17',start_orig_no=1756,end_orig_no=2101,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.5.csv.gz',start_timestamp ='1977-08-02 13:44:09.111000+00:00',end_timestamp='1977-08-02 16:02:35.732000+00:00',start_station='S12',end_station='S17',start_orig_no=1068,end_orig_no=2934,orig_ground_station1=6)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.50.csv.gz',start_timestamp ='1977-08-23 14:55:35.581000+00:00',end_timestamp='1977-08-23 15:09:51.964000+00:00',start_station='S12',end_station='S17',start_orig_no=57344,end_orig_no=1650,orig_ground_station1=9)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.50.csv.gz',start_timestamp ='1977-08-23 19:04:10.039000+00:00',end_timestamp='1977-08-23 21:43:47.757000+00:00',start_station='S12',end_station='S17',start_orig_no=39987,end_orig_no=3492,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.6.csv.gz',start_timestamp ='1977-08-02 16:41:26.691000+00:00',end_timestamp='1977-08-02 19:01:56.757000+00:00',start_station='S12',end_station='S17',start_orig_no=3195,end_orig_no=4125,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.19.8.csv.gz',start_timestamp ='1977-08-04 04:19:28.477000+00:00',end_timestamp='1977-08-04 08:04:57.552000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=3694,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.10.csv.gz',start_timestamp ='1976-04-01 06:09:59.692000+00:00',end_timestamp='1976-04-01 08:08:51.884000+00:00',start_station='S12',end_station='S16',start_orig_no=25558,end_orig_no=2878,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.10.csv.gz',start_timestamp ='1976-04-01 06:09:59.692000+00:00',end_timestamp='1976-04-01 08:08:51.884000+00:00',start_station='S12',end_station='S16',start_orig_no=57468,end_orig_no=2878,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.10.csv.gz',start_timestamp ='1976-04-01 06:16:59.296000+00:00',end_timestamp='1976-04-01 08:08:51.884000+00:00',start_station='S12',end_station='S16',start_orig_no=2136,end_orig_no=2878,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.10.csv.gz',start_timestamp ='1976-04-01 19:52:38.662000+00:00',end_timestamp='1976-04-01 19:52:41.397000+00:00',start_station='S12',end_station='S17',start_orig_no=32768,end_orig_no=32768,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.15.csv.gz',start_timestamp ='1976-04-03 17:24:37.945000+00:00',end_timestamp='1976-04-03 21:23:58.687000+00:00',start_station='S12',end_station='S17',start_orig_no=33578,end_orig_no=1801,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.15.csv.gz',start_timestamp ='1976-04-03 22:20:29.136000+00:00',end_timestamp='1976-04-03 23:48:58.854000+00:00',start_station='S12',end_station='S17',start_orig_no=2199,end_orig_no=2783,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.18.csv.gz',start_timestamp ='1976-04-05 09:27:12.149000+00:00',end_timestamp='1976-04-05 11:00:58.255000+00:00',start_station='S12',end_station='S17',start_orig_no=2875,end_orig_no=2821,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.21.csv.gz',start_timestamp ='1976-04-06 07:57:03.079000+00:00',end_timestamp='1976-04-06 07:57:04.530000+00:00',start_station='S12',end_station='S17',start_orig_no=2660,end_orig_no=2861,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.21.csv.gz',start_timestamp ='1976-04-06 07:57:03.079000+00:00',end_timestamp='1976-04-06 07:57:04.530000+00:00',start_station='S12',end_station='S17',start_orig_no=2861,end_orig_no=2861,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.21.csv.gz',start_timestamp ='1976-04-06 16:57:26.134000+00:00',end_timestamp='1976-04-06 19:53:53.510000+00:00',start_station='S12',end_station='S17',start_orig_no=2165,end_orig_no=3333,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.27.csv.gz',start_timestamp ='1976-04-09 00:20:32.178000+00:00',end_timestamp='1976-04-09 00:52:50.680000+00:00',start_station='S12',end_station='S17',start_orig_no=1439,end_orig_no=1651,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.29.csv.gz',start_timestamp ='1976-04-09 22:22:42.622000+00:00',end_timestamp='1976-04-09 22:22:45.189000+00:00',start_station='S12',end_station='S17',start_orig_no=2864,end_orig_no=2864,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.29.csv.gz',start_timestamp ='1976-04-10 05:19:38.106000+00:00',end_timestamp='1976-04-10 06:54:54.803000+00:00',start_station='S12',end_station='S17',start_orig_no=629,end_orig_no=1258,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.31.csv.gz',start_timestamp ='1976-04-10 14:35:07.540000+00:00',end_timestamp='1976-04-10 15:59:56.376000+00:00',start_station='S12',end_station='S17',start_orig_no=2293,end_orig_no=2854,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.40.csv.gz',start_timestamp ='1976-04-14 18:30:46.761000+00:00',end_timestamp='1976-04-14 19:47:50.697000+00:00',start_station='S12',end_station='S17',start_orig_no=2259,end_orig_no=2769,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.45.csv.gz',start_timestamp ='1976-04-18 06:00:32.519000+00:00',end_timestamp='1976-04-18 07:06:55.263000+00:00',start_station='S12',end_station='S15',start_orig_no=3048,end_orig_no=3486,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.47.csv.gz',start_timestamp ='1976-04-19 11:53:51.870000+00:00',end_timestamp='1976-04-19 12:24:50.876000+00:00',start_station='S12',end_station='S17',start_orig_no=2865,end_orig_no=1782,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.49.csv.gz',start_timestamp ='1976-04-20 15:01:23.901000+00:00',end_timestamp='1976-04-20 15:30:33.891000+00:00',start_station='S12',end_station='S17',start_orig_no=62241,end_orig_no=4121,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.49.csv.gz',start_timestamp ='1976-04-20 15:01:23.901000+00:00',end_timestamp='1976-04-20 15:30:33.891000+00:00',start_station='S12',end_station='S17',start_orig_no=65280,end_orig_no=4121,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.50.csv.gz',start_timestamp ='1976-04-21 12:12:10.384000+00:00',end_timestamp='1976-04-21 14:41:44.583000+00:00',start_station='S12',end_station='S17',start_orig_no=543,end_orig_no=1533,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.51.csv.gz',start_timestamp ='1976-04-21 15:03:29.292000+00:00',end_timestamp='1976-04-21 16:31:55.779000+00:00',start_station='S12',end_station='S17',start_orig_no=1678,end_orig_no=2263,orig_ground_station1=8)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.52.csv.gz',start_timestamp ='1976-04-22 11:08:45.378000+00:00',end_timestamp='1976-04-22 12:27:06.436000+00:00',start_station='S12',end_station='S15',start_orig_no=143,end_orig_no=1042,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.54.csv.gz',start_timestamp ='1976-04-24 13:40:38.976000+00:00',end_timestamp='1976-04-24 17:40:32.037000+00:00',start_station='S12',end_station='S17',start_orig_no=2856,end_orig_no=1590,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.2.9.csv.gz',start_timestamp ='1976-03-31 15:20:59.188000+00:00',end_timestamp='1976-03-31 17:54:51.205000+00:00',start_station='S12',end_station='S17',start_orig_no=3112,end_orig_no=4130,orig_ground_station1=1)
    
    # YYYY
    # df = reset_all_ground_stations_idx(df, gzip_filename, problem_gzip_filename='wtn.20.10.csv.gz',orig_idx_start=0,orig_idx_end=10)
    df = reset_all_ground_stations_idx(df, gzip_filename, problem_gzip_filename='wtn.20.10.csv.gz',orig_idx_start =109080,orig_idx_end=114160,single_station=None)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.10.csv.gz',start_timestamp ='1977-08-31 14:14:40.653000+00:00',end_timestamp='1977-08-31 14:29:34.955000+00:00',start_station='S12',end_station='S14',start_orig_no=4019,end_orig_no=4103,orig_ground_station1=4)
    
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.14.csv.gz',start_timestamp ='1977-09-02 18:55:55.368000+00:00',end_timestamp='1977-09-02 20:55:58.245000+00:00',start_station='S12',end_station='S14',start_orig_no=10,end_orig_no=1002,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.17.csv.gz',start_timestamp ='1977-09-04 16:42:44.599000+00:00',end_timestamp='1977-09-04 17:14:58.745000+00:00',start_station='S12',end_station='S17',start_orig_no=2058,end_orig_no=2324,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.29.csv.gz',start_timestamp ='1977-09-14 13:56:27.413000+00:00',end_timestamp='1977-09-14 17:17:10.092000+00:00',start_station='S12',end_station='S17',start_orig_no=24576,end_orig_no=2200,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.37.csv.gz',start_timestamp ='1977-09-18 17:44:47.329000+00:00',end_timestamp='1977-09-18 17:44:56.878000+00:00',start_station='S12',end_station='S17',start_orig_no=2484,end_orig_no=2865,orig_ground_station1=2)
    # XFAIL
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.39.csv.gz',start_timestamp ='1977-09-19 10:55:37.702000+00:00',end_timestamp='1977-09-19 11:51:53.679000+00:00',start_station='S12',end_station='S17',start_orig_no=1791,end_orig_no=2256,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.40.csv.gz',start_timestamp ='1977-09-19 21:22:18.213000+00:00',end_timestamp='1977-09-19 22:24:52.673000+00:00',start_station='S12',end_station='S17',start_orig_no=2861,end_orig_no=996,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.46.csv.gz',start_timestamp ='1977-09-23 21:18:13.988000+00:00',end_timestamp='1977-09-23 22:18:57.305000+00:00',start_station='S12',end_station='S17',start_orig_no=2471,end_orig_no=2973,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.50.csv.gz',start_timestamp ='1977-09-26 09:25:49.431000+00:00',end_timestamp='1977-09-26 09:39:58.575000+00:00',start_station='S12',end_station='S17',start_orig_no=2922,end_orig_no=3037,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.51.csv.gz',start_timestamp ='1977-09-27 00:33:20.363000+00:00',end_timestamp='1977-09-27 02:59:54.843000+00:00',start_station='S12',end_station='S17',start_orig_no=2653,end_orig_no=3866,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.52.csv.gz',start_timestamp ='1977-09-27 16:45:58.829000+00:00',end_timestamp='1977-09-27 16:54:54.645000+00:00',start_station='S12',end_station='S17',start_orig_no=2399,end_orig_no=2472,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.52.csv.gz',start_timestamp ='1977-09-27 19:03:41.847000+00:00',end_timestamp='1977-09-27 19:04:56.598000+00:00',start_station='S12',end_station='S17',start_orig_no=34578,end_orig_no=3584,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.54.csv.gz',start_timestamp ='1977-09-28 09:59:47.523000+00:00',end_timestamp='1977-09-28 09:59:56.648000+00:00',start_station='S12',end_station='S17',start_orig_no=2853,end_orig_no=2859,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.8.csv.gz',start_timestamp ='1977-08-30 00:19:52.532000+00:00',end_timestamp='1977-08-30 02:14:56.670000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=3153,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.20.8.csv.gz',start_timestamp ='1977-08-30 04:28:41.611000+00:00',end_timestamp='1977-08-30 06:59:54.258000+00:00',start_station='S12',end_station='S17',start_orig_no=33,end_orig_no=1283,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.10.csv.gz',start_timestamp ='1976-05-27 15:18:13.128000+00:00',end_timestamp='1976-05-27 16:54:56.193000+00:00',start_station='S12',end_station='S17',start_orig_no=2849,end_orig_no=3649,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.10.csv.gz',start_timestamp ='1976-05-27 19:44:08.344000+00:00',end_timestamp='1976-05-27 23:23:55.784000+00:00',start_station='S12',end_station='S17',start_orig_no=1437,end_orig_no=3268,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.10.csv.gz',start_timestamp ='1976-05-27 20:41:33.324000+00:00',end_timestamp='1976-05-27 23:23:55.784000+00:00',start_station='S12',end_station='S17',start_orig_no=1924,end_orig_no=3268,orig_ground_station1=4)



    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.12.csv.gz',start_timestamp ='1976-05-29 18:34:04.729000+00:00',end_timestamp='1976-05-29 18:52:54.368000+00:00',start_station='S12',end_station='S17',start_orig_no=598,end_orig_no=753,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.13.csv.gz',start_timestamp ='1976-06-11 17:55:57.380000+00:00',end_timestamp='1976-06-11 18:39:53.605000+00:00',start_station='S12',end_station='S17',start_orig_no=3351,end_orig_no=3714,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.13.csv.gz',start_timestamp ='1976-06-11 18:39:47.362000+00:00',end_timestamp='1976-06-11 18:39:53.605000+00:00',start_station='S12',end_station='S17',start_orig_no=3714,end_orig_no=3714,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.13.csv.gz',start_timestamp ='1976-06-11 23:55:45.423000+00:00',end_timestamp='1976-06-12 03:40:58.453000+00:00',start_station='S12',end_station='S17',start_orig_no=2187,end_orig_no=4057,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.13.csv.gz',start_timestamp ='1976-06-12 00:05:40.731000+00:00',end_timestamp='1976-06-12 03:40:58.453000+00:00',start_station='S12',end_station='S17',start_orig_no=2275,end_orig_no=4057,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.15.csv.gz',start_timestamp ='1976-06-15 10:51:02.244000+00:00',end_timestamp='1976-06-15 13:41:53.453000+00:00',start_station='S12',end_station='S17',start_orig_no=1756,end_orig_no=3170,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.15.csv.gz',start_timestamp ='1976-06-15 11:44:39.086000+00:00',end_timestamp='1976-06-15 13:41:53.453000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3170,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.15.csv.gz',start_timestamp ='1976-06-15 11:44:46.332000+00:00',end_timestamp='1976-06-15 13:41:53.453000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=3170,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.17.csv.gz',start_timestamp ='1976-06-16 09:23:51.994000+00:00',end_timestamp='1976-06-16 09:23:58.412000+00:00',start_station='S12',end_station='S17',start_orig_no=3492,end_orig_no=3492,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.17.csv.gz',start_timestamp ='1976-06-16 09:23:51.994000+00:00',end_timestamp='1976-06-16 09:23:58.412000+00:00',start_station='S12',end_station='S17',start_orig_no=3492,end_orig_no=3492,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.17.csv.gz',start_timestamp ='1976-06-16 09:26:26.547000+00:00',end_timestamp='1976-06-16 13:40:58.734000+00:00',start_station='S12',end_station='S17',start_orig_no=70,end_orig_no=2177,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.17.csv.gz',start_timestamp ='1976-06-16 13:32:56.915000+00:00',end_timestamp='1976-06-16 13:33:03.564000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2200,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.17.csv.gz',start_timestamp ='1976-06-16 13:32:56.915000+00:00',end_timestamp='1976-06-16 13:33:03.564000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2200,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.17.csv.gz',start_timestamp ='1976-06-16 13:38:48.304000+00:00',end_timestamp='1976-06-16 17:08:53.850000+00:00',start_station='S12',end_station='S17',start_orig_no=74,end_orig_no=1813,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 04:13:12.588000+00:00',end_timestamp='1976-06-18 06:55:59.040000+00:00',start_station='S12',end_station='S17',start_orig_no=6148,end_orig_no=3005,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 07:57:46.659000+00:00',end_timestamp='1976-06-18 07:57:52.970000+00:00',start_station='S12',end_station='S14',start_orig_no=3552,end_orig_no=3552,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 08:03:07.868000+00:00',end_timestamp='1976-06-18 14:04:53.203000+00:00',start_station='S12',end_station='S17',start_orig_no=66,end_orig_no=3061,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 12:20:48.281000+00:00',end_timestamp='1976-06-18 14:04:53.203000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3061,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 12:20:48.281000+00:00',end_timestamp='1976-06-18 14:04:53.203000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3061,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 12:20:48.281000+00:00',end_timestamp='1976-06-18 14:04:53.203000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3061,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.18.csv.gz',start_timestamp ='1976-06-18 12:20:48.281000+00:00',end_timestamp='1976-06-18 14:04:53.203000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3061,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.19.csv.gz',start_timestamp ='1976-06-21 11:16:18.845000+00:00',end_timestamp='1976-06-21 18:14:52.956000+00:00',start_station='S12',end_station='S17',start_orig_no=635,end_orig_no=4110,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.19.csv.gz',start_timestamp ='1976-06-21 11:52:15.431000+00:00',end_timestamp='1976-06-21 18:14:52.956000+00:00',start_station='S12',end_station='S17',start_orig_no=59804,end_orig_no=4110,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.19.csv.gz',start_timestamp ='1976-06-21 13:30:46.120000+00:00',end_timestamp='1976-06-21 18:14:52.956000+00:00',start_station='S12',end_station='S17',start_orig_no=1758,end_orig_no=4110,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.19.csv.gz',start_timestamp ='1976-06-21 18:14:46.290000+00:00',end_timestamp='1976-06-21 18:14:52.956000+00:00',start_station='S12',end_station='S17',start_orig_no=4110,end_orig_no=4110,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.19.csv.gz',start_timestamp ='1976-06-21 18:14:46.290000+00:00',end_timestamp='1976-06-21 18:14:52.956000+00:00',start_station='S12',end_station='S17',start_orig_no=4110,end_orig_no=4110,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.19.csv.gz',start_timestamp ='1976-06-21 18:10:52.639000+00:00',end_timestamp='1976-06-21 20:21:52.410000+00:00',start_station='S12',end_station='S17',start_orig_no=1,end_orig_no=1085,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.19.csv.gz',start_timestamp ='1976-06-21 18:19:27.032000+00:00',end_timestamp='1976-06-21 20:21:52.410000+00:00',start_station='S12',end_station='S17',start_orig_no=72,end_orig_no=1085,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.2.csv.gz',start_timestamp ='1977-09-29 09:58:27.716000+00:00',end_timestamp='1977-09-29 11:24:54.509000+00:00',start_station='S12',end_station='S17',start_orig_no=2982,end_orig_no=3697,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.2.csv.gz',start_timestamp ='1977-09-29 17:55:48.105000+00:00',end_timestamp='1977-09-29 17:55:54.622000+00:00',start_station='S12',end_station='S17',start_orig_no=2076,end_orig_no=2076,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.2.csv.gz',start_timestamp ='1977-09-29 17:55:48.105000+00:00',end_timestamp='1977-09-29 17:55:54.622000+00:00',start_station='S12',end_station='S17',start_orig_no=2076,end_orig_no=2076,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.2.csv.gz',start_timestamp ='1977-09-29 17:55:48.105000+00:00',end_timestamp='1977-09-29 17:55:54.622000+00:00',start_station='S12',end_station='S17',start_orig_no=2076,end_orig_no=2076,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.2.csv.gz',start_timestamp ='1977-09-29 17:55:48.105000+00:00',end_timestamp='1977-09-29 17:55:54.622000+00:00',start_station='S12',end_station='S17',start_orig_no=2076,end_orig_no=2076,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 05:39:13.758000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=1507,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 05:48:46.110000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=1586,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 05:59:23.666000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=1674,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 05:59:52.646000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=1678,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 06:00:07.136000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=1680,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 06:04:35.199000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=1717,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 06:37:18.583000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=1988,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:02:52.100000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2852,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 06:58:39.129000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 06:58:39.129000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2201,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.20.csv.gz',start_timestamp ='1976-06-22 07:04:12.397000+00:00',end_timestamp='1976-06-22 09:00:14.254000+00:00',start_station='S12',end_station='S17',start_orig_no=2247,end_orig_no=3207,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.21.csv.gz',start_timestamp ='1976-06-22 17:37:45.130000+00:00',end_timestamp='1976-06-22 18:39:55.821000+00:00',start_station='S12',end_station='S17',start_orig_no=3373,end_orig_no=3887,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.21.csv.gz',start_timestamp ='1976-06-22 17:42:05.949000+00:00',end_timestamp='1976-06-22 18:39:55.821000+00:00',start_station='S12',end_station='S17',start_orig_no=3409,end_orig_no=3887,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.22.csv.gz',start_timestamp ='1976-09-02 12:19:34.072000+00:00',end_timestamp='1976-09-02 14:43:51.103000+00:00',start_station='S12',end_station='S14',start_orig_no=3318,end_orig_no=4273,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 14:24:30.203000+00:00',end_timestamp='1977-03-23 16:48:56.050000+00:00',start_station='S12',end_station='S14',start_orig_no=1702,end_orig_no=2658,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 14:34:27.948000+00:00',end_timestamp='1977-03-23 16:48:56.050000+00:00',start_station='S12',end_station='S14',start_orig_no=1768,end_orig_no=2658,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 14:47:26.828000+00:00',end_timestamp='1977-03-23 16:48:56.050000+00:00',start_station='S12',end_station='S14',start_orig_no=1854,end_orig_no=2658,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 14:56:12.119000+00:00',end_timestamp='1977-03-23 16:48:56.050000+00:00',start_station='S12',end_station='S14',start_orig_no=1912,end_orig_no=2658,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 15:39:39.857000+00:00',end_timestamp='1977-03-23 16:48:56.050000+00:00',start_station='S12',end_station='S14',start_orig_no=2200,end_orig_no=2658,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 15:39:39.857000+00:00',end_timestamp='1977-03-23 16:48:56.050000+00:00',start_station='S12',end_station='S14',start_orig_no=2200,end_orig_no=2658,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 15:39:48.914000+00:00',end_timestamp='1977-03-23 16:48:56.050000+00:00',start_station='S12',end_station='S14',start_orig_no=2201,end_orig_no=2658,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 15:50:41+00:00',end_timestamp='1977-03-23 16:48:56.050000+00:00',start_station='S12',end_station='S14',start_orig_no=2273,end_orig_no=2658,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 22:05:37.955000+00:00',end_timestamp='1977-03-24 00:57:50.281000+00:00',start_station='S12',end_station='S14',start_orig_no=2200,end_orig_no=3340,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 22:05:37.955000+00:00',end_timestamp='1977-03-24 00:57:50.281000+00:00',start_station='S12',end_station='S14',start_orig_no=2200,end_orig_no=3340,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 22:05:47.012000+00:00',end_timestamp='1977-03-24 00:57:50.281000+00:00',start_station='S12',end_station='S14',start_orig_no=2201,end_orig_no=3340,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.25.csv.gz',start_timestamp ='1977-03-23 22:16:30.041000+00:00',end_timestamp='1977-03-24 00:57:50.281000+00:00',start_station='S12',end_station='S14',start_orig_no=2272,end_orig_no=3340,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.27.csv.gz',start_timestamp ='1977-03-25 17:26:47.434000+00:00',end_timestamp='1977-03-25 20:44:50.996000+00:00',start_station='S12',end_station='S14',start_orig_no=64799,end_orig_no=3383,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.28.csv.gz',start_timestamp ='1977-04-25 19:34:39.746000+00:00',end_timestamp='1977-04-25 22:59:55.589000+00:00',start_station='S12',end_station='S17',start_orig_no=47,end_orig_no=1746,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.28.csv.gz',start_timestamp ='1977-04-25 23:22:22.780000+00:00',end_timestamp='1977-04-26 02:39:54.597000+00:00',start_station='S12',end_station='S17',start_orig_no=2236,end_orig_no=3871,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.28.csv.gz',start_timestamp ='1977-04-26 02:40:48.804000+00:00',end_timestamp='1977-04-26 04:59:54.274000+00:00',start_station='S12',end_station='S17',start_orig_no=45,end_orig_no=1196,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.33.csv.gz',start_timestamp ='1977-09-07 05:53:46.222000+00:00',end_timestamp='1977-09-07 06:21:53.715000+00:00',start_station='S12',end_station='S17',start_orig_no=1787,end_orig_no=2019,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.33.csv.gz',start_timestamp ='1977-09-07 05:53:46.222000+00:00',end_timestamp='1977-09-07 06:21:53.715000+00:00',start_station='S12',end_station='S17',start_orig_no=1787,end_orig_no=2019,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.33.csv.gz',start_timestamp ='1977-09-07 05:54:00.712000+00:00',end_timestamp='1977-09-07 06:21:53.715000+00:00',start_station='S12',end_station='S17',start_orig_no=1789,end_orig_no=2019,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.33.csv.gz',start_timestamp ='1977-09-07 06:02:06.141000+00:00',end_timestamp='1977-09-07 06:21:53.715000+00:00',start_station='S12',end_station='S17',start_orig_no=1856,end_orig_no=2019,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.33.csv.gz',start_timestamp ='1977-09-07 06:40:15.154000+00:00',end_timestamp='1977-09-07 08:40:36.526000+00:00',start_station='S12',end_station='S17',start_orig_no=2047,end_orig_no=3043,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.33.csv.gz',start_timestamp ='1977-09-07 07:02:55.926000+00:00',end_timestamp='1977-09-07 08:40:36.526000+00:00',start_station='S12',end_station='S17',start_orig_no=2235,end_orig_no=3043,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.33.csv.gz',start_timestamp ='1977-09-07 07:09:34.412000+00:00',end_timestamp='1977-09-07 08:40:36.526000+00:00',start_station='S12',end_station='S17',start_orig_no=2290,end_orig_no=3043,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.33.csv.gz',start_timestamp ='1977-09-07 07:18:23.312000+00:00',end_timestamp='1977-09-07 08:40:36.526000+00:00',start_station='S12',end_station='S17',start_orig_no=2363,end_orig_no=3043,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 02:29:09.719000+00:00',end_timestamp='1977-09-12 02:55:52.220000+00:00',start_station='S12',end_station='S17',start_orig_no=2131,end_orig_no=2364,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 02:29:09.719000+00:00',end_timestamp='1977-09-12 02:55:52.220000+00:00',start_station='S12',end_station='S17',start_orig_no=2131,end_orig_no=2364,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 02:30:36.667000+00:00',end_timestamp='1977-09-12 02:55:52.220000+00:00',start_station='S12',end_station='S17',start_orig_no=2143,end_orig_no=2364,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 02:46:42.144000+00:00',end_timestamp='1977-09-12 02:55:52.220000+00:00',start_station='S12',end_station='S17',start_orig_no=2289,end_orig_no=2364,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 02:46:56.635000+00:00',end_timestamp='1977-09-12 02:55:52.220000+00:00',start_station='S12',end_station='S17',start_orig_no=2291,end_orig_no=2364,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 03:05:34.260000+00:00',end_timestamp='1977-09-12 04:36:36.441000+00:00',start_station='S12',end_station='S17',start_orig_no=2494,end_orig_no=3247,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 03:05:41.506000+00:00',end_timestamp='1977-09-12 04:36:36.441000+00:00',start_station='S12',end_station='S17',start_orig_no=2495,end_orig_no=3247,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 03:28:16.434000+00:00',end_timestamp='1977-09-12 04:36:36.441000+00:00',start_station='S12',end_station='S17',start_orig_no=2682,end_orig_no=3247,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 03:30:48.592000+00:00',end_timestamp='1977-09-12 04:36:36.441000+00:00',start_station='S12',end_station='S17',start_orig_no=2703,end_orig_no=3247,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.35.csv.gz',start_timestamp ='1977-09-12 03:39:37.521000+00:00',end_timestamp='1977-09-12 04:36:36.441000+00:00',start_station='S12',end_station='S17',start_orig_no=2776,end_orig_no=3247,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.36.csv.gz',start_timestamp ='1977-09-12 04:36:35.634000+00:00',end_timestamp='1977-09-12 06:04:52.614000+00:00',start_station='S12',end_station='S17',start_orig_no=1006,end_orig_no=3978,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.36.csv.gz',start_timestamp ='1977-09-12 11:01:50.294000+00:00',end_timestamp='1977-09-12 13:44:57.636000+00:00',start_station='S12',end_station='S17',start_orig_no=100,end_orig_no=1450,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.36.csv.gz',start_timestamp ='1977-09-12 11:03:02.750000+00:00',end_timestamp='1977-09-12 13:44:57.636000+00:00',start_station='S12',end_station='S17',start_orig_no=110,end_orig_no=1450,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.36.csv.gz',start_timestamp ='1977-09-12 11:03:02.750000+00:00',end_timestamp='1977-09-12 13:44:57.636000+00:00',start_station='S12',end_station='S17',start_orig_no=110,end_orig_no=1450,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.36.csv.gz',start_timestamp ='1977-09-12 11:03:02.750000+00:00',end_timestamp='1977-09-12 13:44:57.636000+00:00',start_station='S12',end_station='S17',start_orig_no=110,end_orig_no=1450,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.36.csv.gz',start_timestamp ='1977-09-12 11:03:02.750000+00:00',end_timestamp='1977-09-12 13:44:57.636000+00:00',start_station='S12',end_station='S17',start_orig_no=110,end_orig_no=1450,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.36.csv.gz',start_timestamp ='1977-09-12 11:07:30.837000+00:00',end_timestamp='1977-09-12 13:44:57.636000+00:00',start_station='S12',end_station='S17',start_orig_no=147,end_orig_no=1450,orig_ground_station1=1)
    
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.37.csv.gz',start_timestamp ='1977-09-13 09:20:01.224000+00:00',end_timestamp='1977-09-13 12:59:57.862000+00:00',start_station='S12',end_station='S17',start_orig_no=3105,end_orig_no=4322,orig_ground_station1=1)
    
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.37.csv.gz',start_timestamp ='1977-09-13 09:49:50.282000+00:00',end_timestamp='1977-09-13 12:59:57.862000+00:00',start_station='S12',end_station='S17',start_orig_no=3352,end_orig_no=4322,orig_ground_station1=1)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.37.csv.gz',start_timestamp ='1977-09-13 12:21:56.112000+00:00',end_timestamp='1977-09-13 12:59:57.862000+00:00',start_station='S12',end_station='S17',start_orig_no=4008,end_orig_no=4322,orig_ground_station1=1)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.37.csv.gz',start_timestamp ='1977-09-13 12:47:03.197000+00:00',end_timestamp='1977-09-13 12:59:57.862000+00:00',start_station='S12',end_station='S17',start_orig_no=4216,end_orig_no=4322,orig_ground_station1=1)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.37.csv.gz',start_timestamp ='1977-09-13 14:04:10.730000+00:00',end_timestamp='1977-09-13 14:30:58.259000+00:00',start_station='S12',end_station='S17',start_orig_no=609,end_orig_no=830,orig_ground_station1=2)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.37.csv.gz',start_timestamp ='1977-09-13 15:39:55.886000+00:00',end_timestamp='1977-09-13 19:11:18.968000+00:00',start_station='S12',end_station='S17',start_orig_no=1494,end_orig_no=3265,orig_ground_station1=2)
    # 
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.38.csv.gz',start_timestamp ='1977-09-14 06:42:41.609000+00:00',end_timestamp='1977-09-14 08:39:55.800000+00:00',start_station='S12',end_station='S17',start_orig_no=2089,end_orig_no=3059,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.38.csv.gz',start_timestamp ='1977-09-14 06:44:08.557000+00:00',end_timestamp='1977-09-14 08:39:55.800000+00:00',start_station='S12',end_station='S17',start_orig_no=2101,end_orig_no=3059,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.38.csv.gz',start_timestamp ='1977-09-14 06:56:05.872000+00:00',end_timestamp='1977-09-14 08:39:55.800000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3059,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.38.csv.gz',start_timestamp ='1977-09-14 06:56:05.872000+00:00',end_timestamp='1977-09-14 08:39:55.800000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3059,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.38.csv.gz',start_timestamp ='1977-09-14 06:56:05.872000+00:00',end_timestamp='1977-09-14 08:39:55.800000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3059,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.38.csv.gz',start_timestamp ='1977-09-14 07:00:55.697000+00:00',end_timestamp='1977-09-14 08:39:55.800000+00:00',start_station='S12',end_station='S17',start_orig_no=2240,end_orig_no=3059,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.38.csv.gz',start_timestamp ='1977-09-14 07:50:40.286000+00:00',end_timestamp='1977-09-14 08:39:55.800000+00:00',start_station='S12',end_station='S17',start_orig_no=2652,end_orig_no=3059,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.38.csv.gz',start_timestamp ='1977-09-14 11:32:51.583000+00:00',end_timestamp='1977-09-14 12:17:53.190000+00:00',start_station='S12',end_station='S17',start_orig_no=3092,end_orig_no=3464,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.39.csv.gz',start_timestamp ='1977-09-15 05:20:33.003000+00:00',end_timestamp='1977-09-15 08:24:55.095000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3726,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.39.csv.gz',start_timestamp ='1977-09-15 05:20:33.003000+00:00',end_timestamp='1977-09-15 08:24:55.095000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=3726,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.39.csv.gz',start_timestamp ='1977-09-15 05:29:14.082000+00:00',end_timestamp='1977-09-15 08:24:55.095000+00:00',start_station='S12',end_station='S17',start_orig_no=2272,end_orig_no=3726,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:39.051000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3843,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:24.560000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3841,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:10.069000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3839,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:02.823000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3838,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:02.823000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3838,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:02.823000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3838,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:10.069000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3839,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:10.069000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3839,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:10.069000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3839,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:10.069000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3839,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:10.069000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3839,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:17.314000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3840,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:17.314000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3840,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:17.314000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3840,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:24.560000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3841,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:24.560000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3841,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:24.560000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3841,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:24.560000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3841,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:24.560000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3841,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:24.560000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3841,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:24.560000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3841,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:31.805000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3842,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:31.805000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3842,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:39.051000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3843,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.40.csv.gz',start_timestamp ='1977-09-18 02:39:46.297000+00:00',end_timestamp='1977-09-18 02:39:52.782000+00:00',start_station='S12',end_station='S17',start_orig_no=3844,end_orig_no=3844,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.41.csv.gz',start_timestamp ='1977-09-21 02:59:50.684000+00:00',end_timestamp='1977-09-21 03:59:58.383000+00:00',start_station='S12',end_station='S17',start_orig_no=496,end_orig_no=993,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.41.csv.gz',start_timestamp ='1977-09-21 06:23:53.958000+00:00',end_timestamp='1977-09-21 08:59:56.759000+00:00',start_station='S12',end_station='S17',start_orig_no=2864,end_orig_no=3331,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.42.csv.gz',start_timestamp ='1977-09-22 16:03:54.918000+00:00',end_timestamp='1977-09-22 20:44:53.515000+00:00',start_station='S12',end_station='S17',start_orig_no=224,end_orig_no=2550,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.42.csv.gz',start_timestamp ='1977-09-22 21:13:41.395000+00:00',end_timestamp='1977-09-22 23:59:57.617000+00:00',start_station='S12',end_station='S17',start_orig_no=2835,end_orig_no=4211,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.44.csv.gz',start_timestamp ='1977-09-24 01:58:39.571000+00:00',end_timestamp='1977-09-24 02:34:52.451000+00:00',start_station='S12',end_station='S17',start_orig_no=1822,end_orig_no=2121,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.44.csv.gz',start_timestamp ='1977-09-24 02:09:46.122000+00:00',end_timestamp='1977-09-24 02:34:52.451000+00:00',start_station='S12',end_station='S17',start_orig_no=1914,end_orig_no=2121,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.44.csv.gz',start_timestamp ='1977-09-24 02:34:22.315000+00:00',end_timestamp='1977-09-24 02:38:48.525000+00:00',start_station='S12',end_station='S17',start_orig_no=2164,end_orig_no=2200,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.44.csv.gz',start_timestamp ='1977-09-24 03:11:55.541000+00:00',end_timestamp='1977-09-24 04:14:06.932000+00:00',start_station='S12',end_station='S14',start_orig_no=277,end_orig_no=792,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 11:20:44.910000+00:00',end_timestamp='1977-09-24 11:29:54.407000+00:00',start_station='S12',end_station='S17',start_orig_no=2645,end_orig_no=2720,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 13:12:31.453000+00:00',end_timestamp='1977-09-24 13:37:51.979000+00:00',start_station='S12',end_station='S17',start_orig_no=3321,end_orig_no=3530,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 13:34:08.327000+00:00',end_timestamp='1977-09-24 16:59:53.329000+00:00',start_station='S12',end_station='S17',start_orig_no=12,end_orig_no=1715,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 13:59:22.553000+00:00',end_timestamp='1977-09-24 16:59:53.329000+00:00',start_station='S12',end_station='S17',start_orig_no=221,end_orig_no=1715,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 14:29:19.337000+00:00',end_timestamp='1977-09-24 16:59:53.329000+00:00',start_station='S12',end_station='S17',start_orig_no=469,end_orig_no=1715,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 15:22:19.936000+00:00',end_timestamp='1977-09-24 16:59:53.329000+00:00',start_station='S12',end_station='S17',start_orig_no=908,end_orig_no=1715,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 16:09:53.901000+00:00',end_timestamp='1977-09-24 16:59:53.329000+00:00',start_station='S12',end_station='S17',start_orig_no=1302,end_orig_no=1715,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 16:38:52.725000+00:00',end_timestamp='1977-09-24 16:59:53.329000+00:00',start_station='S12',end_station='S17',start_orig_no=1542,end_orig_no=1715,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 16:57:04.923000+00:00',end_timestamp='1977-09-24 16:59:53.329000+00:00',start_station='S12',end_station='S17',start_orig_no=1694,end_orig_no=1715,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 17:03:01.746000+00:00',end_timestamp='1977-09-24 18:21:52.264000+00:00',start_station='S12',end_station='S17',start_orig_no=1781,end_orig_no=2433,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 17:32:15.057000+00:00',end_timestamp='1977-09-24 18:21:52.264000+00:00',start_station='S12',end_station='S17',start_orig_no=2023,end_orig_no=2433,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.45.csv.gz',start_timestamp ='1977-09-24 18:20:04.109000+00:00',end_timestamp='1977-09-24 18:21:52.264000+00:00',start_station='S12',end_station='S17',start_orig_no=2419,end_orig_no=2433,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 07:54:56.943000+00:00',end_timestamp='1977-09-25 08:35:10.038000+00:00',start_station='S12',end_station='S14',start_orig_no=3359,end_orig_no=3692,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 08:22:43.311000+00:00',end_timestamp='1977-09-25 08:35:10.038000+00:00',start_station='S12',end_station='S14',start_orig_no=3589,end_orig_no=3692,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 08:37:53.180000+00:00',end_timestamp='1977-09-25 11:33:55.792000+00:00',start_station='S12',end_station='S17',start_orig_no=59,end_orig_no=1516,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 08:46:20.335000+00:00',end_timestamp='1977-09-25 11:33:55.792000+00:00',start_station='S12',end_station='S17',start_orig_no=129,end_orig_no=1516,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 09:14:06.702000+00:00',end_timestamp='1977-09-25 11:33:55.792000+00:00',start_station='S12',end_station='S17',start_orig_no=359,end_orig_no=1516,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 09:25:49.474000+00:00',end_timestamp='1977-09-25 11:33:55.792000+00:00',start_station='S12',end_station='S17',start_orig_no=456,end_orig_no=1516,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 09:38:44.697000+00:00',end_timestamp='1977-09-25 11:33:55.792000+00:00',start_station='S12',end_station='S17',start_orig_no=563,end_orig_no=1516,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 09:55:10.027000+00:00',end_timestamp='1977-09-25 11:33:55.792000+00:00',start_station='S12',end_station='S17',start_orig_no=699,end_orig_no=1516,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 10:18:13.836000+00:00',end_timestamp='1977-09-25 11:33:55.792000+00:00',start_station='S12',end_station='S17',start_orig_no=890,end_orig_no=1516,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 10:51:11.741000+00:00',end_timestamp='1977-09-25 11:33:55.792000+00:00',start_station='S12',end_station='S17',start_orig_no=1163,end_orig_no=1516,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 11:26:34.549000+00:00',end_timestamp='1977-09-25 11:33:55.792000+00:00',start_station='S12',end_station='S17',start_orig_no=1456,end_orig_no=1516,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 11:44:34.064000+00:00',end_timestamp='1977-09-25 13:32:36.254000+00:00',start_station='S12',end_station='S14',start_orig_no=1646,end_orig_no=2540,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 11:44:41.309000+00:00',end_timestamp='1977-09-25 13:32:36.254000+00:00',start_station='S12',end_station='S14',start_orig_no=1647,end_orig_no=2540,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 12:01:13.884000+00:00',end_timestamp='1977-09-25 13:32:36.254000+00:00',start_station='S12',end_station='S14',start_orig_no=1784,end_orig_no=2540,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 12:17:02.989000+00:00',end_timestamp='1977-09-25 13:32:36.254000+00:00',start_station='S12',end_station='S14',start_orig_no=1915,end_orig_no=2540,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 12:41:33.740000+00:00',end_timestamp='1977-09-25 13:32:36.254000+00:00',start_station='S12',end_station='S14',start_orig_no=2118,end_orig_no=2540,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 13:14:17.156000+00:00',end_timestamp='1977-09-25 13:32:36.254000+00:00',start_station='S12',end_station='S14',start_orig_no=2389,end_orig_no=2540,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 13:30:35.242000+00:00',end_timestamp='1977-09-25 13:32:36.254000+00:00',start_station='S12',end_station='S14',start_orig_no=2524,end_orig_no=2540,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 13:48:42.610000+00:00',end_timestamp='1977-09-25 15:59:56.534000+00:00',start_station='S12',end_station='S17',start_orig_no=2715,end_orig_no=3801,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 14:09:07.029000+00:00',end_timestamp='1977-09-25 15:59:56.534000+00:00',start_station='S12',end_station='S17',start_orig_no=2884,end_orig_no=3801,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 14:15:09.283000+00:00',end_timestamp='1977-09-25 15:59:56.534000+00:00',start_station='S12',end_station='S17',start_orig_no=2934,end_orig_no=3801,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 14:38:34.830000+00:00',end_timestamp='1977-09-25 15:59:56.534000+00:00',start_station='S12',end_station='S17',start_orig_no=3128,end_orig_no=3801,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 14:39:54.525000+00:00',end_timestamp='1977-09-25 15:59:56.534000+00:00',start_station='S12',end_station='S17',start_orig_no=3139,end_orig_no=3801,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 15:03:56.296000+00:00',end_timestamp='1977-09-25 15:59:56.534000+00:00',start_station='S12',end_station='S17',start_orig_no=3338,end_orig_no=3801,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 15:27:06.748000+00:00',end_timestamp='1977-09-25 15:59:56.534000+00:00',start_station='S12',end_station='S17',start_orig_no=3530,end_orig_no=3801,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 15:57:32.509000+00:00',end_timestamp='1977-09-25 15:59:56.534000+00:00',start_station='S12',end_station='S17',start_orig_no=3782,end_orig_no=3801,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 16:17:55.720000+00:00',end_timestamp='1977-09-25 17:19:29.966000+00:00',start_station='S12',end_station='S17',start_orig_no=191,end_orig_no=700,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 16:44:07.903000+00:00',end_timestamp='1977-09-25 17:19:29.966000+00:00',start_station='S12',end_station='S17',start_orig_no=408,end_orig_no=700,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 16:58:51.804000+00:00',end_timestamp='1977-09-25 17:19:29.966000+00:00',start_station='S12',end_station='S17',start_orig_no=530,end_orig_no=700,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.46.csv.gz',start_timestamp ='1977-09-25 17:12:37.744000+00:00',end_timestamp='1977-09-25 17:19:29.966000+00:00',start_station='S12',end_station='S17',start_orig_no=644,end_orig_no=700,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.48.csv.gz',start_timestamp ='1977-09-29 19:48:00.575000+00:00',end_timestamp='1977-09-29 19:48:02.426000+00:00',start_station='S12',end_station='S17',start_orig_no=1772,end_orig_no=1772,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.48.csv.gz',start_timestamp ='1977-09-29 23:24:44.873000+00:00',end_timestamp='1977-09-30 02:44:56.016000+00:00',start_station='S12',end_station='S17',start_orig_no=2671,end_orig_no=4328,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.7.csv.gz',start_timestamp ='1976-04-23 08:29:59.140000+00:00',end_timestamp='1976-04-23 09:29:54.069000+00:00',start_station='S12',end_station='S17',start_orig_no=1234,end_orig_no=1630,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.7.csv.gz',start_timestamp ='1976-04-23 08:38:35.668000+00:00',end_timestamp='1976-04-23 09:29:54.069000+00:00',start_station='S12',end_station='S17',start_orig_no=1291,end_orig_no=1630,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 04:45:59.067000+00:00',end_timestamp='1976-04-30 06:25:58.070000+00:00',start_station='S12',end_station='S17',start_orig_no=3167,end_orig_no=3609,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:11:20.249000+00:00',end_timestamp='1976-04-30 06:25:58.070000+00:00',start_station='S12',end_station='S17',start_orig_no=3513,end_orig_no=3609,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:20:50.791000+00:00',end_timestamp='1976-04-30 06:25:58.070000+00:00',start_station='S12',end_station='S17',start_orig_no=3576,end_orig_no=3609,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:23:24.747000+00:00',end_timestamp='1976-04-30 06:25:58.070000+00:00',start_station='S12',end_station='S17',start_orig_no=3593,end_orig_no=3609,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:24:19.085000+00:00',end_timestamp='1976-04-30 06:25:58.070000+00:00',start_station='S12',end_station='S17',start_orig_no=3599,end_orig_no=3609,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:24:45.065000+00:00',end_timestamp='1976-04-30 08:17:02.203000+00:00',start_station='S12',end_station='S17',start_orig_no=6,end_orig_no=749,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:27:00.908000+00:00',end_timestamp='1976-04-30 08:17:02.203000+00:00',start_station='S12',end_station='S17',start_orig_no=21,end_orig_no=749,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:29:43.920000+00:00',end_timestamp='1976-04-30 08:17:02.203000+00:00',start_station='S12',end_station='S17',start_orig_no=39,end_orig_no=749,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:30:11.089000+00:00',end_timestamp='1976-04-30 08:17:02.203000+00:00',start_station='S12',end_station='S17',start_orig_no=42,end_orig_no=749,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:32:54.101000+00:00',end_timestamp='1976-04-30 08:17:02.203000+00:00',start_station='S12',end_station='S17',start_orig_no=60,end_orig_no=749,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:33:48.438000+00:00',end_timestamp='1976-04-30 08:17:02.203000+00:00',start_station='S12',end_station='S17',start_orig_no=66,end_orig_no=749,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 06:33:57.494000+00:00',end_timestamp='1976-04-30 08:17:02.203000+00:00',start_station='S12',end_station='S17',start_orig_no=67,end_orig_no=749,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 07:09:07.590000+00:00',end_timestamp='1976-04-30 08:17:02.203000+00:00',start_station='S12',end_station='S17',start_orig_no=300,end_orig_no=749,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.8.csv.gz',start_timestamp ='1976-04-30 07:55:00.678000+00:00',end_timestamp='1976-04-30 08:17:02.203000+00:00',start_station='S12',end_station='S17',start_orig_no=604,end_orig_no=749,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:44.420000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2199,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:03:53.477000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2200,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:13:14.961000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2262,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:58:13.715000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2560,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:59:44.277000+00:00',end_timestamp='1976-05-02 04:59:53.268000+00:00',start_station='S12',end_station='S17',start_orig_no=2570,end_orig_no=2570,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 04:57:55.006000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2576,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 05:00:28.963000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2593,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 05:01:41.412000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2601,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 05:02:26.693000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2606,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 05:06:04.043000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2630,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 05:06:49.324000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2635,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 05:07:43.661000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2641,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 05:08:12.037000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2644,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 05:09:32.336000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2653,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 05:13:27.798000+00:00',end_timestamp='1976-05-02 06:13:58.765000+00:00',start_station='S12',end_station='S17',start_orig_no=2679,end_orig_no=3079,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.21.9.csv.gz',start_timestamp ='1976-05-02 09:36:47.890000+00:00',end_timestamp='1976-05-02 09:37:41.896000+00:00',start_station='S12',end_station='S17',start_orig_no=1390,end_orig_no=1395,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.12.csv.gz',start_timestamp ='1976-05-03 12:36:12.653000+00:00',end_timestamp='1976-05-03 12:54:54.790000+00:00',start_station='S12',end_station='S17',start_orig_no=1559,end_orig_no=1682,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.13.csv.gz',start_timestamp ='1976-05-03 21:18:06.024000+00:00',end_timestamp='1976-05-03 21:32:56.368000+00:00',start_station='S12',end_station='S17',start_orig_no=2858,end_orig_no=2297,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.13.csv.gz',start_timestamp ='1976-05-03 21:18:06.024000+00:00',end_timestamp='1976-05-03 21:32:56.368000+00:00',start_station='S12',end_station='S17',start_orig_no=213,end_orig_no=2297,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.14.csv.gz',start_timestamp ='1976-05-04 11:09:58.979000+00:00',end_timestamp='1976-05-04 11:46:57.685000+00:00',start_station='S12',end_station='S17',start_orig_no=2761,end_orig_no=3005,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.15.csv.gz',start_timestamp ='1976-05-05 11:01:02.021000+00:00',end_timestamp='1976-05-05 11:35:55.817000+00:00',start_station='S12',end_station='S12',start_orig_no=3261,end_orig_no=3468,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.15.csv.gz',start_timestamp ='1976-05-05 11:05:34.915000+00:00',end_timestamp='1976-05-05 11:35:55.817000+00:00',start_station='S12',end_station='S12',start_orig_no=3268,end_orig_no=3468,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.16.csv.gz',start_timestamp ='1976-05-06 01:58:46.472000+00:00',end_timestamp='1976-05-06 02:52:22.493000+00:00',start_station='S12',end_station='S17',start_orig_no=64783,end_orig_no=1614,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.19.csv.gz',start_timestamp ='1976-05-07 18:52:56.670000+00:00',end_timestamp='1976-05-07 23:20:25.554000+00:00',start_station='S12',end_station='S17',start_orig_no=1883,end_orig_no=3652,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.2.csv.gz',start_timestamp ='1976-04-26 13:10:42.315000+00:00',end_timestamp='1976-04-26 13:43:17.901000+00:00',start_station='S12',end_station='S17',start_orig_no=3252,end_orig_no=3467,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.20.csv.gz',start_timestamp ='1976-05-08 02:07:39.140000+00:00',end_timestamp='1976-05-08 02:32:53.876000+00:00',start_station='S12',end_station='S17',start_orig_no=11663,end_orig_no=756,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.20.csv.gz',start_timestamp ='1976-05-08 12:14:26.877000+00:00',end_timestamp='1976-05-08 14:26:56.200000+00:00',start_station='S12',end_station='S15',start_orig_no=1805,end_orig_no=2689,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.20.csv.gz',start_timestamp ='1976-05-08 12:43:51.630000+00:00',end_timestamp='1976-05-08 14:26:56.200000+00:00',start_station='S12',end_station='S15',start_orig_no=2007,end_orig_no=2689,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.21.csv.gz',start_timestamp ='1976-05-08 22:18:06.769000+00:00',end_timestamp='1976-05-09 00:53:52.427000+00:00',start_station='S12',end_station='S17',start_orig_no=3136,end_orig_no=4167,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.23.csv.gz',start_timestamp ='1976-05-09 21:22:53.268000+00:00',end_timestamp='1976-05-10 00:38:56.831000+00:00',start_station='S12',end_station='S17',start_orig_no=2420,end_orig_no=3718,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.24.csv.gz',start_timestamp ='1976-05-10 22:15:10.427000+00:00',end_timestamp='1976-05-11 00:19:50.176000+00:00',start_station='S12',end_station='S17',start_orig_no=3042,end_orig_no=3867,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.25.csv.gz',start_timestamp ='1976-05-11 09:26:39.742000+00:00',end_timestamp='1976-05-11 11:00:52.350000+00:00',start_station='S12',end_station='S17',start_orig_no=1528,end_orig_no=2150,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.25.csv.gz',start_timestamp ='1976-05-11 09:26:39.742000+00:00',end_timestamp='1976-05-11 11:00:52.350000+00:00',start_station='S12',end_station='S17',start_orig_no=1062,end_orig_no=2150,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.27.csv.gz',start_timestamp ='1976-05-12 07:45:21.391000+00:00',end_timestamp='1976-05-12 11:19:58.917000+00:00',start_station='S12',end_station='S17',start_orig_no=336,end_orig_no=1756,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.27.csv.gz',start_timestamp ='1976-05-12 12:21:50.897000+00:00',end_timestamp='1976-05-12 12:21:53.321000+00:00',start_station='S12',end_station='S17',start_orig_no=2871,end_orig_no=2871,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.28.csv.gz',start_timestamp ='1976-05-12 12:22:14.444000+00:00',end_timestamp='1976-05-12 12:22:16.868000+00:00',start_station='S12',end_station='S17',start_orig_no=238,end_orig_no=238,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.29.csv.gz',start_timestamp ='1976-05-13 08:16:13.636000+00:00',end_timestamp='1976-05-13 09:08:53.520000+00:00',start_station='S12',end_station='S17',start_orig_no=3123,end_orig_no=3471,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.29.csv.gz',start_timestamp ='1976-05-13 10:57:42.792000+00:00',end_timestamp='1976-05-13 12:54:51.232000+00:00',start_station='S12',end_station='S17',start_orig_no=16,end_orig_no=1523,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.3.csv.gz',start_timestamp ='1976-04-26 16:43:37.611000+00:00',end_timestamp='1976-04-26 19:24:52.564000+00:00',start_station='S12',end_station='S17',start_orig_no=887,end_orig_no=1953,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.3.csv.gz',start_timestamp ='1976-04-26 16:43:37.611000+00:00',end_timestamp='1976-04-26 19:24:52.564000+00:00',start_station='S12',end_station='S17',start_orig_no=62991,end_orig_no=1953,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.31.csv.gz',start_timestamp ='1976-05-14 15:58:51.995000+00:00',end_timestamp='1976-05-14 16:25:52.162000+00:00',start_station='S12',end_station='S17',start_orig_no=2408,end_orig_no=2586,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.32.csv.gz',start_timestamp ='1976-05-15 01:35:30.103000+00:00',end_timestamp='1976-05-15 05:40:57.363000+00:00',start_station='S12',end_station='S17',start_orig_no=2875,end_orig_no=3825,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.32.csv.gz',start_timestamp ='1976-05-15 13:13:00.797000+00:00',end_timestamp='1976-05-15 15:30:51.025000+00:00',start_station='S12',end_station='S17',start_orig_no=956,end_orig_no=1867,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.33.csv.gz',start_timestamp ='1976-05-16 05:23:20.472000+00:00',end_timestamp='1976-05-16 06:40:54.827000+00:00',start_station='S12',end_station='S17',start_orig_no=3314,end_orig_no=3827,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.35.csv.gz',start_timestamp ='1976-05-17 06:42:26.065000+00:00',end_timestamp='1976-05-17 07:54:52.722000+00:00',start_station='S12',end_station='S17',start_orig_no=3167,end_orig_no=3646,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.37.csv.gz',start_timestamp ='1976-05-18 07:36:58.554000+00:00',end_timestamp='1976-05-18 07:54:55.379000+00:00',start_station='S12',end_station='S17',start_orig_no=3267,end_orig_no=3385,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.38.csv.gz',start_timestamp ='1976-05-18 19:55:28.999000+00:00',end_timestamp='1976-05-18 20:22:56.700000+00:00',start_station='S12',end_station='S17',start_orig_no=2440,end_orig_no=2621,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.4.csv.gz',start_timestamp ='1976-04-27 16:15:44.354000+00:00',end_timestamp='1976-04-27 20:40:55.688000+00:00',start_station='S12',end_station='S17',start_orig_no=2531,end_orig_no=4287,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.40.csv.gz',start_timestamp ='1976-05-19 19:57:18.105000+00:00',end_timestamp='1976-05-19 22:56:54.343000+00:00',start_station='S12',end_station='S17',start_orig_no=2432,end_orig_no=3621,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.42.csv.gz',start_timestamp ='1976-05-20 23:07:01.293000+00:00',end_timestamp='1976-05-21 01:23:56.681000+00:00',start_station='S12',end_station='S17',start_orig_no=1066,end_orig_no=2199,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.43.csv.gz',start_timestamp ='1976-05-21 18:10:28.698000+00:00',end_timestamp='1976-05-21 19:10:53.490000+00:00',start_station='S12',end_station='S17',start_orig_no=2867,end_orig_no=2699,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.44.csv.gz',start_timestamp ='1976-05-22 06:40:49.316000+00:00',end_timestamp='1976-05-22 08:37:01.357000+00:00',start_station='S12',end_station='S17',start_orig_no=2856,end_orig_no=3161,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.49.csv.gz',start_timestamp ='1976-05-24 13:05:38.554000+00:00',end_timestamp='1976-05-24 15:16:28.596000+00:00',start_station='S12',end_station='S17',start_orig_no=622,end_orig_no=1715,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.49.csv.gz',start_timestamp ='1976-05-24 13:04:59.424000+00:00',end_timestamp='1976-05-24 15:16:28.596000+00:00',start_station='S12',end_station='S17',start_orig_no=627,end_orig_no=1715,orig_ground_station1=4)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.5.csv.gz',start_timestamp ='1976-04-27 21:23:52.934000+00:00',end_timestamp='1976-04-27 22:26:50.258000+00:00',start_station='S12',end_station='S17',start_orig_no=287,end_orig_no=702,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.50.csv.gz',start_timestamp ='1976-05-24 19:43:21.165000+00:00',end_timestamp='1976-05-24 19:43:56.578000+00:00',start_station='S12',end_station='S17',start_orig_no=3926,end_orig_no=3930,orig_ground_station1=4)
    
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.54.csv.gz',start_timestamp ='1976-05-26 15:54:50.936000+00:00',end_timestamp='1976-05-26 15:54:53.117000+00:00',start_station='S12',end_station='S17',start_orig_no=280,end_orig_no=280,orig_ground_station1=1)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.54.csv.gz',start_timestamp ='1976-05-26 18:21:09.412000+00:00',end_timestamp='1976-05-27 00:24:54.693000+00:00',start_station='S12',end_station='S17',start_orig_no=1260,end_orig_no=4302,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.54.csv.gz',start_timestamp ='1976-05-26 21:36:07.017000+00:00',end_timestamp='1976-05-27 00:24:54.693000+00:00',start_station='S12',end_station='S17',start_orig_no=2905,end_orig_no=4302,orig_ground_station1=4)
    
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.3.55.csv.gz',start_timestamp ='1976-05-27 03:43:15.943000+00:00',end_timestamp='1976-05-27 04:58:59.901000+00:00',start_station='S12',end_station='S17',start_orig_no=1295,end_orig_no=2200,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.11.csv.gz',start_timestamp ='1976-06-02 03:34:01.878000+00:00',end_timestamp='1976-06-02 05:14:58.447000+00:00',start_station='S12',end_station='S17',start_orig_no=1791,end_orig_no=2626,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.19.csv.gz',start_timestamp ='1976-06-06 06:43:30.360000+00:00',end_timestamp='1976-06-06 06:53:52.676000+00:00',start_station='S12',end_station='S17',start_orig_no=3198,end_orig_no=3283,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.22.csv.gz',start_timestamp ='1976-06-07 15:09:59.695000+00:00',end_timestamp='1976-06-07 15:29:54.291000+00:00',start_station='S12',end_station='S17',start_orig_no=3926,end_orig_no=4090,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.23.csv.gz',start_timestamp ='1976-06-07 20:32:51.494000+00:00',end_timestamp='1976-06-07 23:45:55.395000+00:00',start_station='S12',end_station='S17',start_orig_no=2871,end_orig_no=3792,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.23.csv.gz',start_timestamp ='1976-06-07 20:32:51.494000+00:00',end_timestamp='1976-06-07 23:45:55.395000+00:00',start_station='S12',end_station='S17',start_orig_no=318,end_orig_no=3792,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.27.csv.gz',start_timestamp ='1976-06-09 05:43:34.362000+00:00',end_timestamp='1976-06-09 06:15:18.798000+00:00',start_station='S12',end_station='S17',start_orig_no=887,end_orig_no=1095,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.3.csv.gz',start_timestamp ='1976-05-29 00:01:06.069000+00:00',end_timestamp='1976-05-29 01:53:52.668000+00:00',start_station='S12',end_station='S17',start_orig_no=3099,end_orig_no=4032,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.31.csv.gz',start_timestamp ='1976-06-10 08:13:47.155000+00:00',end_timestamp='1976-06-10 10:59:05.122000+00:00',start_station='S12',end_station='S17',start_orig_no=1654,end_orig_no=3022,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.32.csv.gz',start_timestamp ='1976-06-10 18:22:12.202000+00:00',end_timestamp='1976-06-10 18:26:54.230000+00:00',start_station='S12',end_station='S17',start_orig_no=4009,end_orig_no=4047,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.33.csv.gz',start_timestamp ='1976-06-11 09:26:03.055000+00:00',end_timestamp='1976-06-11 11:27:53.039000+00:00',start_station='S12',end_station='S17',start_orig_no=3058,end_orig_no=4066,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.43.csv.gz',start_timestamp ='1976-06-17 19:25:35.711000+00:00',end_timestamp='1976-06-17 21:52:53.953000+00:00',start_station='S12',end_station='S17',start_orig_no=2478,end_orig_no=3697,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.44.csv.gz',start_timestamp ='1976-06-18 17:51:15.670000+00:00',end_timestamp='1976-06-18 20:01:46.822000+00:00',start_station='S12',end_station='S17',start_orig_no=1916,end_orig_no=2048,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.44.csv.gz',start_timestamp ='1976-06-18 20:01:45.185000+00:00',end_timestamp='1976-06-18 20:01:46.822000+00:00',start_station='S12',end_station='S17',start_orig_no=2048,end_orig_no=2048,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.44.csv.gz',start_timestamp ='1976-06-19 00:25:23.612000+00:00',end_timestamp='1976-06-19 00:25:37.468000+00:00',start_station='S12',end_station='S17',start_orig_no=2199,end_orig_no=2200,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.44.csv.gz',start_timestamp ='1976-06-19 00:25:23.612000+00:00',end_timestamp='1976-06-19 00:25:37.468000+00:00',start_station='S12',end_station='S17',start_orig_no=2199,end_orig_no=2200,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.48.csv.gz',start_timestamp ='1976-06-20 13:00:28.689000+00:00',end_timestamp='1976-06-20 16:54:55.628000+00:00',start_station='S12',end_station='S17',start_orig_no=2861,end_orig_no=3987,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.48.csv.gz',start_timestamp ='1976-06-20 13:00:28.689000+00:00',end_timestamp='1976-06-20 16:54:55.628000+00:00',start_station='S12',end_station='S17',start_orig_no=2080,end_orig_no=3987,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.48.csv.gz',start_timestamp ='1976-06-20 16:36:49.323000+00:00',end_timestamp='1976-06-20 16:54:55.628000+00:00',start_station='S12',end_station='S17',start_orig_no=3838,end_orig_no=3987,orig_ground_station1=8)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.50.csv.gz',start_timestamp ='1976-06-21 05:18:01.766000+00:00',end_timestamp='1976-06-21 08:48:53.882000+00:00',start_station='S12',end_station='S17',start_orig_no=2867,end_orig_no=3210,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.51.csv.gz',start_timestamp ='1976-06-23 00:46:59.568000+00:00',end_timestamp='1976-06-23 02:59:55.839000+00:00',start_station='S12',end_station='S14',start_orig_no=3098,end_orig_no=4198,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.52.csv.gz',start_timestamp ='1976-06-23 06:57:17.469000+00:00',end_timestamp='1976-06-23 10:17:58.421000+00:00',start_station='S12',end_station='S17',start_orig_no=2037,end_orig_no=3698,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.53.csv.gz',start_timestamp ='1976-06-23 22:27:19.421000+00:00',end_timestamp='1976-06-23 22:54:53.302000+00:00',start_station='S12',end_station='S17',start_orig_no=2872,end_orig_no=2427,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.54.csv.gz',start_timestamp ='1976-06-24 13:42:54.979000+00:00',end_timestamp='1976-06-24 15:17:49.178000+00:00',start_station='S12',end_station='S17',start_orig_no=1415,end_orig_no=2200,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.4.7.csv.gz',start_timestamp ='1976-05-31 02:55:48.589000+00:00',end_timestamp='1976-05-31 03:49:54.022000+00:00',start_station='S12',end_station='S17',start_orig_no=2974,end_orig_no=3421,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.1.csv.gz',start_timestamp ='1976-06-25 08:35:50.278000+00:00',end_timestamp='1976-06-25 10:29:58.446000+00:00',start_station='S12',end_station='S17',start_orig_no=2854,end_orig_no=3797,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.10.csv.gz',start_timestamp ='1976-06-30 01:28:42.912000+00:00',end_timestamp='1976-06-30 02:57:56.019000+00:00',start_station='S12',end_station='S17',start_orig_no=3083,end_orig_no=3821,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.13.csv.gz',start_timestamp ='1976-07-01 10:45:08.802000+00:00',end_timestamp='1976-07-01 11:28:58.259000+00:00',start_station='S12',end_station='S15',start_orig_no=3148,end_orig_no=3510,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.16.csv.gz',start_timestamp ='1976-07-02 17:08:13.313000+00:00',end_timestamp='1976-07-02 18:30:55.642000+00:00',start_station='S12',end_station='S17',start_orig_no=2940,end_orig_no=3624,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.18.csv.gz',start_timestamp ='1976-07-04 02:19:01.111000+00:00',end_timestamp='1976-07-04 02:43:52.709000+00:00',start_station='S12',end_station='S17',start_orig_no=3147,end_orig_no=3352,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.20.csv.gz',start_timestamp ='1976-07-04 19:03:37.433000+00:00',end_timestamp='1976-07-04 20:36:57.590000+00:00',start_station='S12',end_station='S17',start_orig_no=3541,end_orig_no=4313,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.21.csv.gz',start_timestamp ='1976-07-05 02:53:16.743000+00:00',end_timestamp='1976-07-05 03:10:44.849000+00:00',start_station='S12',end_station='S12',start_orig_no=936,end_orig_no=1080,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.25.csv.gz',start_timestamp ='1976-07-07 00:51:23.780000+00:00',end_timestamp='1976-07-07 02:44:55.956000+00:00',start_station='S12',end_station='S17',start_orig_no=960,end_orig_no=1898,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.27.csv.gz',start_timestamp ='1976-07-07 21:23:12.496000+00:00',end_timestamp='1976-07-07 23:29:54.492000+00:00',start_station='S12',end_station='S17',start_orig_no=1157,end_orig_no=1711,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.30.csv.gz',start_timestamp ='1976-07-08 19:53:55.639000+00:00',end_timestamp='1976-07-08 19:54:52.791000+00:00',start_station='S12',end_station='S17',start_orig_no=4025,end_orig_no=4032,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.35.csv.gz',start_timestamp ='1976-07-10 10:09:09.319000+00:00',end_timestamp='1976-07-10 11:10:52.033000+00:00',start_station='S12',end_station='S17',start_orig_no=1024,end_orig_no=1533,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.36.csv.gz',start_timestamp ='1976-07-11 05:31:09.623000+00:00',end_timestamp='1976-07-11 05:54:58.738000+00:00',start_station='S12',end_station='S17',start_orig_no=1834,end_orig_no=2029,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.40.csv.gz',start_timestamp ='1976-07-12 18:36:13.078000+00:00',end_timestamp='1976-07-12 20:14:51.838000+00:00',start_station='S12',end_station='S17',start_orig_no=3139,end_orig_no=3955,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.49.csv.gz',start_timestamp ='1976-07-15 12:58:58.672000+00:00',end_timestamp='1976-07-15 13:43:52.852000+00:00',start_station='S12',end_station='S17',start_orig_no=3356,end_orig_no=3727,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.50.csv.gz',start_timestamp ='1976-07-15 23:33:50.835000+00:00',end_timestamp='1976-07-16 03:39:55.359000+00:00',start_station='S12',end_station='S17',start_orig_no=1538,end_orig_no=3574,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.53.csv.gz',start_timestamp ='1976-07-17 13:46:48.422000+00:00',end_timestamp='1976-07-17 13:54:53.195000+00:00',start_station='S12',end_station='S17',start_orig_no=3167,end_orig_no=3233,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.55.csv.gz',start_timestamp ='1976-07-18 05:29:43.601000+00:00',end_timestamp='1976-07-18 05:29:45.067000+00:00',start_station='S12',end_station='S17',start_orig_no=16644,end_orig_no=16644,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.55.csv.gz',start_timestamp ='1976-07-18 12:18:34.246000+00:00',end_timestamp='1976-07-18 13:49:56.824000+00:00',start_station='S12',end_station='S17',start_orig_no=3421,end_orig_no=4187,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.55.csv.gz',start_timestamp ='1976-07-18 12:18:36.058000+00:00',end_timestamp='1976-07-18 13:49:56.824000+00:00',start_station='S12',end_station='S17',start_orig_no=33892,end_orig_no=4187,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.55.csv.gz',start_timestamp ='1976-07-18 13:15:32.113000+00:00',end_timestamp='1976-07-18 13:49:56.824000+00:00',start_station='S12',end_station='S17',start_orig_no=3903,end_orig_no=4187,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.6.csv.gz',start_timestamp ='1976-06-27 22:46:48.293000+00:00',end_timestamp='1976-06-28 02:29:54.294000+00:00',start_station='S12',end_station='S17',start_orig_no=1642,end_orig_no=3488,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.6.csv.gz',start_timestamp ='1976-06-28 06:50:32.241000+00:00',end_timestamp='1976-06-28 06:50:34.091000+00:00',start_station='S12',end_station='S17',start_orig_no=2866,end_orig_no=2866,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.5.9.csv.gz',start_timestamp ='1976-06-29 02:35:45.827000+00:00',end_timestamp='1976-06-29 02:55:56.090000+00:00',start_station='S12',end_station='S17',start_orig_no=3074,end_orig_no=3241,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.14.csv.gz',start_timestamp ='1976-07-24 20:45:46.371000+00:00',end_timestamp='1976-07-24 22:55:56.272000+00:00',start_station='S12',end_station='S17',start_orig_no=3125,end_orig_no=4202,orig_ground_station1=4)
    
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.16.csv.gz',start_timestamp ='1976-07-25 23:01:59.187000+00:00',end_timestamp='1976-07-26 00:50:24.205000+00:00',start_station='S12',end_station='S14',start_orig_no=2540,end_orig_no=3436,orig_ground_station1=5)
    
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.17.csv.gz',start_timestamp ='1976-07-26 05:19:59.572000+00:00',end_timestamp='1976-07-26 07:24:57.634000+00:00',start_station='S12',end_station='S17',start_orig_no=2206,end_orig_no=3240,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.19.csv.gz',start_timestamp ='1976-07-27 05:36:06.668000+00:00',end_timestamp='1976-07-27 06:43:53.006000+00:00',start_station='S12',end_station='S17',start_orig_no=2876,end_orig_no=4256,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.2.csv.gz',start_timestamp ='1976-07-19 05:39:50.229000+00:00',end_timestamp='1976-07-19 05:39:52.416000+00:00',start_station='S12',end_station='S17',start_orig_no=2870,end_orig_no=2870,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.22.csv.gz',start_timestamp ='1976-07-28 11:34:48.950000+00:00',end_timestamp='1976-07-28 13:39:55.284000+00:00',start_station='S12',end_station='S17',start_orig_no=26511,end_orig_no=2752,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.23.csv.gz',start_timestamp ='1976-07-29 09:24:46.633000+00:00',end_timestamp='1976-07-29 09:24:48.557000+00:00',start_station='S12',end_station='S17',start_orig_no=2875,end_orig_no=2875,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.25.csv.gz',start_timestamp ='1976-07-30 00:46:00.922000+00:00',end_timestamp='1976-07-30 02:54:58.157000+00:00',start_station='S12',end_station='S17',start_orig_no=3091,end_orig_no=4158,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.27.csv.gz',start_timestamp ='1976-07-30 18:00:28.788000+00:00',end_timestamp='1976-07-30 19:24:55.024000+00:00',start_station='S12',end_station='S17',start_orig_no=2856,end_orig_no=2898,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.29.csv.gz',start_timestamp ='1976-07-31 15:27:08.591000+00:00',end_timestamp='1976-07-31 15:27:10.526000+00:00',start_station='S12',end_station='S17',start_orig_no=2867,end_orig_no=2867,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.29.csv.gz',start_timestamp ='1976-07-31 21:25:47.980000+00:00',end_timestamp='1976-07-31 21:26:18.743000+00:00',start_station='S12',end_station='S17',start_orig_no=2038,end_orig_no=2860,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.29.csv.gz',start_timestamp ='1976-07-31 21:26:16.960000+00:00',end_timestamp='1976-07-31 21:26:18.743000+00:00',start_station='S12',end_station='S17',start_orig_no=2860,end_orig_no=2860,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.3.csv.gz',start_timestamp ='1976-07-19 14:04:39.308000+00:00',end_timestamp='1976-07-19 14:44:58.707000+00:00',start_station='S12',end_station='S17',start_orig_no=2187,end_orig_no=2520,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.33.csv.gz',start_timestamp ='1976-08-02 10:06:58.984000+00:00',end_timestamp='1976-08-02 11:55:53.348000+00:00',start_station='S12',end_station='S17',start_orig_no=3131,end_orig_no=4032,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.34.csv.gz',start_timestamp ='1976-08-02 22:51:48.645000+00:00',end_timestamp='1976-08-02 23:55:56.224000+00:00',start_station='S12',end_station='S14',start_orig_no=8310,end_orig_no=1549,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.35.csv.gz',start_timestamp ='1976-08-03 18:25:53.315000+00:00',end_timestamp='1976-08-03 18:28:54.378000+00:00',start_station='S12',end_station='S17',start_orig_no=2170,end_orig_no=2194,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.35.csv.gz',start_timestamp ='1976-08-03 18:25:53.315000+00:00',end_timestamp='1976-08-03 18:28:54.378000+00:00',start_station='S12',end_station='S17',start_orig_no=2170,end_orig_no=2194,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.4.csv.gz',start_timestamp ='1976-07-20 07:36:02.795000+00:00',end_timestamp='1976-07-20 08:56:59.316000+00:00',start_station='S12',end_station='S17',start_orig_no=2272,end_orig_no=2940,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.45.csv.gz',start_timestamp ='1976-08-07 16:36:44.875000+00:00',end_timestamp='1976-08-07 17:45:55.367000+00:00',start_station='S12',end_station='S17',start_orig_no=1873,end_orig_no=2445,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.46.csv.gz',start_timestamp ='1976-08-08 14:24:59.328000+00:00',end_timestamp='1976-08-08 16:24:54.113000+00:00',start_station='S12',end_station='S17',start_orig_no=257,end_orig_no=1249,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.52.csv.gz',start_timestamp ='1976-08-11 09:29:21.632000+00:00',end_timestamp='1976-08-11 10:54:53.024000+00:00',start_station='S12',end_station='S17',start_orig_no=826,end_orig_no=1532,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.6.9.csv.gz',start_timestamp ='1976-07-22 22:13:30.683000+00:00',end_timestamp='1976-07-22 22:24:58.805000+00:00',start_station='S12',end_station='S17',start_orig_no=3184,end_orig_no=3278,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.11.csv.gz',start_timestamp ='1976-08-17 09:01:31.871000+00:00',end_timestamp='1976-08-17 13:44:57.419000+00:00',start_station='S12',end_station='S14',start_orig_no=1940,end_orig_no=3823,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.11.csv.gz',start_timestamp ='1976-08-17 10:35:59.319000+00:00',end_timestamp='1976-08-17 13:44:57.419000+00:00',start_station='S12',end_station='S14',start_orig_no=2572,end_orig_no=3823,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.13.csv.gz',start_timestamp ='1976-08-18 12:45:21.767000+00:00',end_timestamp='1976-08-18 13:24:54.015000+00:00',start_station='S12',end_station='S14',start_orig_no=2526,end_orig_no=2787,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.17.csv.gz',start_timestamp ='1976-08-20 12:43:35.589000+00:00',end_timestamp='1976-08-20 13:44:52.042000+00:00',start_station='S12',end_station='S14',start_orig_no=3397,end_orig_no=3802,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.19.csv.gz',start_timestamp ='1976-08-22 01:22:47.032000+00:00',end_timestamp='1976-08-22 02:43:52.799000+00:00',start_station='S12',end_station='S14',start_orig_no=2202,end_orig_no=2737,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.20.csv.gz',start_timestamp ='1976-08-22 06:11:10.441000+00:00',end_timestamp='1976-08-22 06:59:57.508000+00:00',start_station='S12',end_station='S14',start_orig_no=63104,end_orig_no=659,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.29.csv.gz',start_timestamp ='1976-08-26 22:28:41.857000+00:00',end_timestamp='1976-08-26 23:54:52.230000+00:00',start_station='S12',end_station='S14',start_orig_no=2315,end_orig_no=2885,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.30.csv.gz',start_timestamp ='1976-08-27 17:02:30.861000+00:00',end_timestamp='1976-08-27 18:17:58.314000+00:00',start_station='S12',end_station='S14',start_orig_no=774,end_orig_no=1273,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.33.csv.gz',start_timestamp ='1976-08-29 18:03:11.720000+00:00',end_timestamp='1976-08-29 18:31:15.379000+00:00',start_station='S12',end_station='S14',start_orig_no=2472,end_orig_no=2657,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.34.csv.gz',start_timestamp ='1976-08-29 18:31:23.519000+00:00',end_timestamp='1976-08-29 22:23:09.896000+00:00',start_station='S12',end_station='S14',start_orig_no=28,end_orig_no=4194,orig_ground_station1=1)
    # not required
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.40.csv.gz',start_timestamp ='1976-09-02 18:35:14.484000+00:00',end_timestamp='1976-09-02 20:57:51.067000+00:00',start_station='S12',end_station='S14',start_orig_no=635,end_orig_no=2689,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.40.csv.gz',start_timestamp ='1976-09-03 05:28:46.562000+00:00',end_timestamp='1976-09-03 07:06:52.694000+00:00',start_station='S12',end_station='S14',start_orig_no=2546,end_orig_no=3195,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.41.csv.gz',start_timestamp ='1976-09-03 22:04:17.398000+00:00',end_timestamp='1976-09-03 22:31:54.559000+00:00',start_station='S12',end_station='S14',start_orig_no=2008,end_orig_no=2190,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.45.csv.gz',start_timestamp ='1976-09-05 16:53:15.470000+00:00',end_timestamp='1976-09-05 16:53:17.838000+00:00',start_station='S12',end_station='S14',start_orig_no=713,end_orig_no=15,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.45.csv.gz',start_timestamp ='1976-09-05 16:53:17.281000+00:00',end_timestamp='1976-09-05 16:53:17.838000+00:00',start_station='S12',end_station='S14',start_orig_no=15,end_orig_no=15,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.49.csv.gz',start_timestamp ='1976-09-07 01:51:42.430000+00:00',end_timestamp='1976-09-07 03:24:58.427000+00:00',start_station='S12',end_station='S14',start_orig_no=2555,end_orig_no=3172,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.56.csv.gz',start_timestamp ='1976-09-11 17:09:45.104000+00:00',end_timestamp='1976-09-11 17:56:51.644000+00:00',start_station='S12',end_station='S14',start_orig_no=50305,end_orig_no=3137,orig_ground_station1=6)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.7.56.csv.gz',start_timestamp ='1976-09-12 05:23:02.484000+00:00',end_timestamp='1976-09-12 07:58:56.200000+00:00',start_station='S12',end_station='S14',start_orig_no=808,end_orig_no=1840,orig_ground_station1=1)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.17.csv.gz',start_timestamp ='1976-09-21 02:00:49.480000+00:00',end_timestamp='1976-09-21 05:28:57.428000+00:00',start_station='S12',end_station='S14',start_orig_no=106,end_orig_no=1484,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.19.csv.gz',start_timestamp ='1976-09-22 02:31:59.478000+00:00',end_timestamp='1976-09-22 06:07:58.208000+00:00',start_station='S12',end_station='S14',start_orig_no=246,end_orig_no=1676,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.21.csv.gz',start_timestamp ='1976-09-23 03:19:40.930000+00:00',end_timestamp='1976-09-23 03:19:51.508000+00:00',start_station='S12',end_station='S14',start_orig_no=1572,end_orig_no=2876,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.24.csv.gz',start_timestamp ='1976-09-24 16:10:59.475000+00:00',end_timestamp='1976-09-24 17:43:57.178000+00:00',start_station='S12',end_station='S14',start_orig_no=1962,end_orig_no=2577,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.26.csv.gz',start_timestamp ='1976-09-25 17:07:20.868000+00:00',end_timestamp='1976-09-25 18:29:50.965000+00:00',start_station='S12',end_station='S14',start_orig_no=2238,end_orig_no=2784,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.27.csv.gz',start_timestamp ='1976-09-26 05:24:36.860000+00:00',end_timestamp='1976-09-26 09:57:57.172000+00:00',start_station='S12',end_station='S14',start_orig_no=2425,end_orig_no=4235,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.29.csv.gz',start_timestamp ='1976-09-27 07:47:59.816000+00:00',end_timestamp='1976-09-27 10:55:53.660000+00:00',start_station='S12',end_station='S14',start_orig_no=2386,end_orig_no=3630,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.31.csv.gz',start_timestamp ='1976-09-28 07:33:59.006000+00:00',end_timestamp='1976-09-28 10:16:50.045000+00:00',start_station='S12',end_station='S14',start_orig_no=2425,end_orig_no=3503,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.37.csv.gz',start_timestamp ='1976-10-01 18:07:20.464000+00:00',end_timestamp='1976-10-01 20:29:52.036000+00:00',start_station='S12',end_station='S14',start_orig_no=3087,end_orig_no=4029,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.37.csv.gz',start_timestamp ='1976-10-02 02:25:48.715000+00:00',end_timestamp='1976-10-02 04:59:54.474000+00:00',start_station='S12',end_station='S14',start_orig_no=2401,end_orig_no=3421,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.4.csv.gz',start_timestamp ='1976-09-14 04:47:15.064000+00:00',end_timestamp='1976-09-14 09:21:50.299000+00:00',start_station='S12',end_station='S14',start_orig_no=1778,end_orig_no=3595,orig_ground_station1=7)
    # manually changed to correct record
    # 2352,"1976-10-07T02:31:31.812000Z",10,"S15",0,7,516,516,429,516,516,429,516,516,429,515,516,429,"00011",47,"1110001001000011101101"

    # new!!
    df = reset_all_ground_stations_idx(df, gzip_filename, problem_gzip_filename='wtn.8.46.csv.gz',orig_idx_start =70512,orig_idx_end=109616,single_station=None)

    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.43.csv.gz',start_timestamp ='1976-10-05 02:37:50.565000+00:00',end_timestamp='1976-10-05 04:52:55.104000+00:00',start_station='S15',end_station='S14',start_orig_no=2352,end_orig_no=3274,orig_ground_station1=2)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.46.csv.gz',start_timestamp ='1976-10-07T02:31:31.812000Z',end_timestamp='1976-10-07 04:09:54.967000+00:00',start_station='S12',end_station='S14',start_orig_no=2352,end_orig_no=3003,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.47.csv.gz',start_timestamp ='1976-10-07 16:53:23.203000+00:00',end_timestamp='1976-10-07 18:42:54.475000+00:00',start_station='S12',end_station='S12',start_orig_no=2201,end_orig_no=2926,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.5.csv.gz',start_timestamp ='1976-09-14 22:15:59.429000+00:00',end_timestamp='1976-09-15 00:14:54.651000+00:00',start_station='S12',end_station='S14',start_orig_no=2537,end_orig_no=3324,orig_ground_station1=5)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.50.csv.gz',start_timestamp ='1976-10-09 03:25:28.538000+00:00',end_timestamp='1976-10-09 03:52:56.271000+00:00',start_station='S12',end_station='S14',start_orig_no=2572,end_orig_no=2753,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.50.csv.gz',start_timestamp ='1976-10-09 03:25:31.342000+00:00',end_timestamp='1976-10-09 03:52:56.271000+00:00',start_station='S12',end_station='S14',start_orig_no=2572,end_orig_no=2753,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.50.csv.gz',start_timestamp ='1976-10-09 02:29:58.827000+00:00',end_timestamp='1976-10-09 03:52:56.271000+00:00',start_station='S12',end_station='S14',start_orig_no=2201,end_orig_no=2753,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.50.csv.gz',start_timestamp ='1976-10-09 03:25:28.538000+00:00',end_timestamp='1976-10-09 03:52:56.271000+00:00',start_station='S12',end_station='S14',start_orig_no=2572,end_orig_no=2753,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.50.csv.gz',start_timestamp ='1976-10-09 03:25:31.342000+00:00',end_timestamp='1976-10-09 03:52:56.271000+00:00',start_station='S12',end_station='S14',start_orig_no=2572,end_orig_no=2753,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.53.csv.gz',start_timestamp ='1976-10-10 02:55:36.888000+00:00',end_timestamp='1976-10-10 03:00:55.693000+00:00',start_station='S12',end_station='S16',start_orig_no=2301,end_orig_no=2333,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.54.csv.gz',start_timestamp ='1976-10-10 22:58:41.444000+00:00',end_timestamp='1976-10-10 22:58:42.434000+00:00',start_station='S12',end_station='S16',start_orig_no=2850,end_orig_no=2850,orig_ground_station1=7)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.54.csv.gz',start_timestamp ='1976-10-11 06:39:59.511000+00:00',end_timestamp='1976-10-11 08:39:55.110000+00:00',start_station='S12',end_station='S16',start_orig_no=373,end_orig_no=968,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.55.csv.gz',start_timestamp ='1976-10-12 12:09:24.907000+00:00',end_timestamp='1976-10-12 14:07:56.550000+00:00',start_station='S12',end_station='S16',start_orig_no=2974,end_orig_no=3562,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.6.csv.gz',start_timestamp ='1976-09-15 08:08:19.094000+00:00',end_timestamp='1976-09-15 09:23:55.455000+00:00',start_station='S12',end_station='S14',start_orig_no=3015,end_orig_no=3515,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.8.7.csv.gz',start_timestamp ='1976-09-16 06:28:59.392000+00:00',end_timestamp='1976-09-16 09:32:58.273000+00:00',start_station='S12',end_station='S14',start_orig_no=3028,end_orig_no=4246,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.1.csv.gz',start_timestamp ='1976-10-12 19:44:21.863000+00:00',end_timestamp='1976-10-12 21:29:48.628000+00:00',start_station='S12',end_station='S16',start_orig_no=1415,end_orig_no=1938,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.10.csv.gz',start_timestamp ='1976-10-20 16:30:54.196000+00:00',end_timestamp='1976-10-20 19:29:48.661000+00:00',start_station='S12',end_station='S16',start_orig_no=2191,end_orig_no=3079,orig_ground_station1=3)
    # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.13.csv.gz',start_timestamp ='1976-10-22 15:05:11.055000+00:00',end_timestamp='1976-10-22 16:59:56.890000+00:00',start_station='S12',end_station='S16',start_orig_no=3473,end_orig_no=4041,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.14.csv.gz',start_timestamp ='1976-10-22 18:47:52.802000+00:00',end_timestamp='1976-10-22 19:29:48.612000+00:00',start_station='S12',end_station='S16',start_orig_no=1,end_orig_no=759,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.14.csv.gz',start_timestamp ='1976-10-23 05:01:38.399000+00:00',end_timestamp='1976-10-23 07:19:49.741000+00:00',start_station='S12',end_station='S16',start_orig_no=1476,end_orig_no=2172,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.16.csv.gz',start_timestamp ='1976-10-24 20:28:46.590000+00:00',end_timestamp='1976-10-24 22:59:55.824000+00:00',start_station='S12',end_station='S16',start_orig_no=2946,end_orig_no=3696,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.17.csv.gz',start_timestamp ='1976-10-25 16:55:47.657000+00:00',end_timestamp='1976-10-25 19:29:56.782000+00:00',start_station='S12',end_station='S16',start_orig_no=975,end_orig_no=1740,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.20.csv.gz',start_timestamp ='1976-10-27 09:15:48.444000+00:00',end_timestamp='1976-10-27 10:58:49.928000+00:00',start_station='S12',end_station='S16',start_orig_no=2231,end_orig_no=2742,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.25.csv.gz',start_timestamp ='1976-11-01 03:58:49.397000+00:00',end_timestamp='1976-11-01 05:29:46.865000+00:00',start_station='S12',end_station='S16',start_orig_no=1852,end_orig_no=2303,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.25.csv.gz',start_timestamp ='1976-11-01 10:00:18.481000+00:00',end_timestamp='1976-11-01 11:19:47.329000+00:00',start_station='S12',end_station='S16',start_orig_no=3758,end_orig_no=4152,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.26.csv.gz',start_timestamp ='1976-11-01 19:06:52.305000+00:00',end_timestamp='1976-11-01 21:02:58.512000+00:00',start_station='S12',end_station='S16',start_orig_no=2201,end_orig_no=2777,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.29.csv.gz',start_timestamp ='1976-11-03 09:05:49.233000+00:00',end_timestamp='1976-11-03 09:44:50.739000+00:00',start_station='S12',end_station='S16',start_orig_no=2630,end_orig_no=2823,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.34.csv.gz',start_timestamp ='1976-11-07 01:21:15.496000+00:00',end_timestamp='1976-11-07 05:00:50.292000+00:00',start_station='S12',end_station='S16',start_orig_no=2,end_orig_no=1105,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.34.csv.gz',start_timestamp ='1976-11-07 01:33:21.815000+00:00',end_timestamp='1976-11-07 05:00:50.292000+00:00',start_station='S12',end_station='S16',start_orig_no=75,end_orig_no=1105,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.34.csv.gz',start_timestamp ='1976-11-07 10:32:51.320000+00:00',end_timestamp='1976-11-07 12:12:15.569000+00:00',start_station='S12',end_station='S16',start_orig_no=1707,end_orig_no=2200,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.35.csv.gz',start_timestamp ='1976-11-08 04:27:20.436000+00:00',end_timestamp='1976-11-08 07:47:58.882000+00:00',start_station='S12',end_station='S16',start_orig_no=582,end_orig_no=1578,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.37.csv.gz',start_timestamp ='1976-11-09 03:38:04.151000+00:00',end_timestamp='1976-11-09 04:10:51.648000+00:00',start_station='S12',end_station='S16',start_orig_no=2038,end_orig_no=2200,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.38.csv.gz',start_timestamp ='1976-11-10 05:32:14.672000+00:00',end_timestamp='1976-11-10 08:59:55.519000+00:00',start_station='S12',end_station='S16',start_orig_no=2245,end_orig_no=3276,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.39.csv.gz',start_timestamp ='1976-11-11 04:24:58.879000+00:00',end_timestamp='1976-11-11 08:14:47.828000+00:00',start_station='S12',end_station='S16',start_orig_no=2275,end_orig_no=3416,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.40.csv.gz',start_timestamp ='1976-11-11 11:03:58.156000+00:00',end_timestamp='1976-11-11 15:01:58.046000+00:00',start_station='S12',end_station='S16',start_orig_no=1,end_orig_no=1990,orig_ground_station1=4)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.44.csv.gz',start_timestamp ='1976-11-13 05:00:51.408000+00:00',end_timestamp='1976-11-13 08:29:53.156000+00:00',start_station='S12',end_station='S14',start_orig_no=2757,end_orig_no=4141,orig_ground_station1=7)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.46.csv.gz',start_timestamp ='1976-11-14 06:28:44.059000+00:00',end_timestamp='1976-11-14 08:29:55.566000+00:00',start_station='S12',end_station='S14',start_orig_no=2563,end_orig_no=3365,orig_ground_station1=7)
    # XFAIL
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.50.csv.gz',start_timestamp ='1976-11-16 14:12:38.040000+00:00',end_timestamp='1976-11-16 14:49:54.908000+00:00',start_station='S12',end_station='S14',start_orig_no=2511,end_orig_no=2757,orig_ground_station1=8)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.8.csv.gz',start_timestamp ='1976-10-19 01:46:47.100000+00:00',end_timestamp='1976-10-19 03:14:55.676000+00:00',start_station='S12',end_station='S16',start_orig_no=1554,end_orig_no=1991,orig_ground_station1=9)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.9.9.csv.gz',start_timestamp ='1976-10-20 02:11:27.685000+00:00',end_timestamp='1976-10-20 04:15:55.892000+00:00',start_station='S12',end_station='S16',start_orig_no=1629,end_orig_no=1980,orig_ground_station1=9)




    # main tapes
    df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='pse.a12.1.12.csv.gz', 
              orig_idx_start=1, orig_idx_end=1377,single_station='S12')
    rec_drop_damaged_total += rec_drop_damaged1
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.103.csv.gz',start_timestamp ='1970-02-27 15:44:50.445000+00:00',end_timestamp='1970-02-27 16:16:30.466000+00:00',start_station='S12',end_station='S12',start_orig_no=32,end_orig_no=43,orig_ground_station1=10)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.15.csv.gz',start_timestamp ='1969-12-04 02:15:04.559000+00:00',end_timestamp='1969-12-04 03:37:29.890000+00:00',start_station='S12',end_station='S12',start_orig_no=264,end_orig_no=295,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.21.csv.gz',start_timestamp ='1969-12-09 17:46:19.581000+00:00',end_timestamp='1969-12-09 23:09:14.485000+00:00',start_station='S12',end_station='S12',start_orig_no=77,end_orig_no=196,orig_ground_station1=11)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.21.csv.gz',start_timestamp ='1969-12-10 06:49:33.560000+00:00',end_timestamp='1969-12-10 08:00:37.844000+00:00',start_station='S12',end_station='S12',start_orig_no=366,end_orig_no=392,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.27.csv.gz',start_timestamp ='1969-12-15 20:09:40.256000+00:00',end_timestamp='1969-12-15 21:01:50.696000+00:00',start_station='S12',end_station='S12',start_orig_no=129,end_orig_no=149,orig_ground_station1=2)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.35.csv.gz',start_timestamp ='1969-12-23 20:33:22.454000+00:00',end_timestamp='1969-12-23 21:04:42.576000+00:00',start_station='S12',end_station='S12',start_orig_no=138,end_orig_no=149,orig_ground_station1=3)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.37.csv.gz',start_timestamp ='1969-12-25 20:52:33.500000+00:00',end_timestamp='1969-12-25 21:30:37.543000+00:00',start_station='S12',end_station='S12',start_orig_no=145,end_orig_no=159,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.39.csv.gz',start_timestamp ='1969-12-27 22:50:55.426000+00:00',end_timestamp='1969-12-27 23:03:54.888000+00:00',start_station='S12',end_station='S12',start_orig_no=189,end_orig_no=193,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.41.csv.gz',start_timestamp ='1969-12-29 14:29:29.465000+00:00',end_timestamp='1969-12-29 16:05:59.547000+00:00',start_station='S12',end_station='S12',start_orig_no=4,end_orig_no=39,orig_ground_station1=11)
# check
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.46.csv.gz',start_timestamp ='1970-01-03 18:52:02.258000+00:00',end_timestamp='1970-01-03 20:01:53.484000+00:00',start_station='S12',end_station='S12',start_orig_no=100,end_orig_no=126,orig_ground_station1=10)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.46.csv.gz',start_timestamp ='1970-01-04 10:55:44.436000+00:00',end_timestamp='1970-01-04 11:44:29.599000+00:00',start_station='S12',end_station='S12',start_orig_no=456,end_orig_no=474,orig_ground_station1=2)

    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.5.csv.gz',start_timestamp ='1969-11-24 08:05:02.424000+00:00',end_timestamp='1969-11-24 13:20:56.401000+00:00',start_station='S12',end_station='S12',start_orig_no=394,end_orig_no=510,orig_ground_station1=11)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.68.csv.gz',start_timestamp ='1970-01-26 09:39:36.473000+00:00',end_timestamp='1970-01-26 13:36:29.699000+00:00',start_station='S12',end_station='S12',start_orig_no=422,end_orig_no=509,orig_ground_station1=11)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.77.csv.gz',start_timestamp ='1970-02-02 02:11:05.459000+00:00',end_timestamp='1970-02-02 03:36:49.401000+00:00',start_station='S12',end_station='S12',start_orig_no=262,end_orig_no=293,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.83.csv.gz',start_timestamp ='1970-02-08 08:22:36.501000+00:00',end_timestamp='1970-02-08 09:49:40.744000+00:00',start_station='S12',end_station='S12',start_orig_no=398,end_orig_no=430,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.84.csv.gz',start_timestamp ='1970-02-09 01:57:34.682000+00:00',end_timestamp='1970-02-09 10:04:54.224000+00:00',start_station='S12',end_station='S12',start_orig_no=257,end_orig_no=437,orig_ground_station1=5)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.1.9.csv.gz',start_timestamp ='1969-11-27 19:44:34.485000+00:00',end_timestamp='1969-11-27 22:56:04.749000+00:00',start_station='S12',end_station='S12',start_orig_no=120,end_orig_no=190,orig_ground_station1=5)

    df = reset_all_ground_stations_idx(df, gzip_filename, problem_gzip_filename='pse.a12.2.84.csv.gz',orig_idx_start=0,orig_idx_end=143371,single_station=None)
    df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='pse.a12.2.126.csv.gz',start_timestamp ='1970-09-15 08:33:51.466000+00:00',end_timestamp='1970-09-15 08:34:04.749000+00:00',start_station='S12',end_station='S12',start_orig_no=403,end_orig_no=403,orig_ground_station1=11)




    return df, rec_drop_damaged_total

def all_drop_idx(df, gzip_filename, problem_gzip_filename,
      orig_idx_start,orig_idx_end,single_station=None):

    rec_drop_damaged = 0
    if problem_gzip_filename in gzip_filename:

        to_drop = list(range(orig_idx_start,orig_idx_end+1))

        len_before = len(df)
        if single_station is not None:
            df.drop(df[(df.orig_idx.isin(to_drop)) & (df.orig_station == single_station)].index, inplace=True)
        else:
            df.drop(df[df.orig_idx.isin(to_drop)].index, inplace=True)
        df.reset_index(inplace=True,drop=True)
        len_after = len(df)

        rec_drop_damaged = len_before - len_after

        logging.info('WARNING: {} Removing {} known damaged record(s)'.format(
          problem_gzip_filename,rec_drop_damaged))   

        print('WARNING: {} Removing {} known damaged record(s)'.format(
          problem_gzip_filename,rec_drop_damaged))   

    return df, rec_drop_damaged

def all_drop_data(df, gzip_filename, problem_gzip_filename,
      start_timestamp,end_timestamp,start_station,end_station,corr_ground_station1,corr_ground_station2=None,start_orig_no=None,end_orig_no=None,single_station=None):

    # manually mark data for deletion
    # this method is deprecated - use all_drop_idx which is a lot easier to 
    # use

    rec_drop_damaged = 0

    if problem_gzip_filename in gzip_filename:
        if corr_ground_station2 is None:
            corr_ground_station2 = corr_ground_station1


     
        if pd.isna(start_station):
            start_station = 'S0'
        if pd.isna(end_station):
            end_station = 'S0'
        if pd.isna(corr_ground_station1):
            corr_ground_station1 = -1
        if pd.isna(corr_ground_station2):
            corr_ground_station2 = -1  

        # make a temp datafile, where the nulls are replaced
        df_not_na = df
        df_not_na['corr_ground_station'].fillna(-1,inplace=True)
        df_not_na['orig_station'].fillna('S0',inplace=True)

        start_i_error = False
        end_i_error = False
        start_i_orig_no_error = False
        end_i_orig_no_error = False
        if start_orig_no is None:
            try:
                start_i = (df_not_na[(df_not_na['orig_timestamp'] == start_timestamp) & 
                  (df_not_na['orig_station'] == start_station) & 
                  (df_not_na['corr_ground_station'] == corr_ground_station1)].index[0])
            except IndexError:
                logging.info('Error in start_i_error')
                start_i_error = True

            try:
                end_i = (df_not_na[(df_not_na['orig_timestamp'] == end_timestamp) &
                  (df_not_na['orig_station'] == end_station) &
                  (df_not_na['corr_ground_station'] == corr_ground_station2)].index[0])
            except IndexError:
                logging.info('Error in end_i_error')
                end_i_error = True
        else:
            try:
                start_i = (df_not_na[(df_not_na['orig_timestamp'] == start_timestamp) & 
                  (df_not_na['orig_station'] == start_station) & 
                  (df_not_na['corr_ground_station'] == corr_ground_station1) & 
                  (df_not_na['orig_no'] == start_orig_no)
                    ].index[0])
            except IndexError:
                logging.info('Error in start_i_orig_no_error')
                start_i_orig_no_error = True

            try:
                end_i = (df_not_na[(df_not_na['orig_timestamp'] == end_timestamp) &
                  (df_not_na['orig_station'] == end_station) &
                  (df_not_na['corr_ground_station'] == corr_ground_station2) &
                  (df_not_na['orig_no'] == end_orig_no)
                    ].index[0])
            except IndexError:
                logging.info('Error in end_i_orig_no_error')
                end_i_orig_no_error = True

        if start_i_error or end_i_error or start_i_orig_no_error or end_i_orig_no_error:
            logging.info('Exception occurred finding data in this datafile:\n{}'.format(df.iloc[0:400].to_string()))
            logging.info('EXCEPTION: IndexError {} {} {} {} {} {} {} {} {}'.format(problem_gzip_filename, start_timestamp,end_timestamp,start_orig_no,end_orig_no,start_station,end_station,corr_ground_station1,corr_ground_station2))
            raise IndexError


        first_slice = df[0:start_i]
        end_slice = df[end_i+1:]

        before = len(df) 
        if single_station is None:
            df = pd.concat([first_slice, end_slice], ignore_index = True)
        else:
            middle_slice = df[start_i:end_i+2].copy()
            middle_slice.drop(middle_slice[(middle_slice['orig_station'] == single_station)].index , inplace=True)
            middle_slice.reset_index(inplace=True,drop=True)
            a1 = pd.concat([first_slice, middle_slice], ignore_index = True)
            df = pd.concat([a1, end_slice], ignore_index = True)
            
        after = len(df)

        rec_drop_damaged = before-after
        logging.info('WARNING: {} Removing {} known damaged record(s)'.format(
          problem_gzip_filename,rec_drop_damaged))
        

    return df, rec_drop_damaged

def timestamp_adjust(df, gzip_filename, problem_gzip_filename,old_timestamp,new_timestamp,station1,corr_ground_station1):
    if problem_gzip_filename in gzip_filename:
        start_i = (df[(df['orig_timestamp'] == old_timestamp) & 
          (df['orig_station'] == station1) & 
          (df['corr_ground_station'] == corr_ground_station1)].index[0])
        df.at[start_i,'corr_timestamp'] = new_timestamp
        df.at[start_i,'clock_flag'] = 0
        logging.info('WARNING: Adjusting timestamp for known damaged timestamp\n{}'.format(
          df.iloc[start_i:start_i+1].to_string()))
    
    return df

def get_new_ground_station(ground_station1):

    if ground_station1 == -1:
        logging.info('EXCEPTION: Ground station = -1')
        raise Exception

    if ground_station1 < 99:
        low = ground_station1*100
        high = (ground_station1+1)*100
    else:
        low = ground_station1+1
        a = ground_station1 // 100
        high = (a + 1)*100
    
    # find values greater than low and smaller than high
    gs = [i for i in config.extra_ground_stations if (i >= low) and (i < high)]
    if len(gs) == 0:
        new_ground_station = low
    else: 
        new_ground_station = gs[-1] +1 
    config.extra_ground_stations.append(new_ground_station)

    
    # print('ground_station1 ', ground_station1,new_ground_station)
    return new_ground_station

def reset_all_ground_stations_idx(df, gzip_filename, problem_gzip_filename,
        orig_idx_start,orig_idx_end,single_station=None):

    if problem_gzip_filename in gzip_filename:
        logging.info('resetting some ground stations')
        print('resetting some ground stations')

        ground_station=df.iloc[orig_idx_start].orig_ground_station
        # get the new ground station
        new_ground_station = get_new_ground_station(ground_station1=ground_station)
        # new_ground_station = 44

        to_change = list(range(orig_idx_start,orig_idx_end+1))

        # logging.info(df.head().to_string())
        len_before = len(df)
        if single_station is not None:

            # df.loc[(df['1st']=='a') & (df['2nd']==2), '2nd'] = 9
            # df.drop(df[(df.orig_idx.isin(to_drop)) & (df.orig_station == single_station)].index, inplace=True)
            df.at[df.orig_idx.isin(to_change), 'corr_ground_station'] = new_ground_station
        else:
            df.at[df.orig_idx.isin(to_change), 'corr_ground_station'] = new_ground_station
        df.reset_index(inplace=True,drop=True)
        len_after = len(df)

        # logging.info(df.head(20).to_string())

        # rec_drop_damaged = len_before - len_after

    return df



def reset_all_ground_stations(df, gzip_filename, problem_gzip_filename,
      start_timestamp,end_timestamp,start_station,end_station,orig_ground_station1,start_orig_no=None,end_orig_no=None,orig_ground_station2=None,new_ground_station=None):

    if problem_gzip_filename in gzip_filename:
        orig_len = len(df)


        if start_orig_no is None or end_orig_no is None:
            # TODO sort this out 
            print('Set orig no')
            exit()

        if new_ground_station is None:
            
            # get the new ground station
            new_ground_station = get_new_ground_station(ground_station1=orig_ground_station1)

        if orig_ground_station2 is None:
            orig_ground_station2 = orig_ground_station1

        if pd.isna(start_station):
            start_station = 'S0'
        if pd.isna(end_station):
            end_station = 'S0'
        if pd.isna(orig_ground_station1):
            orig_ground_station1 = -1
        if pd.isna(orig_ground_station2):
            orig_ground_station2 = -1  

    # # make a temp datafile, where the nulls are replaced
    # df_not_na = df
    # df_not_na['corr_ground_station'].fillna(-1,inplace=True)
    # df_not_na['orig_station'].fillna('S0',inplace=True)
        # try:
        start_i = (df[(df['orig_timestamp'] == start_timestamp) & 
          (df['orig_station'] == start_station) & 
          (df['orig_ground_station'] == orig_ground_station1) &
          (df['orig_no'] == start_orig_no)
          ].index[0])
        end_i = (df[(df['orig_timestamp'] == end_timestamp) & 
          (df['orig_station'] == end_station) & 
          (df['orig_ground_station'] == orig_ground_station2) &
          (df['orig_no'] == end_orig_no)
          ].index[0])
        # except IndexError:
        #     logging.info('EXCEPTION: IndexError {} {} {} {} {} {} {}'.format(problem_gzip_filename, start_timestamp,end_timestamp,start_station,end_station,corr_ground_station1,corr_ground_station2))
        #     raise IndexError

        # logging.info(start_i)
        # logging.info(end_i)
        df_fixed = df.iloc[start_i:end_i+1].copy()
        # t1 = datetime.now()
        df_fixed.loc[:, 'corr_ground_station'] = new_ground_station
        # df_fixed.loc[(df_fixed['orig_ground_station'] == corr_ground_station1), 'orig_ground_station'] = new_ground_station

        # logging.info(df_fixed.head().to_string())
        # logging.info('{},{}'.format(start_i,end_i))
        # # 305640,355739


        first_slice = df[0:start_i]
        second_slice = df[end_i+1:]
        
        df2 = pd.concat([first_slice, df_fixed], ignore_index = True)
        df = pd.concat([df2, second_slice], ignore_index = True)

        end_len = len(df)
        if orig_len != end_len:
            logging.info('EXCEPTION: Problem with reset_all_ground_stations start_i={} end_i={}'.format(start_i,end_i,len(df_fixed)))
            raise Exception
        if len(df_fixed) == 0:
            logging.info('EXCEPTION: Problem with reset_all_ground_stations start_i={} end_i={} df_fixed={}'.format(start_i,end_i,len(df_fixed)))
            raise Exception

    return df

def all_station_order(df,gzip_filename):

    logging.debug('DEBUG: all_station_order')


    # get the station order from the last 5 records (seems to go wrong less 
    # frequently)
    station_list = df['orig_station'].iloc[-5:].to_list()
    station_dict = OrderedDict()
    # add them to the dictionary backwards
    for i, x in enumerate(station_list[::-1]):
        station_dict[x] = i

    # remake the list 
    station_order = list(station_dict.keys())
    # reorder the list
    station_order = station_order[::-1]

    if 'S0' in station_order:
        # try at the beginning
        station_list = df['orig_station'].iloc[:5].to_list()
        station_dict = OrderedDict()
        # add them to the dictionary 
        for i, x in enumerate(station_list):
            station_dict[x] = i

        # remake the list 
        station_order = list(station_dict.keys())
        if 'S0' in station_order:
            logging.info('WARNING: Unable to determine station order {}'.format(station_order))
            station_order.remove('S0')

    logging.info('INFO: Station order={}'.format(
      station_order))


    config.station_order = station_order

# 
# # ['S16', 'S17', 'S12', 'S15', 'S16']
# 
# # ['S17', 'S12', 'S15', 'S16']
# # 17 12 15 16 17

# def all_station_order(df,gzip_filename):
#     if ('wtn.1.41.csv') in gzip_filename:
#         station_order = ['S12', 'S15', 'S16', 'S17']
#         logging.info('WARNING: Setting station order manually={}'.format(
#           station_order))
# 
#     elif ('wtn.1.44.csv') in gzip_filename:
#         station_order = ['S12', 'S15', 'S16', 'S17']
#         logging.info('WARNING: Setting station order manually={}'.format(
#           station_order))
# 
#     elif ('wtn.1.50.csv') in gzip_filename:
#         station_order = ['S12', 'S15', 'S16', 'S17']
#         logging.info('WARNING: Setting station order manually={}'.format(
#           station_order))
# 
#     elif ('wtn.10.14.csv') in gzip_filename:
#         station_order = ['S12', 'S15', 'S16', 'S14']
#         logging.info('WARNING: Setting station order manually={}'.format(
#           station_order))
# 
#     elif ('wtn.10.18.csv') in gzip_filename:
#         station_order = ['S12', 'S16', 'S15', 'S14']
#         logging.info('WARNING: Setting station order manually={}'.format(
#           station_order))
# 
# 
# 
#     elif ('wtn.10.20.csv') in gzip_filename:
#         station_order = ['S12', 'S16', 'S15', 'S14']
#         logging.info('WARNING: Setting station order manually={}'.format(
#           station_order))
# 
#     elif ('wtn.10.25.csv') in gzip_filename:
#         station_order = ['S12', 'S16', 'S15', 'S14']
#         logging.info('WARNING: Setting station order manually={}'.format(
#           station_order))
# 
#     else: 
# 
#         station_list = list(np.unique(df['orig_station'].to_numpy()))
#         if 'S0' in station_list:
#             station_list.remove('S0')
#         if len(station_list) == 5:
#             station_order = ['S12', 'S15', 'S16', 'S14', 'S17']
#         else:
#             station_order = []
#             for i in range(0,6):
#                 station = df['orig_station'].iloc[i]
#                 if station not in station_order:
#                     station_order.append(station)
# 
#             logging.info('WARNING: Station order={}'.format(
#               station_order))
# 
#     return station_order

def all_station_check(df,station_order):
    # depreciated, just too slow
    # check that stations are recorded to the right frame group and 
    # correct if necessary 
        
    len_stations = len(station_order)
    sta_dict = {}

    for i, sta in enumerate(station_order):
        frame_group = list(range(i+1,61,len_stations))
        sta_dict[sta] = frame_group

    for i, sta in enumerate(station_order):
        # filter on isin frame group 
        df_st = df[df['orig_frame'].isin(sta_dict[sta])]

        idx_list = df_st[df_st['orig_station'] != sta].index.tolist()
        idx_list.extend(df_st[pd.isna(df_st['orig_station'])].index.tolist())
        percent = 100*len(idx_list)/len(df)

        if len(idx_list) > 0 and percent < 0.01: 
            logging.info('WARNING: Correcting station for these records\n{}'.format(df.iloc[idx_list].to_string()))    
            for idx in idx_list:
                df.at[idx,'orig_station'] = sta
        elif percent >= 0.01:          
            logging.info('EXCEPTION: Too many stations recorded incorrectly?\n{}'.format(df.iloc[idx_list].head(100).to_string()))
            raise Exception

    return df

def all_delete_damaged(df):

    # rec_damaged_sync_total = ((df['sync'] != '1110001001000011101101') | (df['clock_flag'].isna()) | (df['orig_timestamp'].isna()) | (df['frame'] == -99)
    #   | ((df['bit_synchronizer'] == '00101') & (df['clock_flag'] == 1))).sum()

    no_to_shift = len(config.station_order)

    rec_damaged_sync = df[(df['sync'] != '1110001001000011101101') 
        | (df['clock_flag'].isna()) 
        | (df['orig_timestamp'].isna()) 
        | (df['orig_station'] == 'S0') 
        | (df['frame'] == -99)
        | (df['sync'].shift(-no_to_shift) != '1110001001000011101101')].index.to_list()

    rec_damaged_sync_total = len(rec_damaged_sync)

    if rec_damaged_sync_total > 0:
        logging.info('WARNING: Removing {} damaged record(s) - sync code is incorrect or key information is missing'.format(
          rec_damaged_sync_total))
        # df.drop(df[((df['sync'] != '1110001001000011101101') | (df['clock_flag'].isna()) | (df['orig_timestamp'].isna()) | (df['frame'] == -99)  
        #   | ((df['bit_synchronizer'] == '00101') & (df['clock_flag'] == 1)))].index , inplace=True)
        df.drop(df.index[rec_damaged_sync], inplace=True)
        df.reset_index(inplace=True,drop=True)

    return df, rec_damaged_sync_total

def all_drop_duplicates(df):

    df_duplicates = df[df.duplicated(subset =['orig_timestamp',
      'orig_station','corr_ground_station', 'frame'], keep='first')]
    rec_duplicates = len(df_duplicates)

    if rec_duplicates > 0: 
        logging.info('WARNING: Found {} duplicates'.format(rec_duplicates))
        perfect_dup = df[df.duplicated(subset =['orig_timestamp',
          'orig_station','corr_ground_station',
          'orig_mh1_1','orig_mh2_1','orig_mhz_1',
          'orig_mh1_2','orig_mh2_2','orig_mhz_2',
          'orig_mh1_3','orig_mh2_3','orig_mhz_3',
          'orig_mh1_4','orig_mh2_4','orig_mhz_4',
          'frame'
          ], keep='first')]
        perfect_dup_count = len(perfect_dup)
        logging.info('WARNING: Found {} perfect duplicates'.format(perfect_dup_count))
        if perfect_dup_count != rec_duplicates:
            logging.info('FAIL: More duplicates than perfect duplicates'.format(perfect_dup_count))

        # df_dup = df[df.duplicated(subset =['orig_timestamp',
        #   'orig_station','corr_ground_station','frame'
        #   ], keep='first')]

        # view the df_duplicates
        # logging.info('Duplicates \n{}'.format(perfect_dup.to_string()))

        # df_duplicates_keep = df[df.duplicated(subset =['orig_timestamp',
        #   'orig_station','corr_ground_station', 'frame'], keep=False)]

        # remove the duplicates
        # before = len(df)
        # # logging.info('WARNING: Found {} duplicates\n{}'.format(rec_duplicates, df_duplicates_keep.to_string()))
        # df = df.drop_duplicates(subset =['orig_timestamp',
        #   'orig_station','corr_ground_station'], keep='first')
        # after = len(df)
        # logging.info('WARNING: Found {} duplicates'.format(len(perfect_dup_keep_first)))

        df.reset_index(inplace=True,drop=True)

    return df, rec_duplicates

def all_station_duplicates_v2(df,gzip_filename):
    logging.debug('DEBUG: all_station_duplicates_v2')

    # this code relies on the stations being regularly spaced
    # so it will not work if called after deleting any records 

    # if we find a lot of duplicate times at different stations
    # it is very likely that there is a problem
    rec_station_duplicates = 0

    # in wtn.10.29.csv there is a problem where S14 copies S12 for a while
    # test_rec_station_duplicates = df.duplicated(subset =['orig_timestamp',
    #   'corr_ground_station','frame'], keep='first').sum()

    df_filter = df.copy()

    # create a hash on the combined columns
    df_filter['hash'] = pd.util.hash_pandas_object(df[['corr_ground_station',
        'orig_mh1_1','orig_mh2_1','orig_mhz_1',
        'orig_mh1_2','orig_mh2_2','orig_mhz_2',
        'orig_mh1_3','orig_mh2_3','orig_mhz_3',
        'orig_mh1_4','orig_mh2_4','orig_mhz_4',
        'frame']], 
        index=False)

    len_stations = len(config.station_order)    
    if len_stations == 0:
        logging.info('EXCEPTION: Station order not set')
        raise Exception

    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0

    df_filter['shift1'] = False
    df_filter['shift2'] = False
    df_filter['shift3'] = False
    df_filter['shift4'] = False

    # shift the hash to look for differences 
    if len_stations > 1:
        # df_filter['hash'] = df_filter['hash'].astype('int')
        df_filter['shift1'] = (df_filter['hash'].shift(1, fill_value=0) == df_filter.hash)
        count1 = (df_filter.shift1 == True).sum()
        if count1 > 0:
            df_filter['shift1a'] = (df_filter['shift1'].shift(-1, fill_value=0) == True)
            df_filter1 = df_filter[(df_filter.shift1 == True) | (df_filter.shift1a == True)]
            logging.info('SEVERE: Found {} repeated lines\n'.format(count1*2))
            logging.debug(df_filter1.to_string())

    if len_stations > 2:
        df_filter['shift2'] = (df_filter['hash'].shift(2, fill_value=0) == df_filter.hash)
        count2 = (df_filter.shift2 == True).sum()
        if count2 > 0:
            df_filter['shift2a'] = (df_filter['shift2'].shift(-2, fill_value=0) == True)
            df_filter2 = df_filter[(df_filter.shift2 == True) | (df_filter.shift2a == True)]
            logging.info('SEVERE: Found {} repeated lines\n'.format(count2*2))   
            logging.debug(df_filter2.to_string())    

    if len_stations > 3:
        df_filter['shift3'] = (df_filter['hash'].shift(3, fill_value=0) == df_filter.hash)
        count3 = (df_filter.shift3 == True).sum()
        if count3 > 0:
            df_filter['shift3a'] = (df_filter['shift3'].shift(-3, fill_value=0) == True)
            df_filter3 = df_filter[(df_filter.shift3 == True) | (df_filter.shift3a == True)]
            logging.info('SEVERE: Found {} repeated lines\n'.format(count3*2)) 
            logging.debug(df_filter3.to_string())       

    if len_stations > 4:
        df_filter['shift4'] = (df_filter['hash'].shift(4, fill_value=0) == df_filter.hash)
        count4 = (df_filter.shift4 == True).sum()
        if count4 > 0:
            df_filter['shift4a'] = (df_filter['shift4'].shift(-4, fill_value=0) == True)
            df_filter4 = df_filter[(df_filter.shift4 == True) | (df_filter.shift4a == True)]
            logging.info('SEVERE: Found {} repeated lines\n'.format(count4*2))   
            logging.debug(df_filter4.to_string())     

    rec_station_duplicates = count1 + count2 + count3 + count4
    if rec_station_duplicates > 0:
        # find the records to drop (use the ones where shift == True, because that means that we can keep the first ones)
        df_drop = df_filter[(df_filter.shift1 == True) | (df_filter.shift2 == True) | (df_filter.shift3 == True) | (df_filter.shift4 == True)]
        df.drop(df_drop.index,inplace=True)
        logging.info('WARNING: Dropping {} repeated lines'.format(len(df_drop)))
        logging.debug('-----------')
        logging.debug('WARNING: Dropping {} repeated lines \n{}'.format(len(df_drop),df_drop.to_string()))
        logging.debug('-----------')
        df.reset_index(inplace=True,drop=True)

        rec_station_duplicates = len(df_drop)

        # if 'test_import.py' not in gzip_filename:
        #     logging.info('EXCEPTION: Found repeated lines')      
        #     raise Exception

    return df, rec_station_duplicates    



  

def all_flat_seismogram(df,gzip_filename):
    logging.debug('DEBUG: all_flat_seismogram')

    # this code relies on the stations being regularly spaced
    # so it will not work if called after deleting any records 

    # if we find a lot of duplicate times at different stations
    # it is very likely that there is a problem
    # flat_seismogram_count = 0
   

    df_filter = df.copy()
    
    # create a hash on the combined columns, without the frame
    df_filter['hash'] = pd.util.hash_pandas_object(df[['corr_ground_station',
        'orig_mh1_1','orig_mh2_1','orig_mhz_1',
        'orig_mh1_2','orig_mh2_2','orig_mhz_2',
        'orig_mh1_3','orig_mh2_3','orig_mhz_3',
        'orig_mh1_4','orig_mh2_4','orig_mhz_4']], 
        index=False)

    len_stations = len(config.station_order)    
    if len_stations == 0:
        logging.info('EXCEPTION: Station order not set')
        raise Exception

    # number of (perfect) cumulative records before there is a warning
    max_test = 5
    # for station in ['S16']:
    
    local_station_order = config.station_order.copy()
    if 'S17' in local_station_order:
        local_station_order.remove('S17')
    for station in local_station_order:
        df_filter1 = df_filter[df_filter.orig_station == station].copy()
        # shift the hash 
        df_filter1['shift_neg'] = (df_filter1['hash'].shift(-1, fill_value=0) == df_filter1.hash)
        df_filter1['shift_pos'] = (df_filter1['hash'].shift(1, fill_value=0) == df_filter1.hash)
        
        # df_gst_e.reset_index(inplace=True,drop=True)

        df_filter1['shift_found'] = np.where((df_filter1.shift_pos == True) | (df_filter1.shift_neg == True), True, False)
        # get a cummulative sum of the hashes - suggests that it is really flat
        df_filter1['cumsum_shift'] = df_filter1['shift_found'] * (df_filter1['shift_found'].groupby((df_filter1['shift_found'] != df_filter1['shift_found'].shift()).cumsum()).cumcount() +1)        
        max1 = df_filter1['cumsum_shift'].max()
        if max1 > max_test:
            logging.info('WARNING: Many consecutive identical records found for station {} max={}'.format(station,max1))
        
        # s=df_filter1.groupby('hash').hash.count()
        # total_hash_count = len(s)
        # 100 records/100 hashes = good 
        # 1 hash = bad        
        
        # unique_rec_percentage = 100 * total_hash_count/len(df_filter1)
        # if unique_rec_percentage < 90:
        #    logging.info('SEVERE WARNING: Many identical records found for station {}'.format(station))

        # if (unique_rec_percentage < 90) or (max1 > max_test):
        #    logging.info('SEVERE WARNING: Many identical records found for station {} - full record\n{}'.format(station,df_filter1.to_string()))
        # else:
        #    logging.info('PASSED: Mainly unique records {}'.format(station))

        # something like this will get the number of repeated hashes     
        # s1 = s[s>2]
        # logging.info(s1.count())
        
    # if 'test_import.py' not in gzip_filename:
    #     logging.info('EXCEPTION: Found repeated lines')      
    #     raise Exception

def all_drop_station_duplicates(df,gzip_filename):
    rec_station_duplicates = 0
# depricated!!!
#     # remember to use the correct index when working on this method!
# 
#     logging.debug('DEBUG: all_drop_station_duplicates')
# 
#     # if we find a lot of duplicate times at different stations
#     # there may be a problem 
#     rec_station_duplicates = 0
# 
#     # in wtn.10.29.csv there is a problem where S14 copies S12 for a while
#     # test_rec_station_duplicates = df.duplicated(subset =['orig_timestamp',
#     #   'corr_ground_station','frame'], keep='first').sum()
# 
# 
# 
#     df_filter = df.copy()
# 
#     # filter to remove 'S17', because there are too many duplicates
#     df_filter = df_filter[(df_filter.orig_station != 'S17')]
# 
#     # logging.info(df_filter.head(100).to_string())
# 
#     t1 = datetime.now()
#     idx_duplicate = df_filter[df_filter.duplicated(subset =[
#       'corr_ground_station',
#       'orig_mh1_1','orig_mh2_1','orig_mhz_1',
#       'orig_mh1_2','orig_mh2_2','orig_mhz_2',
#       'orig_mh1_3','orig_mh2_3','orig_mhz_3',
#       'orig_mh1_4','orig_mh2_4','orig_mhz_4',
#       'frame'
#     ], keep='first')].index.tolist()
#     t2 = datetime.now()
#     logging.info('t2 : {}'.format((t2-t1).total_seconds()))
# 
#     if len(idx_duplicate) > 0:
#         idx_list_final_dup = []
# 
#         last_found_list = []
#         first_found_list = []
# 
#         df_filter_first = df.iloc[idx_duplicate].copy()
# 
#         # logging.info('first \n{}'.format(df_filter_first.to_string()))
# 
#         df_filter_first['corr_timestamp_low'] = df_filter_first['corr_timestamp'] - pd.Timedelta(seconds=0.5)
#         df_filter_first['corr_timestamp_high'] = df_filter_first['corr_timestamp'] + pd.Timedelta(seconds=0.5)
#         array = df_filter_first[['corr_timestamp_low', 'corr_timestamp_high', 'frame', 'orig_station','corr_timestamp', 'orig_mh1_1',  'orig_mh2_1',  'orig_mhz_1',  'orig_mh1_2',  'orig_mh2_2',  'orig_mhz_2',  'orig_mh1_3',  'orig_mh2_3',  'orig_mhz_3',  'orig_mh1_4',  'orig_mh2_4',  'orig_mhz_4']].to_numpy()
#         index_array = df_filter_first.index.to_numpy()
# 
#         # logging.info('len df_filter_first {}'.format(len(df_filter_first)))
# 
#         t1 = datetime.now()
#         for i in range(0,len(df_filter_first)):
#         # for i in range(68,69):
# 
#             corr_timestamp_low = array[i,0]
#             corr_timestamp_high = array[i,1]
#             frame = array[i,2]
#             orig_station = array[i,3]
#             corr_timestamp = array[i,4]
#             orig_mh1_1 = array[i,5]
#             orig_mh2_1 = array[i,6]  
#             orig_mhz_1 = array[i,7]  
#             orig_mh1_2 = array[i,8] 
#             orig_mh2_2 = array[i,9] 
#             orig_mhz_2 = array[i,10]
#             orig_mh1_3 = array[i,11] 
#             orig_mh2_3 = array[i,12] 
#             orig_mhz_3 = array[i,13] 
#             orig_mh1_4 = array[i,14] 
#             orig_mh2_4 = array[i,15]  
#             orig_mhz_4 = array[i,16]
#             idx = index_array[i]
# 
#             # logging.info('{} {} {}'.format(corr_timestamp,idx,orig_station))
# 
#             # t_a = datetime.now()
# 
# 
#             idx_list_dup = df[(df.corr_timestamp > corr_timestamp_low) & 
#                 (df.corr_timestamp < corr_timestamp_high) & 
#                 (df.frame == frame) & 
#                 (df.orig_station != orig_station) &
#                 (df.orig_mh1_1 == orig_mh1_1) &
#                 (df.orig_mh2_1 == orig_mh2_1) &
#                 (df.orig_mhz_1 == orig_mhz_1) &
#                 (df.orig_mh1_2 == orig_mh1_2) &
#                 (df.orig_mh2_2 == orig_mh2_2) &
#                 (df.orig_mhz_2 == orig_mhz_2) &
#                 (df.orig_mh1_3 == orig_mh1_3) &
#                 (df.orig_mh2_3 == orig_mh2_3) &
#                 (df.orig_mhz_3 == orig_mhz_3) &
#                 (df.orig_mh1_4 == orig_mh1_4) &
#                 (df.orig_mh2_4 == orig_mh2_4) &
#                 (df.orig_mhz_4 == orig_mhz_4)
#                 ].index.tolist()
# 
#             if i == 0:
#                 t2 = datetime.now()
#                 logging.info('time after first iteration : {}'.format((t2-t1).total_seconds()))
#                 logging.info('len idx_list_dup {}'.format(len(idx_list_dup)))
#             if i % 100 == 0:
#                 t2 = datetime.now()
#                 logging.info('time so far {}: {}'.format(i,(t2-t1).total_seconds()))
#                 logging.info('len idx_list_dup {}: {}'.format(i,len(idx_list_dup)))
# 
#             # DON'T FORGET: idx_list_dup applies to df 
# 
#             # t_b = datetime.now()
#             # logging.info('t each old: {}'.format((t_b-t_a).total_seconds()))
# 
#             # check for a duplicate, which means that it is probably recording to the wrong ground station
#             if len(idx_list_dup) > 0:
#                 # logging.info(df.iloc[idx_list_dup])
#                 last_found_list.append(idx)
#                 first_found_list.extend(idx_list_dup)
# 
#         t3 = datetime.now()
#         logging.info('t end: {}'.format((t3-t1).total_seconds()))
# 
#         last_found_list.sort()
#         if len(last_found_list) > 0:
# 
#             logging.debug('WARNING: {} records found that are copied from another station to be deleted\n{}'.format(len(last_found_list),df.iloc[last_found_list].to_string()))
#             logging.debug('WARNING: {} records found that are copied to another station (will not be deleted)\n{}'.format(len(first_found_list),df.iloc[first_found_list].to_string()))
# 
#             df.drop(df.index[last_found_list],inplace=True)
#             df.reset_index(inplace=True,drop=True)
#             rec_station_duplicates = len(last_found_list)
# 
# 
# 
    return df, rec_station_duplicates

########################
# example of vectorized search and update
# def soc_iter(TEAM,home,away,ftr):
#     df['Draws'] = 'No_Game'
#     df.loc[((home == TEAM) & (ftr == 'D')) | ((away == TEAM) & (ftr == 'D')), 'Draws'] = 'Draw'
#     df.loc[((home == TEAM) & (ftr != 'D')) | ((away == TEAM) & (ftr != 'D')), 'Draws'] = 'No_Draw'
# 
# 
# df['Draws'] = soc_iter('Arsenal', df['HomeTeam'].values, df['AwayTeam'].values,df['FTR'.values])

# split the dataframe into station and ground station
def all_split(df,gzip_filename):

    df_st_list = []

    # find where the ground station changes
    gs_range = (df['corr_ground_station'].ne(df['corr_ground_station'].astype('int').shift())).to_frame().apply(lambda x: x.index[x].tolist())
    gs_range_first = gs_range['corr_ground_station'].tolist()

    # when starting the tapes (orig_no == 1), the ground station has shifted
    # (this may be that the ground station is not recorded correctly)
    # or the trace contains duplicates because of errors 
    df['orig_no_shift'] = df['orig_no'].shift(1)
    gs_range_additional1 = df[(df['orig_no'] == 1) & (df['orig_no_shift'] > 2)].index.tolist()

    gs_range_additional1 = set(gs_range_additional1)
    gs_range_additional = list(gs_range_additional1.difference(gs_range_first))
    gs_range_additional.sort()

    if len(gs_range_additional) > 0:
        logging.info('Possible different ground stations?:')
        for a in gs_range_additional:  
            logging.info('\n{}'.format(df.iloc[a-1:a+2].to_string()))

    gs_range_first.extend(gs_range_additional)
    # gs_range_first = list(set(gs_range_first))
    gs_range_first.sort()
    df.drop(['orig_no_shift'],axis=1,inplace=True)     
    # 
    gs_range_last = []
    for i, first in enumerate(gs_range_first[1:]):
        gs_range_last.append(first-1)
    gs_range_last.append(len(df)-1)

    for a, b in zip(gs_range_first, gs_range_last):
        # check for possible errors in the ground stations 
        # logging.info(gs_range_additional)
        if a in gs_range_additional:
            # logging.info(a)
            # prev_ground_station = df['corr_ground_station'].iloc[a-1]
            current_ground_station = df['corr_ground_station'].iloc[a]
            # logging.info(current_ground_station)
            # logging.info(extra_ground_station)
            new_ground_station = get_new_ground_station(ground_station1=current_ground_station)
            logging.info('WARNING: Ground station - {} Changing ground station to {}.'.format(
                current_ground_station, new_ground_station))
            df['corr_ground_station'].iloc[a:b+1] = new_ground_station
        # logging.info('{} {}'.format(a, b))
        # logging.info('\n{}'.format(df.iloc[a:a+1].to_string()))
        # logging.info('\n{}'.format(df.iloc[b:b+1].to_string()))
        # range_diff = b - a
        # if (range_diff) < 10:
        #     logging.info('Too few records for the ground station: {}'.format(range_diff))

    logging.info('Ground stations - First and last records')
    
    for a, b in zip(gs_range_first, gs_range_last):
        # logging.info('{} {}'.format(a, b))
        logging.info('n={}\n{}'.format(b-a,df.iloc[a:a+1].to_string()))
        logging.info('\n{}'.format(df.iloc[b:b+1].to_string()))
        # range_diff = b - a
        # if (range_diff) < 10:
        #     logging.info('Too few records for the ground station: {}'.format(range_diff))

    for a, b in zip(gs_range_first, gs_range_last):
        # get each ground station 
        df_gst = df.iloc[a:b+1].copy(deep=True)

        ground_station = df_gst['corr_ground_station'].iloc[0]

        # step 6 - test for a small number of records - probably invalid
        gst_len = len(df_gst)
        if gst_len < 60: 
            # if the records are errors, use 'continue'
            # if they are fine, use 'pass'
            # this one has been checked 
            # if 'wtn.10.18.csv' in gzip_filename:
            #     logging.info('WARNING: Only {} record(s) found for ground station {}'.format(gst_len,ground_station))
            #     continue
            # if 'wtn.10.39.csv' in gzip_filename:
            #     logging.info('WARNING: Only {} record(s) found for ground station {}'.format(gst_len,ground_station))
            #     continue
            if 'wtn.11.7.csv' in gzip_filename:
                pass
            else:
                logging.info('WARNING: Only {} record(s) found for ground station {} (records will not be used)'.format(gst_len,ground_station))
                continue
                # raise Exception

        grouped_st = df_gst.groupby('orig_station')
        # group by station
        for st in grouped_st.groups.keys():

            if st == 'S17':
                continue

            df_st = grouped_st.get_group(st)

            # make a new index
            df_st.reset_index(inplace=True,drop=True)
            global total
            total += len(df_st)

            # get the record for the ground station (the whole ground station,
            # not just for this record)

            last_station = df_gst['orig_station'].iloc[-1]
            last_timestamp = df_gst['orig_timestamp'].iloc[-1]
            
            st_dict = {
            'df_gst' : df_st,
            'last_station' : last_station,
            'last_timestamp' : last_timestamp
            }
            # add to list of groups 
            df_st_list.append(st_dict)


    return df_st_list

def station_consec_frames_inner(df_gst,valid_idx,valid_frame):
    # make a column of consecutive frames 
    total_len = len(df_gst)
    frames = np.arange(90)
    frame_check = np.tile(frames,total_len+1)

    mod1 = valid_idx % 90
    start_frame = valid_frame - mod1
    if start_frame < 90:
        start_frame += 90

    fr = frame_check[start_frame:start_frame+total_len]
    df_gst['frame_tile'] = fr

    return df_gst 


def station_fix_bad_timestamps(df_bad,df_dropped,df_orig):

    compliance_passed = True

    # YYYY

    # logging.info('Before\n{}'.format(df_bad.to_string()))

    # logging.info(df_bad.to_string())

    # logging.info('df_bad\n{}'.format(df_bad.to_string()))
    # initial_len = len(df_bad)
    # if initial_len < 3:
    #     logging.info('EXCEPTION: station_fix_bad requires at least three records')
    #     raise Exception

    rec_adjusted_timestamps = 0

    df_bad.reset_index(inplace=True,drop=True)

    logging.debug('df_bad')
    logging.debug(df_bad.to_string())

    corr_ground_station  = df_bad.corr_ground_station.iloc[0]
    orig_station = df_bad.orig_station.iloc[0]

    # time_diff has to be calculated from the last good record to the next good record 
    time_diff = (df_bad['corr_timestamp'].iloc[-1] -  df_bad['corr_timestamp'].iloc[0]).total_seconds()
    # count the frame increments (from second to end)
    sum_diff = df_bad['frame_change'].iloc[1:].sum()

    # logging.info('very bad')
    # logging.info(df_bad.to_string())

    # if sum_diff == 0:
    #     corr_ground_station  = df_test.corr_ground_station.iloc[0]
    #     orig_station = df_test.orig_station.iloc[0]
    #     logging.info('WARNING: Ground station/Station - {} {} Unable to fix the timing because there is no frame change'.format(corr_ground_station, orig_station,df_test.to_string()))
    #     # logging.info('unable')
    #     # logging.info(df_gst.to_string())
    #     # logging.info(df_bad.to_string())
    #     # global df_orig_shown
    #     # if df_orig_shown == False:
    #     #     logging.debug('--------------------------')
    #     #     logging.debug(df_orig.to_string())
    #     #     df_orig_shown = True
    #     #     logging.debug('--------------------------')
    #     compliance_passed = False
    #     # return from the error
    #     return df_bad, df_dropped, rec_adjusted_timestamps, compliance_passed`

    if sum_diff == 0:
        # logging.debug('--------------------------')
        # logging.debug(df_orig.to_string())
        # df_orig_shown = True
        # logging.debug('--------------------------')
        # raise Exception 
        logging.info('WARNING: Ground station/Station - {} {} Repeated frame'.format(corr_ground_station, orig_station))
        # logging.info('unable')
        # logging.info(df_gst.to_string())
        # logging.info(df_bad.to_string())
        # global df_orig_shown
        # if df_orig_shown == False:
        #     logging.debug('--------------------------')
        #     logging.debug(df_orig.to_string())
        #     df_orig_shown = True
        #     logging.debug('--------------------------')
        compliance_passed = False
        # return from the error
        return df_bad, df_dropped, rec_adjusted_timestamps, compliance_passed
    

    single_gap = time_diff/(sum_diff)



    # start from the last timestamp that was correct
    orig_timestamp=df_bad['corr_timestamp'].iloc[0]
    frame_change = 0

    # estimate the end time
    end_time = df_bad['corr_timestamp'].iloc[0] + pd.Timedelta(seconds=sum_diff*DELTA*4)
    logging.info('INFO: End time estimate = {}'.format(end_time))

    idx = df_bad.first_valid_index() + 1
    first_i = idx+1
    # print('idx ', idx)
    # logging.info('iloc')
    # logging.info(df_bad.iloc[0:1].to_string())
    # logging.info('loc')
    # logging.info(df_bad.loc[idx:idx+1].to_string())
    # df_bad.at[24,'orig_mh1_1'] = 42
    # df_bad.at[25,'end'] = 'cheese puff'

    # this test is a bit tighter than calculate gaps 
    # .at uses the label position, so we need the real index numbers




    # some tests at the beginning 
    last_idx = len(df_bad) - 1
    indices = [0,last_idx]
    df_test = df_bad.iloc[indices].copy()
    df_test.reset_index(inplace=True,drop=True)
    # print('AA')
    df_test, gaps_long, gaps_8888 = calculate_gaps(df_test)
    
    # logging.info('df_test')
    # logging.info(df_test.to_string())
    # 
    # logging.info(df_test.to_string())


    if df_test.corr_gap_count.iloc[-1] == -8888:

        logging.info('WARNING: Ground station/Station - {} {} Unable to fix the timing'.format(corr_ground_station, orig_station,df_test.to_string()))
        # logging.info('unable')
        # logging.info(df_gst.to_string())
        # logging.info(df_bad.to_string())
        # global df_orig_shown
        # if df_orig_shown == False:
        #     logging.debug('--------------------------')
        #     logging.debug(df_orig.to_string())
        #     df_orig_shown = True
        #     logging.debug('--------------------------')
        compliance_passed = False
        # return from the error
        return df_bad, df_dropped, rec_adjusted_timestamps, compliance_passed

    else:
        pass2 = False

        # logging.info('This is df_test:')
        # logging.info(df_test.to_string())

        if single_gap > 0.6019  and single_gap < 0.6061 and (df_bad['frame_change'] != 0).all():
            for current_i in range(idx, idx+len(df_bad)-2):
                # print(current_i, single_gap)
                # update the timestamps 
                frame_change += df_bad['frame_change'].loc[current_i]
                new_timestamp = orig_timestamp + pd.Timedelta(seconds=(frame_change*single_gap))
                df_bad.at[current_i,'corr_timestamp'] = new_timestamp
                # df_bad.at[current_i,'orig_mh1_1'] = 42
                rec_adjusted_timestamps += 1
    
            # fix any gaps
            df_bad = station_fix_missing_timestamps2(df_bad)
    
            # now make a test to see if it worked
            # print('BB')
            df_bad, gaps_long, gaps_8888 = calculate_gaps(df_bad)
            bad_count = (df_bad['corr_gap_count'] != 1).sum()
            if bad_count > 0:
                pass2 = False
            else:
                pass2 = True
    
        if pass2 == False:
            to_drop = list(range(1,len(df_bad)-1))
            # logging.info('Dropping these {}'.format(len(to_drop)))

            df_dropped1 = df_bad.iloc[to_drop]
            df_dropped = pd.concat([df_dropped, df_dropped1], ignore_index = False)


            # df_bad.drop(df_bad.index[to_drop],inplace=True)
            # df_bad.reset_index(inplace=True,drop=True)
            # df_bad, gaps_long, gaps_8888 = calculate_gaps(df_bad)

            # logging.info('Calculate gaps seems wrong')
            # logging.info(df_bad.to_string())
    
            # logging.info('Pass2 is False 1')
            # logging.info(df_test.to_string())

            # fix any gaps
            df_bad = station_fix_missing_timestamps2(df_test)
    
            # logging.info('Pass2 is False 2')
            # logging.info(df_bad.to_string())

    # logging.info('End\n{}'.format(df_bad.to_string()))

    return df_bad, df_dropped, rec_adjusted_timestamps, compliance_passed

def calculate_delta4(df_gst):

    t1 = datetime.now()

    delta4_arr = np.empty(len(df_gst))
    delta4_arr[:] = np.NaN

    # Estimate the delta4 value 
    # This makes estimates where there are good data and interpolates 
    # into any gaps.
    # df_gst['delta4'] = np.NaN

    idx_list = df_gst[df_gst['frame'] == 0].index.tolist()
    for i, idx_0 in enumerate(idx_list[1:]): 
        idx_0_prev = idx_list[i]
        cumsum_0 = df_gst['cumsum'].iloc[idx_0]
        cumsum_0_prev = df_gst['cumsum'].iloc[idx_0_prev]
        if cumsum_0_prev + 90 == cumsum_0:
            # this is suitable to estimate delta4
            timestamp_0 = df_gst['orig_timestamp'].iloc[idx_0]    
            timestamp_0_prev = df_gst['orig_timestamp'].iloc[idx_0_prev]   
            time_diff = (timestamp_0-timestamp_0_prev)
            delta4 = (time_diff/90).total_seconds()
            # check delta4 is reasonable 
            if (delta4 > (DELTA*4)*1.01) or (delta4 < (DELTA*4)*0.99):
                logging.info('WARNING: delta4 not suitable id=idx_0 delta4={}'.format(idx_0, delta4))
                print('WARNING: delta4 not suitable id=idx_0 delta4={}'.format(idx_0, delta4))
            else: 
                delta4_arr[idx_0_prev:idx_0] = delta4


    # t2 = datetime.now()
    # logging.info('timeit delta1 v2: {}'.format((t2-t1).total_seconds()))

    # t1 = datetime.now()    

    # make sure we've found some delta4 values  
    delta4_count = (~(np.isnan(delta4_arr))).sum()
    
    if delta4_count  < 1:
        # if not, try again, ignoring clock flag (cumsum is only acculated when 
        # cumsum is not set)
        for i, idx_0 in enumerate(idx_list[1:]): 
            idx_0_prev = idx_list[i]
            corr_gap_count_test = (df_gst['corr_gap_count'].iloc[idx_0_prev:idx_0] == 1).all()
            if idx_0_prev + 90 == idx_0 and corr_gap_count_test:
                # this is suitable to estimate delta4
                timestamp_0 = df_gst['orig_timestamp'].iloc[idx_0]    
                timestamp_0_prev = df_gst['orig_timestamp'].iloc[idx_0_prev]   
                time_diff = (timestamp_0-timestamp_0_prev)
                delta4 = (time_diff/90).total_seconds()
                # print(idx_0_prev,idx_0,time_diff.total_seconds(),delta4)
                delta4_arr[idx_0_prev:idx_0] = delta4
                # check delta4 is reasonable 
                if (delta4 > (DELTA*4)*1.01) or (delta4 < (DELTA*4)*0.99):
                    logging.info('WARNING: delta4 not suitable id=idx_0 delta4={}'.format(idx_0, delta4))
                    print('WARNING: delta4 not suitable id=idx_0 delta4={}'.format(idx_0, delta4))
                else: 
                    delta4_arr[idx_0_prev:idx_0] = delta4

    # t2 = datetime.now()
    # logging.info('timeit delta3 : {}'.format((t2-t1).total_seconds()))


    delta4_count = (~(np.isnan(delta4_arr))).sum()
    df_gst['delta4'] = delta4_arr
    # if something is still not found, use the nominal version. 
    if delta4_count  > 0:
        df_gst['delta4'].interpolate(method='linear', inplace=True)
        df_gst['delta4'].backfill(inplace=True)
    else: 
        orig_station = df_gst['orig_station'].iloc[0]
        corr_ground_station = df_gst['corr_ground_station'].iloc[0]
        logging.info('WARNING: Ground station/Station - {} {} Using the nominal value of DELTA.'.format(corr_ground_station, orig_station))
        df_gst['delta4'] = DELTA*4

    return df_gst


def calculate_cumsum2(df_gst):

    orig_station = df_gst['orig_station'].iloc[0]
    corr_ground_station = df_gst['corr_ground_station'].iloc[0]

    # print('CC')
    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

    # find the 'good' records (this now includes records that have been fixed)
    # and doesn't take into account the clock flag
    df_gst['good2'] = np.where((df_gst['corr_gap_count'] == 1 ), True, False)
    # ignore records that are pretty close to the tolerance 
    # this will only allow the code to connect up the beginning and end of 
    # a cumulative section which is joined by something which is out of tolerance
    df_gst['good2'] = np.where((df_gst['corr_gap'] > config.lower_tolerance ) & (df_gst['corr_gap'] < config.higher_tolerance ) & (df_gst['frame_change'] == 1 ),True, df_gst['good2'] )



 # if corr_gap > 0.5538 and corr_gap < 0.6538 and frame_change == 1:

    df_gst['cumsum2'] = df_gst['good2'] * (df_gst['good2'].groupby((df_gst['good2'] != df_gst['good2'].shift()).cumsum()).cumcount() + 1)

    # check for any -7777 errors
    # print('DD')
    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
    count_7777 = (df_gst['corr_gap_count'] == -7777).sum()
    if count_7777 > 0:
        # find the 'good' records, this time including -7777 errors
        df_gst['good2'] = np.where((df_gst['corr_gap_count'] == 1 ) | (df_gst['corr_gap_count'] == -7777 ), True, False)
        # df_gst['good'] = np.where((df_gst['clock_flag'] == 0 ), df_gst['good'] , False)
        df_gst['cumsum2'] = df_gst['good2'] * (df_gst['good2'].groupby((df_gst['good2'] != df_gst['good2'].shift()).cumsum()).cumcount() + 1)    

    # logging.debug('Cumsum2--------------------------')
    # logging.debug(df_gst.to_string())
    # # df_orig_shown = True
    # logging.debug('--------------------------')

    # XXXX
    
    return df_gst

# def calculate_cumsum(df_gst,df_dropped):
def calculate_cumsum(df_gst):

    orig_station = df_gst['orig_station'].iloc[0]
    corr_ground_station = df_gst['corr_ground_station'].iloc[0]

    # print('CC')
    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

    # find the 'good' records
    df_gst['good'] = np.where((df_gst['corr_gap_count'] == 1 ), True, False)
    df_gst['good'] = np.where((df_gst['clock_flag'] == 0 ), df_gst['good'] , False)
    df_gst['cumsum'] = df_gst['good'] * (df_gst['good'].groupby((df_gst['good'] != df_gst['good'].shift()).cumsum()).cumcount() + 1)

    # check for any -7777 errors
    # print('DD')
    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
    count_7777 = (df_gst['corr_gap_count'] == -7777).sum()
    if count_7777 > 0:
        # find the 'good' records, this time including -7777 errors
        df_gst['good'] = np.where((df_gst['corr_gap_count'] == 1 ) | (df_gst['corr_gap_count'] == -7777 ), True, False)
        # df_gst['good'] = np.where((df_gst['clock_flag'] == 0 ), df_gst['good'] , False)
        df_gst['cumsum'] = df_gst['good'] * (df_gst['good'].groupby((df_gst['good'] != df_gst['good'].shift()).cumsum()).cumcount() + 1)    

# YYYY
        # logging.debug('--------------------------')
        # logging.debug('XcumulativeX')
        # logging.debug(df_gst.to_string())
        # # df_orig_shown = True
        # logging.debug('--------------------------')

    # valid_idx = None
    # 
    # if valid_idx is None: 
    #     test_i = int(len(df_gst) / 2) 
    #     df_filter = df_gst.iloc[test_i:]
    #     idx_list = df_filter[df_filter['cumsum'].gt(100)].index.tolist()
    #     if len(idx_list) > 0:
    #         # get the first good record after the middle value (cumsum greater than 100)
    #         valid_idx = idx_list[0]
    # 
    # # try successively easier tests 
    # if valid_idx is None: 
    #     idx_list = df_filter[df_filter['cumsum'].gt(10)].index.tolist()
    #     if len(idx_list) > 0:
    #         valid_idx = idx_list[0]
    # 
    # if valid_idx is None: 
    #     # keep trying easier tests - now ignore clock count
    #     df_gst['good'] = np.where((df_gst['corr_gap_count'] == 1 ), True, False)
    #     df_gst['cumsum'] = df_gst['good'] * (df_gst['good'].groupby((df_gst['good'] != df_gst['good'].shift()).cumsum()).cumcount() + 1)
    #     test_i = int(len(df_gst) / 2) 
    #     df_filter = df_gst.iloc[test_i:]
    #     idx_list = df_filter[df_filter['cumsum'].gt(100)].index.tolist()
    #     if len(idx_list) > 0:
    #         # get the first good record after the middle value (cumsum greater than 100)
    #         valid_idx = idx_list[0]
    # 
    # # keep trying successively easier tests - still ignoring clock count 
    # if valid_idx is None: 
    #     max1 = df_gst['cumsum'].max()
    #     if max1 > 3: 
    #         idx_list = df_gst[df_gst['cumsum'] == max1].index.tolist()
    #         if len(idx_list) > 0:
    #             valid_idx = idx_list[0]
    # 
    # if valid_idx is None: 
    # 
    #         # drop the records if there's nothing valid (or it's really 
    #         # short)
    #         df_dropped1 = df_gst
    #         logging.info("WARNING: Ground station/Station - {} {} Dropping {} record(s) because no valid index found ".format(
    #             corr_ground_station, orig_station, len(df_gst)))
    # 
    #         df_dropped = pd.concat([df_dropped, df_dropped1], ignore_index = False)
    #         df_gst = df_gst[0:0]
    #         valid_idx = None
    #         valid_frame = None
    #         valid_time = None
    # 
    # else: 
    # 
    #     cumsum = df_gst['cumsum'].iloc[valid_idx]
    # 
    #     # what was the maximum cumulative sum? 
    #     # a large number indicates that the transmission was working well for 
    #     # part of the time
    #     max1 = df_gst['cumsum'].max()
    # 
    # 
    # # if there's no suitable valid id (due to clock flags), put the middle one in - change valid id to middle anyway? - 30 mins 
    # #  	take out the 2000 thingy, make sure it's still working - 10 mins
    # 
    #     # print('reset valid index for testing')
    #     # valid_idx = 2000
    # 
    #     df_gst['valid_idx'] = pd.NA
    #     df_gst.at[valid_idx, 'valid_idx'] = 'VALID_IDX' 
    #     # print('setting ', valid_idx)
    #     # logging.info(df_gst.iloc[350:370].to_string())
    # 
    #     valid_frame = df_gst['frame'].iloc[valid_idx]
    #     valid_time = df_gst['corr_timestamp'].iloc[valid_idx]
    # 
    #     logging.info("INFORMATION: Ground station/Station - {} {} Using i={}. Consecutive 'good' records valid frame={} valid time={}, length={} cumsum={}".format(
    #         corr_ground_station, orig_station, valid_idx, valid_frame, valid_time, len(df_gst), cumsum))

    return df_gst
    # return df_gst, df_dropped, valid_idx, valid_frame, valid_time 

# def calculate_increasing(df_gst):
# 
#     logging.debug('DEBUG: calculate_increasing')
# 
#     # check for frame gaps which aren't increasing with a frame shift 
#     df_gst['actual_frame_gap_shift'] = df_gst['actual_frame_gap'] - df_gst['actual_frame_gap'].shift(1)  
#     df_gst.at[0,'actual_frame_gap_shift'] = 1
#     # negative frame gaps have actual_frame_gap_shift which is zero or negative
#     df_gst['sample_fit'] = np.where((df_gst['actual_frame_gap_shift'] <1), -4444, df_gst['sample_fit'])
# 
#     # also check for records which are consecutive with the first error 
#     idx_list_begin = (df_gst[df_gst['sample_fit'] == -4444]).index.tolist()
#     if len(idx_list_begin) > 0: 
#         for idx_begin in idx_list_begin:
#             begin_max = df_gst['actual_frame_gap'].iloc[idx_begin-1]
#             df_filter_begin = df_gst[idx_begin:]
#             # search for next record greater than the last good one
#             idx_list = (df_filter_begin[df_filter_begin['actual_frame_gap'] <= begin_max]).index.tolist()
#             # if any are found, the last valid record is the one before 
#             # if len(idx_list_zero) > 0: 
#             #     idx_end = idx_list_zero[0] -1
#             # else: 
#             #     # if not, just go to the end of the record 
#             #     idx_end = len(df_gst) - 1
#             # mark as -4444 error 
#             for i in idx_list:
#                 df_gst.at[i,'sample_fit'] = -4444
# 
#     return df_gst

# def check_non_increasing(df_gst,df_dropped):
# 
#     rec_deleted_non_incr = 0
# 
#     # check for non increasing entries 
#     df_gst = calculate_increasing(df_gst)
# 
#     idx_list = (df_gst[df_gst['sample_fit'] == -4444]).index.tolist()
#     if len(idx_list) > 0:
# 
#         orig_station = df_gst['orig_station'].iloc[0]
#         corr_ground_station = df_gst['corr_ground_station'].iloc[0]
# 
#         # percent to drop 
#         percent_error = 100*(len(idx_list)/len(df_gst))
#         # raise exception if dropping more than 0.5 percent
#         if percent_error > 0.5 and len(idx_list) > 20:
#             logging.info('SEVERE: Ground station/Station - {} {} Dropped more than 0.5% of records ({} of {}) because the timing is not increasing'.format(corr_ground_station, orig_station, len(idx_list),len(df_gst)))
#             print('SEVERE: Ground station/Station - {} {} Dropped more than 0.5% of records ({} of {}) because the timing is not increasing'.format(corr_ground_station, orig_station,len(idx_list),len(df_gst)))
#             logging.info(df_gst.to_string()) 
#         else:
#             logging.info('WARNING: Ground station/Station - {} {} Dropped {} record(s) because the timing is not increasing'.format(corr_ground_station, orig_station, len(idx_list)))
#         rec_deleted_timestamps = len(idx_list)
# 

# 
#         df_dropped1 = df_gst.iloc[idx_list]
#         df_dropped = pd.concat([df_dropped, df_dropped1], ignore_index = False)
# 
#         # drop non fitting records
#         df_gst.drop(df_gst.index[idx_list],inplace=True)
#         df_gst.reset_index(inplace=True,drop=True) 
#         df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
# 
#     return df_gst, df_dropped, rec_deleted_non_incr

# def check_non_unique_entries(df_gst,df_dropped):
# 
#     logging.debug('DEBUG: check_non_unique_entries')
# 
#     rec_deleted_non_uniq = 0
# 
# 
# 
#     # ignore the -4444 errors, because they will be removed anyway
#     # where ??? 
#     # df_filter = df_gst[~((df_gst['sample_fit'] == -2222) | (df_gst['sample_fit'] == -3333) (df_gst['sample_fit'] == -4444))]
#     df_gst['sample_fit'].fillna(0, inplace=True)
#     df_filter = df_gst[df_gst['sample_fit'] == 0]
#     # print(len(df_filter))
# 
#     # check for non unique entries 
#     df_dup = df_filter[df_filter.duplicated(subset =['actual_frame_gap'], keep='last')]
#     # df_dup = df_gst[df_gst.duplicated(subset =['actual_frame_gap'], keep=False)]
# 
#     if len(df_dup) > 0:
#         # print('temporarily raising an exception')
#         # logging.info('EXCEPTION: {} duplicate(s) found.'.format(len(df_dup)))
#         # logging.info(df_dup.to_string())
#         # print('EXCEPTION: {} duplicate(s) found.'.format(len(df_dup)))
#         # raise Exception
# 
#         rec_deleted_non_uniq = len(df_dup)
#         # drop the non unique entries 
# 
#         # view the duplicates - keep all for the view
#         df_dup_keep = df_filter[df_filter.duplicated(subset =['actual_frame_gap'], keep=False)]
#         # df_dup = df_gst[df_gst.duplicated(subset =['actual_frame_gap'], keep='last')]
#         logging.info('WARNING: {} duplicate(s) in actual_frame_gap\n{}'.format(rec_deleted_non_uniq, df_dup_keep.to_string()))
#         print('WARNING: {} duplicate(s) in actual_frame_gap'.format(rec_deleted_non_uniq))
# 
#         # TODO - check that the index works OK 
# 
#         # remove the duplicates - keep the last one 
#         df_dropped1 = df_gst.iloc[df_dup.index]
#         df_dropped = pd.concat([df_dropped, df_dropped1], ignore_index = False)
#         df_gst.drop_duplicates(subset =['actual_frame_gap'], keep='last', inplace=True)
# 
#         df_gst.reset_index(inplace=True,drop=True)
# 
#         # logging.info('WARNING: Removed {} duplicate(s)     in actual_frame_gap'.format(rec_deleted_non_uniq))
#         df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
# 
#     return df_gst, df_dropped, rec_deleted_non_uniq



# fix timestamps if possible 
# first remove any that don't have a gap at all 
# try to do something with the -8888 errors 
    # make begin and end sections

def make_good_dataframes(df_gst,low,high):
    lists = []
    for l, h in zip(low,high):
        l1 = df_gst.iloc[l:h+1].copy()
        l1.reset_index(inplace=True,drop=True)
        lists.append(l1)

    return lists 

def make_bad_dataframes(df_gst,low,high):
    lists = []
    low_bad = low.copy()
    high_bad = high.copy()

    high_bad.insert(0,-1)
    low_bad.append(len(df_gst))

    # add any records between good records to the bad list
    for h, l in zip(high_bad,low_bad):
        # print(h, l)
        l1 = df_gst.iloc[h+1:l].copy()
        # print(l1)
        if len(l1) > 0:
            lists.append(l1)

    return lists 

def find_good_records(lst):

    # search backwards for first value greater than cumsum_final_test
    end_idx = len(lst) - 1
    high = []
    low = []

    while True:

        idx = None
        try:
            idx_backward = next(x[0] for x in enumerate(lst[end_idx::-1]) if x[1] >= config.cumsum_final_test)
            idx = end_idx - idx_backward
            high.append(idx)
        except StopIteration:
            # didn't find any high values 
            pass

        low1 = None
        # if a value larger than the cumsum test was found, now find the first zero 
        if idx is not None:
            end_idx = idx
            low1 = 0
            try:
                idx_backward = next(x[0] for x in enumerate(lst[end_idx::-1]) if x[1] == 0)
                low1 = end_idx - idx_backward
            except StopIteration:
                # didn't find a zero (this can be OK when it is the first record)
                pass
            low.append(low1)

        if low1 is not None:
            end_idx = low1
        else:
            break

        if end_idx == 0:
            break

    low.reverse()
    high.reverse()

    return low, high

def split_good_records(df_gst):


    # logging.debug('split_good_records')
    # logging.debug(df_gst.to_string())
    # logging.debug('still a game??')

    cumsum2 = df_gst['cumsum2'].tolist()

    # logging.debug('--------------------------')
    # logging.debug(df_gst.to_string())
    # # df_orig_shown = True
    # logging.debug('--------------------------')

    low, high = find_good_records(cumsum2)
    good_lists = make_good_dataframes(df_gst,low,high)
    bad_lists = make_bad_dataframes(df_gst,low,high)

    # # remake the begin and end reocrds
    # df_gst['begin_end'] = 'None'
    # idx_list = (df_gst[(df_gst['cumsum2'] == 0)]).index.tolist()



    # if len(idx_list) > 0:
    # 
    #     begin_list = []
    #     end_list = []
    # 
    #     # make the begin and end markers 
    #     for idx in idx_list:
    #         prev_idx = idx - 1
    #         idx_plus = idx + config.cumsum_test + 1
    #         if prev_idx > 0:
    #             if df_gst['cumsum2'].iloc[prev_idx] > config.cumsum_test:
    #                 df_gst.at[prev_idx,'begin_end'] = 'end'
    #                 end_idx_orig = df_gst['orig_idx'].iloc[prev_idx]
    #                 end_list.append(end_idx_orig)
    #                 # print(prev_idx)
    #         if idx_plus < len(df_gst) - 1:
    #             if df_gst['cumsum2'].iloc[idx_plus] == config.cumsum_test + 1:
    #                 df_gst.at[idx,'begin_end'] = 'begin'
    #                 begin_idx_orig = df_gst['orig_idx'].iloc[idx]
    #                 begin_list.append(begin_idx_orig)
    #                 # print(idx)
    #         else:
    #             begin_idx_orig = df_gst['orig_idx'].iloc[len(df_gst) - 1]
    #             begin_list.append(begin_idx_orig)


    # logging.debug('--------------------------')
    # logging.debug(df_gst.to_string())
    # # df_orig_shown = True
    # logging.debug('--------------------------')
        
    # exit()
    # print('change this ')
    # df_gst['days'] = df_gst.corr_timestamp.dt.normalize()
    # gb = df_gst.groupby('days')    
    # df_list = [gb.get_group(x) for x in gb.groups]

    return good_lists, bad_lists

def fix_timestamps(df_gst,df_dropped):

    rec_adjusted_timestamps = 0
    compliance_passed = True

    orig_station = df_gst['orig_station'].iloc[0]
    corr_ground_station = df_gst['corr_ground_station'].iloc[0]


    # valid_idx_list = df_gst[df_gst['valid_idx'] == 'VALID_IDX'].index.tolist()
    # valid_idx = valid_idx_list[0]
    # valid_frame = df_gst['frame'].iloc[valid_idx]
    # valid_time = df_gst['corr_timestamp'].iloc[valid_idx]

    # logging.info('status now?')
    # # 
    # logging.debug('DEBUG all')
    # logging.debug(df_gst.to_string())
    # exit()

    # first remove any with no time gap at all
    idx_list_8888 = (df_gst[(df_gst['corr_gap_count'] == -8888)]).index.tolist()
    if len(idx_list_8888) > 0:
        to_drop = []
        # logging.debug('--------------------------')
        # logging.debug(df_gst.to_string())
        # global df_orig_shown
        # df_orig_shown = True
        # logging.debug('--------------------------')
        for idx in idx_list_8888:
            prev_corr_timestamp = df_gst.corr_timestamp.iloc[idx-1]
            curr_corr_timestamp = df_gst.corr_timestamp.iloc[idx]
            if prev_corr_timestamp == curr_corr_timestamp:
                to_drop.append(idx)

        if len(to_drop) > 0:
            df_dropped1 = df_gst.iloc[to_drop]
            df_dropped = pd.concat([df_dropped, df_dropped1], ignore_index = False)
            logging.info('WARNING: Ground station/Station - {} {} Dropping {} records with repeated timestamps'.format(corr_ground_station, orig_station, len(df_dropped1)))
            # logging.debug('--------------------------')
            # logging.debug(df_orig.to_string())
            # df_orig_shown = True
            # logging.debug('--------------------------')
            df_gst.drop(df_gst.index[to_drop],inplace=True)
            df_gst.reset_index(inplace=True,drop=True) 

            

    # try to do something with the -8888 errors 
    idx_list_8888 = (df_gst[(df_gst['corr_gap_count'] == -8888)]).orig_idx.tolist()
    if len(idx_list_8888) > 0:


        idx_list = (df_gst[(df_gst['cumsum'] == 0)]).index.tolist()
        if len(idx_list) > 0:

            begin_list = []
            end_list = []
                
    # YYYY
            # first make the begin and end markers 
            for idx in idx_list:
                prev_idx = idx - 1
                idx_plus = idx + config.cumsum_test + 1
                if prev_idx > 0:
                    if df_gst['cumsum'].iloc[prev_idx] > config.cumsum_test:
                        df_gst.at[prev_idx,'begin_end'] = 'end'
                        end_idx_orig = df_gst['orig_idx'].iloc[prev_idx]
                        end_list.append(end_idx_orig)
                        # print(prev_idx)
                if idx_plus < len(df_gst) - 1:
                    if df_gst['cumsum'].iloc[idx_plus] == config.cumsum_test + 1:
                        df_gst.at[idx,'begin_end'] = 'begin'
                        begin_idx_orig = df_gst['orig_idx'].iloc[idx]
                        begin_list.append(begin_idx_orig)
                        # print(idx)
                else:
                    begin_idx_orig = df_gst['orig_idx'].iloc[len(df_gst) - 1]
                    begin_list.append(begin_idx_orig)


        # df_orig is used to view the original dataframe if errors are found
        df_orig = df_gst.copy()


        # global df_orig_shown
        # print(df_orig_shown)
        # if df_orig_shown == False:
        #     logging.debug('--------------------------')
        #     logging.debug(df_orig.to_string())
        #     df_orig_shown = True
        #     logging.debug('--------------------------')

        begin_orig_idx = None
        # logging.info(df_gst.to_string())

        # find each -8888 error and try to fix
        # could rewrite this to be a little quicker using the begin and end lists, but this
        # is working, so I'm leaving it alone. 
        for idx_8888 in idx_list_8888:
    

            if begin_orig_idx is not None and idx_8888 <= begin_orig_idx:
                # already dealt with, so continue 
                # print('already dealt with, so continue')
                continue

            # can ignore -8888 if within wider tolerance 
            # if corr_gapcorr_gap

            # logging.info(idx_8888)
            # if idx_8888 == 113173 or idx_8888 == 113173.0:
            #     logging.info('Something found')

            
            # we use the original index, because records get inserted and 
            # deleted 
            orig_idx_list = df_gst[df_gst['orig_idx'] == idx_8888].index
            if len(orig_idx_list) > 0:
                i = orig_idx_list[0]
                # logging.info('a')
            else:
                # if not found, it is already deleted 
                # logging.info('b')
                continue

            # logging.info('pallaver')
            # logging.info(i)
            # logging.info(begin_orig_idx)

            # logging.info('problematic record = {} last begin = {}'.format(idx_8888,begin_orig_idx))
        
            # orig_idx_list = df_gst[df_gst['corr_gap'] == idx_8888].index
            
            # corr_gap = df_gst['corr_gap'].iloc[i]
            # frame_change = df_gst['frame_change'].iloc[i]
            # if the gap is pretty close to what it should be, just ignore it
            # if corr_gap > 0.5538 and corr_gap < 0.6538 and frame_change == 1:
            #     logging.info('WARNING - ignoring error close to tolerance \n{}'.format(df_gst.iloc[i-1:i+1].to_string()))
            #     continue 

            # if i == 113173 or i == 113173.0:
            #     logging.info('The one with the error')
            #     logging.info(df_gst.iloc[i:i+1].to_string())
            # find previous end
            # find next begin

            end = None
            begin = None
            begin_orig_idx = None
            # if i == 10413:
            #     logging.info('just the one')
            #     logging.info(df_gst.to_string())
            df_filter = df_gst[0:i+1]
            # logging.info('for filter before ')
            # logging.info(df_filter.to_string())
            end_list = df_filter[df_filter['begin_end'] == 'end'].index.to_list()
            if len(end_list) > 0:
                end = end_list[-1]

            df_filter = df_gst[i:]
            # logging.info('for filter after ')
            # logging.info(df_filter.to_string())
            begin_list = df_filter[df_filter['begin_end'] == 'begin'].index.to_list()
            if len(begin_list) > 0:
                begin = begin_list[0]
                begin_orig_idx = df_gst['orig_idx'].iloc[begin]
                # logging.info('begin {} begin_orig_idx {}'.format(begin, begin_orig_idx))

            # logging.info('exit here')
            # logging.info(end)
            # logging.info(begin)
            # exit()

            # logging.info('End begin')
            # logging.info(end)
            # logging.info(begin)

            
            if end is not None and begin is not None:
                # logging.info('A')
                # logging.info(end)
                # logging.info(begin)

                # begin_orig_idx = df_gst['orig_idx'].iloc[begin]

                # try to fix the timestamps 
                # take a copy, which we change if possible 
                df_bad = df_gst.iloc[end:begin+1].copy()
                # logging.info(df_bad.to_string())
                # logging.info('my bad')
                df_bad, df_dropped, rec_adjusted_timestamps1, compliance_passed = station_fix_bad_timestamps(df_bad,df_dropped,df_orig)
                if compliance_passed == False:
                    return df_gst, df_dropped, df_orig, rec_adjusted_timestamps, compliance_passed

                # the index will change 
                rec_adjusted_timestamps += rec_adjusted_timestamps1
                
                df_before = df_gst.iloc[:end]
                df_after = df_gst.iloc[begin+1:]

                # concat back together 
                df_gst = pd.concat([df_before, df_bad, df_after], ignore_index = True)
                df_gst.reset_index(inplace=True,drop=True) 
                # print('EE')
                df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

            else:
                # logging.info('B')
                # logging.info('else')

                if begin is None and end is None:
                    to_drop = list(range(0,len(df_gst)))
                    df_dropped1 = df_gst.iloc[to_drop]
                    df_dropped = pd.concat([df_dropped, df_dropped1], ignore_index = False)
                    print('SEVERE: Ground station/Station - {} {} Dropping the section ({} records)'.format(corr_ground_station, orig_station, len(df_dropped1)))
                    logging.info('SEVERE: Ground station/Station - {} {} Dropping the section ({} records)'.format(corr_ground_station, orig_station, len(df_dropped1)))
                    logging.debug('Dropping the section--------------------------')
                    logging.debug(df_orig.to_string())
                    df_orig_shown = True
                    logging.debug('Dropping the section--------------------------')
                    df_gst.drop(df_gst.index[to_drop],inplace=True)
                    # there's no data left, so break out of the loop
                    break 
                    

                    # raise Exception

                elif begin is None:
                    # logging.info('C')
                    # it wasn't found, so delete to the end 
                    to_drop = list(range(end,len(df_gst)))
                    df_dropped1 = df_gst.iloc[to_drop]
                    df_dropped = pd.concat([df_dropped, df_dropped1], ignore_index = False)
                    logging.info('WARNING - dropping {} record(s)'.format(len(df_dropped1)))
                    df_gst.drop(df_gst.index[to_drop],inplace=True)
#             
#             
                elif end is None:
                    # logging.info('before')
                    # logging.info(df_gst.to_string())
                    # it wasn't found, so delete from the beginning
                    to_drop = list(range(0,begin))
                    df_dropped1 = df_gst.iloc[to_drop]
                    df_dropped = pd.concat([df_dropped, df_dropped1], ignore_index = False)
                    logging.info('WARNING - dropping {} record(s)'.format(len(df_dropped1)))
                    df_gst.drop(df_gst.index[to_drop],inplace=True)
                    df_gst.reset_index(inplace=True,drop=True) 
                    # print('FF')
                    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
                    # logging.info('Probably ok, but need to check')
                    # logging.info('after')
                    # logging.info(df_gst.to_string())
                    # raise Exception

    else: # if no gaps were found, we need a reference for df_orig (an empty copy)
        df_orig = df_gst.drop(df_gst.index).copy()

    if len(df_gst) > 0:

        # fix any gaps
        df_gst = station_fix_missing_timestamps2(df_gst)

        # now make a test to see if it worked
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        # print('FF')
    
    return df_gst, df_dropped, df_orig, rec_adjusted_timestamps, compliance_passed

# check the timeseries data and remove any that don't fit
def check_compliance(df_gst,gzip_filename):

    logging.debug('DEBUG: check_compliance')
    compliance_passed = True


    rec_deleted_timestamps = 0 #meaning changed
    rec_fixed_frame = 0
    rec_fixed_simple_timestamp = 0
    rec_adjusted_timestamps = 0
    rec_deleted_non_uniq = 0

    # if gaps > 0:
    orig_station = df_gst['orig_station'].iloc[0]
    corr_ground_station = df_gst['corr_ground_station'].iloc[0]

    # if config.initial:
    #     # use this setting to view large breaks and gaps 
    #     station_test_missing_timestamps(df_gst,gzip_filename)
    #     return df_gst, None, 0, 0, 0, 0, 0
    # else: 
    #     station_test_missing_timestamps(df_gst,gzip_filename)


    # Fix the frame number if it was incorrectly recorded 
    df_gst, rec_fixed_frame, rec_repeated_frame = station_fix_frames(df_gst)

    # Where possible, do a really simple fix on the timestamp
    df_gst, rec_fixed_simple_timestamp = station_simple_fix_timestamp(df_gst)

    # find the cumulative sum of the good records 
    df_gst = calculate_cumsum(df_gst)
    # 
    # logging.debug('Overall--------------------------')
    # logging.debug(df_gst.to_string())
    # # df_orig_shown = True
    # logging.debug('Overall--------------------------')

    if len(df_gst) > 0:

        # estimate the delta as it changes over time
        df_gst = calculate_delta4(df_gst)

        # make empty copy for df_dropped
        df_dropped = df_gst.drop(df_gst.index).copy()

        # fix timestamps if possible 
        df_gst, df_dropped, df_orig, rec_adjusted_timestamps, compliance_passed = fix_timestamps(df_gst, df_dropped)

        rec_deleted_non_uniq = 0

    return df_gst, df_dropped, df_orig, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps, rec_repeated_frame, compliance_passed

def station_simple_fix_timestamp(df_gst):

    logging.debug('DEBUG: station_simple_fix_timestamp')

    rec_fixed_simple_timestamp = 0 
    # if there are frame gaps, that aren't really gaps, get rid of them 
    simple_timestamp_errors = ((df_gst['corr_gap_count']==-8888) & (df_gst['frame_change']==1)).sum()
    if simple_timestamp_errors > 0:
        orig_station = df_gst['orig_station'].iloc[0]
        corr_ground_station = df_gst['corr_ground_station'].iloc[0]

        idx_list = df_gst[(df_gst['corr_gap_count']==-8888) & (df_gst['frame_change'] == 1)].index.tolist()
        for i in idx_list:
            if i > 0 and i < len(df_gst) -1 : 
                next_frame_change = df_gst['frame_change'].iloc[i+1]
                prev_timestamp = df_gst['corr_timestamp'].iloc[i-1]
                current_timestamp = df_gst['corr_timestamp'].iloc[i]
                next_timestamp = df_gst['corr_timestamp'].iloc[i+1]

                single_gap = (current_timestamp -  prev_timestamp).total_seconds()
                double_gap = (next_timestamp -  prev_timestamp).total_seconds()

                if single_gap > 0.6029 and single_gap < 0.6041:
                    # we can ignore it (probably because it has already been fixed from 
                    # previous record)
                    pass

                elif next_frame_change == 1 and double_gap > 1.2059 and double_gap < 1.2081:

                    new_timestamp = prev_timestamp + pd.Timedelta(seconds=double_gap/2)
                    df_gst.at[i,'corr_timestamp'] = new_timestamp
                    rec_fixed_simple_timestamp += 1

        if rec_fixed_simple_timestamp > 0:
            # calculate the gaps and frame change
            # print('GG')
            df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 
            logging.info('WARNING: Ground station/Station - {} {} Simple adjustments to {} timestamp(s)'.format(
                corr_ground_station, orig_station, rec_fixed_simple_timestamp))

    return df_gst, rec_fixed_simple_timestamp

def station_fix_frames(df_gst):
    # Fix the frame number if it was incorrectly recorded

    logging.debug('DEBUG: station_fix_frames')

    rec_fixed_frame = 0
    rec_repeated_frame = 0

    # if there are frames greater than 89, set to 89
    df_gst['frame'] = np.where((df_gst['frame'] >89), 89, df_gst['frame'])


    # calculate the gaps and frame change
    # print('HH')
    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 


    # logging.info('first calculate gaps ')
    # logging.info(df_gst.to_string())
    # exit()


    # sometimes the frame number has been recorded incorrectly
    # it can be corrected if the gaps are correct, and the difference between the previous and 
    # next is 2
    frame_gaps = ((df_gst['frame_change']!=1) & (df_gst['corr_gap_count'] != df_gst['frame_change'])).sum()
    if frame_gaps > 0:
    
        orig_station = df_gst['orig_station'].iloc[0]
        corr_ground_station = df_gst['corr_ground_station'].iloc[0]
    
        idx_list = df_gst[(df_gst['frame_change']!=1) & (df_gst['corr_gap_count'] != df_gst['frame_change'])].index.tolist()
        for i in idx_list:
            if i > 0 and i < len(df_gst) -1 : 

                # logging.info('frame stuff')
                # logging.info(df_gst['orig_idx'].iloc[i])
                prev_frame = df_gst['frame'].iloc[i-1]
                curr_frame = df_gst['frame'].iloc[i]
                next_frame = df_gst['frame'].iloc[i+1]    
                corr_gap = df_gst['corr_gap'].iloc[i]
                frame_change = df_gst['frame_change'].iloc[i]
                prev_corr_gap_count = df_gst['corr_gap'].iloc[i-1]

                if (corr_gap > 0.6029 and corr_gap < 0.6041 and add_or_minus_frame(prev_frame,2) == next_frame):
                    # first check that it's not already been fixed 
                    if add_or_minus_frame(prev_frame,1) != curr_frame:
                        df_gst.at[i,'frame'] = add_or_minus_frame(prev_frame,1)
                        rec_fixed_frame += 1


        if rec_fixed_frame > 0:
            # print('II')
            df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
            logging.info('WARNING: Ground station/Station - {} {} Simple adjustments to {} frame number(s)'.format(
                corr_ground_station, orig_station, rec_fixed_frame))

    idx_list = df_gst[(df_gst['frame_change']==0) & (df_gst['corr_gap_count']==-8888) ].index.tolist()
    to_drop = []
    for i in idx_list:
        if i > 0 and i < len(df_gst) -1 : 
            # prev_frame = df_gst['frame'].iloc[i-1]
            # curr_frame = df_gst['frame'].iloc[i]
            # next_frame = df_gst['frame'].iloc[i+1]    
            # corr_gap = df_gst['corr_gap'].iloc[i]
            # frame_change = df_gst['frame_change'].iloc[i]
            # prev_corr_gap_count = df_gst['corr_gap'].iloc[i-1]
    
            # if the two frames are equal, then it is a repeated frame
            # (the case where the frames are equal, but they are 90 
            # timestamps apart is caught in calculate_gaps())
    
            # delete the PREVIOUS frame 
            to_drop.append(i-1)
            rec_repeated_frame += 1
    
    if len(to_drop) > 0:
        df_gst.drop(df_gst.index[to_drop],inplace=True)
        df_gst.reset_index(inplace=True,drop=True) 
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)


        # 
        # 
        # # finally, if there are still some left, delete records that are next to each other
        # frame_gaps = ((df_gst['frame_change']>1) & (df_gst['corr_gap_count'] != df_gst['frame_change'])).sum()
        # 
        # if frame_gaps > 0:
        #     to_drop2 = []
        #     idx_list = df_gst[(df_gst['frame_change']>1) & (df_gst['corr_gap_count'] != df_gst['frame_change'])].index.tolist()
        #     for a, idx in enumerate(idx_list[1:]):
        #         # check for consecutive entries
        #         prev_idx = idx_list[a]
        #         if idx-1 == prev_idx:
        #             # if they are equal, remove both
        #             to_drop2.append(idx)
        #             rec_repeated_frame += 1
        #             if prev_idx not in to_drop2:
        #                 to_drop2.append(prev_idx)
        #                 rec_repeated_frame += 1
        #     if len(to_drop2) > 0:
        #         df_gst.drop(df_gst.index[to_drop2],inplace=True)
        #         df_gst.reset_index(inplace=True,drop=True)
        #         df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        # 
        # if rec_repeated_frame > 0:
        #     logging.info('WARNING: Ground station/Station - {} {} Removed {} repeated frames and repeated gaps'.format(
        #         corr_ground_station, orig_station, rec_repeated_frame)) 
        # 
        # df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)



    return df_gst, rec_fixed_frame, rec_repeated_frame

def station_test_missing_timestamps(df_gst,gzip_filename):
    '''
    Test for large gaps and negative gaps which might need resets 

    '''
    logging.debug('DEBUG: station_test_missing_timestamps')

    # print('JJ')
    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 

    # gaps_very_long = (((df_gst['corr_gap']>53.1344*10) & (df_gst['corr_gap_count']!=-8888)) | (df_gst['corr_gap']<0)).sum()
    # warn about both long gaps and negative gaps
    gaps_very_long = ((df_gst['corr_gap']>53.1344*10) | (df_gst['corr_gap']<0)).sum()
    if gaps_very_long > 0:
        orig_station = df_gst['orig_station'].iloc[0]
        orig_ground_station = df_gst['orig_ground_station'].iloc[0]
        idx_list = df_gst[(df_gst['corr_gap']>53.1344*10) | (df_gst['corr_gap']<0)].index.tolist()
        problem_gzip_filename = os.path.basename(gzip_filename)
        for i in idx_list:
            if df_gst['corr_gap'].iloc[i] < 0:

                corr_timestamp1 = df_gst['corr_timestamp'].iloc[i]
                corr_timestamp_low = corr_timestamp1 - pd.Timedelta(seconds=1)
                corr_timestamp_high = corr_timestamp1 + pd.Timedelta(seconds=1)
                frame1 = df_gst['frame'].iloc[i]
                orig_idx = df_gst['orig_idx'].iloc[i]

                logging.debug('WARNING: Neg Gap {:0.02f}, Corr_gap {} {} {}'.format(df_gst['corr_gap'].iloc[i],df_gst['corr_gap_count'].iloc[i],corr_timestamp_low,corr_timestamp_high))

                idx_list_dup = df_gst[(df_gst['corr_timestamp'] > corr_timestamp_low) & (df_gst['corr_timestamp'] < corr_timestamp_high) & (df_gst['frame'] == frame1)].index.tolist()
                # check for a duplicate, which means that it is probably recording to the wrong ground station
                if len(idx_list_dup) >= 2:

                    # note that we use the orig_ground_station, because
                    # when this code is used, we need to use the 
                    # original one, not the corrected one
                    # pd.to_datetime(config.last_timestamp)
                    logging.debug('WARNING: Negative Gap {:0.02f}, Corr_gap {}'.format(df_gst['corr_gap'].iloc[i],df_gst['corr_gap_count'].iloc[i]))
                    # logging.debug('''df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='{}',start_timestamp ='{}',end_timestamp='{}',start_station='{}',end_station='{}',start_orig_no={},end_orig_no={},orig_ground_station1={})'''.format(problem_gzip_filename, pd.to_datetime(corr_timestamp1),pd.to_datetime(config.last_timestamp),orig_station,config.last_station,orig_no,config.last_orig_no,orig_ground_station))
                    logging.debug('''df = reset_all_ground_stations_idx(df, gzip_filename, problem_gzip_filename='{}',orig_idx_start ={},orig_idx_end={},single_station=None)'''.format(problem_gzip_filename,orig_idx,config.last_orig_idx))
                    # YYYY
            else:
                # logging.debug('WARNING: Large Gap {:0.02f}, Corr_gap {}'.format(df_gst['corr_gap'].iloc[i],df_gst['corr_gap_count'].iloc[i]))
                logging.debug('WARNING: Large Gap\n{}'.format(df_gst.iloc[i-1:i+1].to_string()))


def station_fix_missing_timestamps2(df_bad):
    '''
    Search for gaps in the record. Blank records are inserted if the 
    following conditions are met:
    Can either be called by the bad record, or for the whole 
    dataframe.
    TODO - what conditions 
    
    '''

    logging.debug('DEBUG: station_fix_missing_timestamps2')

    # gaps_long = (df_bad['corr_gap_count']>1).sum()
    # if gaps_long > 0:
        

    # first make a new index that we can maniuplate
    df_bad['index1'] = df_bad.index

    # gaps_largest = df_gst['corr_gap'].max()

    # logging.info('before bad')
    # logging.info(df_bad.to_string())

    # YYYY

    rec_missing = 0
    # look for large gaps 
    idx_list = df_bad[(df_bad['corr_gap_count'] > 1) | (df_bad['corr_gap_count'] == -9999)].index.tolist()
    # logging.info('errors found ')
    # for i in idx_list:
    #     logging.info(df_bad.iloc[i:i+1].to_string())

    for i in idx_list:

        orig_station = df_bad['orig_station'].iloc[0]
        corr_ground_station = df_bad['corr_ground_station'].iloc[0]
        orig_ground_station = df_bad['orig_ground_station'].iloc[0]

        total_gap = df_bad['corr_gap'].iloc[i]  
        corr_gap_count = df_bad['corr_gap_count'].iloc[i]
        # if corr_gap_count == -9999:
        #     frame_gap = df_bad['guess_count'].iloc[i]
        # else:
        frame_gap = corr_gap_count
        # single_gap = total_gap/corr_gap_count
        single_gap = total_gap/frame_gap
        # single_gap = 

        new_i = i-1 
        corr_timestamp=df_bad['corr_timestamp'].iloc[i-1]
        current_frame = df_bad['frame'].iloc[i-1]
        # current_actual_frame_gap = df_bad['actual_frame_gap'].iloc[i-1]
        # sample_no = df_bad['sample_no'].iloc[i]
        dict_list = []

        # print(frame_change)
        for sample_i in range(0, frame_gap-1):
            # insert new records with index numbers which are evenly split 
            # between the previous record and the next
            missing_interval = 1/frame_gap
            new_i += missing_interval
            
            new_timestamp = corr_timestamp + pd.Timedelta(seconds=(sample_i+1)*single_gap)
            if corr_gap_count == -9999:
                current_frame = pd.NA
            else:
                current_frame = add_or_minus_frame(current_frame,1)
            # current_actual_frame_gap += 1
            
            dict = {
                'index1' : new_i, 
                'orig_no' : pd.NA, 
                'bit_synchronizer' : pd.NA, 
                'orig_timestamp' : pd.NA, 
                'corr_timestamp' : new_timestamp, 
                'orig_frame' : pd.NA, 
                'orig_station' : orig_station, 
                'clock_flag' : pd.NA,
                'orig_ground_station' : orig_ground_station, 
                'orig_mh1_1' : pd.NA, 
                'orig_mh2_1' : pd.NA, 
                'orig_mhz_1' : pd.NA, 
                'orig_mh1_2' : pd.NA, 
                'orig_mh2_2' : pd.NA, 
                'orig_mhz_2' : pd.NA, 
                'orig_mh1_3' : pd.NA, 
                'orig_mh2_3' : pd.NA, 
                'orig_mhz_3' : pd.NA, 
                'orig_mh1_4' : pd.NA, 
                'orig_mh2_4' : pd.NA, 
                'orig_mhz_4' : pd.NA, 
                'sync' : '1110001001000011101101',
                'frame' : current_frame,
                'corr_gap' : single_gap,
                'frame_change' : 1,
                'corr_gap_count' : 1,
                'corr_ground_station' : corr_ground_station,
                # 'actual_frame_gap' : current_actual_frame_gap
                }

            dict_list.append(dict)

            rec_missing += 1

        df_append = pd.DataFrame.from_dict(dict_list)
        # logging.info('after 1\n{}'.format(df_append.to_string()))
        df_bad = df_bad.append(df_append)

    if rec_missing>0:
        # note that we do this at the end
        # print('a')
        df_bad.set_index('index1', inplace=True)
        df_bad = df_bad.sort_index().reset_index(drop=True)
        # logging.info('WARNING: Ground station/Station - {} {} Inserted {} blank record(s)'.format(
        #     corr_ground_station, orig_station, rec_missing))
        # print('KK')
        df_bad, gaps_long, gaps_8888 = calculate_gaps(df_bad)

    # logging.info('after bad\n{}'.format(df_bad.to_string()))

    return df_bad


def add_or_minus_frame(frame,add=1):
    if add < 0:
        new_frame = frame + add
        # print('first guess ', new_frame)

        records = add*-1 // 90
        # print('records ', records)
        # new_frame = frame + add
        # print('first guess ', new_frame)
        if new_frame < 0:
            new_frame += records*90
        # print('new frame ', new_frame)
    elif add == 0:
        new_frame = frame
    elif add > 0:
        new_frame = frame + add
        # print('first guess ', new_frame)
        if new_frame > 89:
            full_frames = new_frame // 90
            # print('full_frames ', full_frames)
            new_frame -= (full_frames)*90

    return new_frame

def minus_frame(old_frame,new_frame):
    diff = new_frame - old_frame
    if diff < 0:
        diff += 90
    return diff

# def station_remove_damaged_timestamps(df_gst):
#     # delete records with gaps that are outside the tolerance
# 
#     rec_damaged_outside_tolerance = 0
# 
#     gaps_8888 = (df_gst['corr_gap_count']==-8888).sum()
# 
#     if gaps_8888 > 0: 
# 
#         # if this happens, raise an error? 
#         # logging.info('Raising an error because we need to remove damaged timestamps')
# 
#         orig_station = df_gst['orig_station'].iloc[0]
#         corr_ground_station = df_gst['corr_ground_station'].iloc[0]
# 
#         idx_list = df_gst[(df_gst['corr_gap_count'] != 1)].index.tolist()
# 
#         gaps_8888 = (df_gst['corr_gap_count']==-8888).sum()
# 
#         if gaps_8888 == len(df_gst) - 1:
#             # mark last record for deletion if only one record left
#             df_gst.at[0,'corr_gap_count'] = -8888
#             gaps_8888 += 1
# 
#         if gaps_8888 > 0:
#             # idx_list = df_gst[df_gst['corr_gap_count'] == -8888].index.tolist()
#             # for i in idx_list:
#             #     logging.debug('BEFORE\n{}'.format(df_gst.iloc[i-5:i+5].to_string()))
# 
#             logging.info('WARNING: Ground station/Station - {} {} Removing {} damaged record(s) - time differences outside normal tolerance or records with many gaps'.format(
#                 corr_ground_station, orig_station, gaps_8888))
#             # drop damaged records
#             df_gst.drop( df_gst[ df_gst['corr_gap_count'] == -8888].index , inplace=True)
# 
#             # reindex
#             df_gst.reset_index(inplace=True,drop=True)
# 
#             rec_damaged_outside_tolerance = gaps_8888
# 
#             # recalculate the gaps (this must be done at the end)
#             df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)    
# 
#     return df_gst, rec_damaged_outside_tolerance

def station_final_check(df_gst,df_orig):
    # Make final checks to check everything is OK with the timestamps

    # recalculate the gaps
    # print('LL')
    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

    # rec_damaged_outside_tolerance = (df_gst['tolerant_gap_select']==False).sum()

    orig_station = df_gst['orig_station'].iloc[0]
    corr_ground_station = df_gst['corr_ground_station'].iloc[0]

    global df_orig_shown

    # gaps in frame seq

    gaps_7777 = (df_gst['corr_gap_count']==-7777).sum()
    if gaps_7777 > 0: 
        logging.info('SEVERE: Ground station/Station - {} {} {} frames are reset to zero (-7777 error)'.format(corr_ground_station, orig_station, gaps_7777))
        print('SEVERE: Ground station/Station - {} {} {} frames are reset to zero (-7777 error)'.format(corr_ground_station, orig_station, gaps_7777))
        # raise Exception

    # idx_list = df_gst[~df_gst['corr_gap_count'].isin([1,-7777])].index.tolist()
    idx_list = df_gst[~((df_gst['frame_change'] == 1) | (df_gst['corr_gap_count'] == -7777) | (df_gst['frame_change'].isna()))].index.tolist()
    
    if len(idx_list) > 0:
        
        logging.info('EXCEPTION: Final Check: Ground station/Station - {} {} Gaps in frame sequence: {}'.format(
            corr_ground_station, orig_station, len(idx_list)))
        print('EXCEPTION: Final Check: Ground station/Station - {} {} Gaps in frame sequence: {}'.format(
            corr_ground_station, orig_station, len(idx_list)))
        for i in idx_list:
            logging.info(df_gst.iloc[i-1:i+2].to_string())
        if df_orig_shown == False:
            logging.debug('--------------------------')
            logging.debug(df_orig.to_string())
            df_orig_shown = True
            logging.debug('--------------------------')


    # idx_list = df_gst[~df_gst['guess_count'].isna()].index.tolist()
    # if len(idx_list) > 0:
    # 
    #     logging.info('SEVERE: Final Check: Ground station/Station - {} {} Some guesses made: {}'.format(
    #         corr_ground_station, orig_station, len(idx_list)))
    #     print('SEVERE: Final Check: Ground station/Station - {} {} Some guesses made: {}'.format(
    #         corr_ground_station, orig_station, len(idx_list)))
    #     for i in idx_list:
    #         logging.info(df_gst.iloc[i-1:i+2].to_string())
    #     if df_orig_shown == False:
    #         logging.debug('--------------------------')
    #         logging.debug(df_orig.to_string())
    #         df_orig_shown = True
    #         logging.debug('--------------------------')


    idx_list = df_gst[df_gst['frame'].isna()].index.tolist()
    if len(idx_list) > 0:
        logging.info('WARNING: Final Check: Ground station/Station - {} {} Record(s) without frame numbers: {}'.format(
            corr_ground_station, orig_station, len(idx_list)))
        logging.info(df_gst.iloc[idx_list].to_string())

    rec_adjusted_timestamps_new = ((df_gst['orig_timestamp'] != df_gst['corr_timestamp']) & (~df_gst['orig_no'].isna())).sum()

    # logging.info('final check')
    # logging.info(df_gst.to_string())

    # if gaps_8888 > 0: 
    # 
    #     idx_list = df_gst[df_gst['corr_gap_count']==-8888].index.tolist()
    #     # idx_list = df_gst[df_gst['tolerant_gap_select'] == False].index.tolist()
    #     # logging.info(idx_list)
    #     for i in idx_list:
    #         logging.debug('WARNING: Final Check - Record outside of tolerance:\n{}'.format(df_gst.iloc[i-1:i+2].to_string()))
    # 
    #     # logging.info(df_gst.to_string())
    # 
    #     # logging.info('EXCEPTION: Final Check: Ground station/Station - {} {} Outside of tolerance: {}'.format(
    #     #     corr_ground_station, orig_station, gaps_8888))
    #     # # logging.info('nasty end\n{}'.format(df_gst.to_string()))
    #     # raise Exception
    # 
    # # logging.info('The End{}\n{}'.format(len(df_gst),df_gst.to_string()))
    # 
    # if gaps_long > 0: 
    # 
    #     idx_list = df_gst[df_gst['corr_gap_count']>1].index.tolist()
    #     # idx_list = df_gst[df_gst['tolerant_gap_select'] == False].index.tolist()
    #     # logging.info(idx_list)
    #     for i in idx_list:
    #         logging.debug('WARNING: Final Check - Record outside of tolerance:\n{}'.format(df_gst.iloc[i-1:i+1].to_string()))
    # 
    #     logging.info('EXCEPTION: Final Check: Ground station/Station - {} {} Unable to insert gaps: {}'.format(
    #         corr_ground_station, orig_station, gaps_long))
    #     raise Exception
    # 
    # rec_negative = (df_gst['corr_gap']<0).sum()
    # if rec_negative > 0:
    #     logging.info('EXCEPTION: Final Check: Ground station/Station - {} {} Negative gaps: {}'.format(
    #         corr_ground_station, orig_station, rec_negative))
    #     # logging.info('not stopping right now')
    #     raise Exception
    # 
    # # re-calculate whether there are any gaps in the frame sequence
    # df_gst['actual_frame_gap_diff'] = df_gst['actual_frame_gap'] - df_gst['actual_frame_gap'].shift(1)
    # df_gst.at[0,'actual_frame_gap_diff'] = 1
    # idx_list = df_gst[df_gst['actual_frame_gap_diff'] != 1].index.tolist()
    # if len(idx_list) > 0:
    #     for idx in idx_list:
    #         logging.info(df_gst.iloc[idx:idx+1].to_string())
    
    # print(df_gst['actual_frame_gap_diff'].iloc[0:10])
    # not_monotonic = (diff != 1).sum() 
    # print(not_monotonic)
    # if not_monotonic > 0:
    #     logging.info(not_monotonic)
        # logging.info('EXCEPTION: Final Check: Ground station/Station - {} {} Gaps in frame sequence: {}'.format(
        #     corr_ground_station, orig_station, rec_negative))
    #     logging.info('EXCEPTION: Final Check: Ground station/Station - {} {} Gaps in frame sequence: {}'.format(
    #         corr_ground_station, orig_station, len(idx_list)))
    # #     # logging.info('not stopping right now')
    #     raise Exception 


    return df_gst

def calculate_frame_change(df_gst):

    df_gst['frame_change'] = df_gst['frame'] - df_gst['frame'].shift(1)  
    df_gst.at[0,'frame_change'] = 1
    df_gst['frame_change'] = np.where((df_gst['frame_change'] <0), df_gst['frame_change']+90, df_gst['frame_change'])


    # this special case includes where the timing has been guessed to put the data
    # back together
    # df_gst['frame_change'] = np.where((df_gst['frame'].isna()), 1, df_gst['frame_change'])

    return df_gst

def calculate_gaps(df_gst):

    gaps_long = 0
    gaps_8888 = 0

    if len(df_gst) > 0:

        corr_gap1= df_gst['corr_gap'].iloc[0]
        frame_change1= df_gst['frame_change'].iloc[0]
        corr_gap_count1= df_gst['corr_gap_count'].iloc[0]

        # df_gst['frame'] = to_Int64(df_gst['frame'])
        # Calcuate the gaps between timestamps, and create boolean columns which 
        # determine which values are within the tolerance

        # calculate the corrected gap based on the corrected timestamp column
        df_gst['corr_gap'] = df_gst['corr_timestamp'].diff().dt.total_seconds()

        # calculate which records are within the tolerance 
        # df_gst['tolerant_gap'] = np.where((df_gst['corr_gap'] > 0.5999) & (df_gst['corr_gap'] < 0.6081), True, False)
        # Set the first record to True (it is currently False, because there)
        # is nothing to compare it to)
        # df_gst.at[0,'tolerant_gap'] = True

        # tolerant_gap_select is the same as tolerant_gap, except that gaps 
        # which are longer than 60 s are excluded
        # df_gst['tolerant_gap_select'] = np.where((df_gst['corr_gap'] >= 60.000), True, df_gst['tolerant_gap'])

        # logging.info('calculate FRAME{}'.format(type(df_gst['frame'].iloc[20])))

        df_gst = calculate_frame_change(df_gst)
        # logging.info('frame ')
        
        # default for correct gap count 
        df_gst['corr_gap_count'] = df_gst['frame_change']
        # logging.info('frame 1')

        # check most of the records
        # df_gst['corr_gap_count'] = np.where((df_gst['frame_change']) == 1 & ((df_gst['corr_gap'] <  0.6009) |  (df_gst['corr_gap'] <  0.6071)), df_gst['corr_gap_count'],  -8888) 
        df_gst['corr_gap_count'] = np.where((df_gst['frame_change']) == 1 & (df_gst['corr_gap'] <  config.high_tolerance), df_gst['corr_gap_count'],  -8888)

        # small and negative gaps, error = -8888
        df_gst['corr_gap_count'] = np.where((df_gst['corr_gap'] <= config.low_tolerance ) , -8888, df_gst['corr_gap_count'])

        # make checks to see whether the gaps larger than 1 are correct 
        idx_list = df_gst[df_gst['frame_change'] > 1].index.tolist()
        for i in idx_list:
            single_gap = df_gst['corr_gap'].iloc[i]/df_gst['frame_change'].iloc[i]
            if single_gap >= config.low_tolerance and single_gap <= config.high_tolerance: 
                df_gst.at[i,'corr_gap_count'] = df_gst['frame_change'].iloc[i]
            # special test for weird condition when 
            # the frame resets itself
            frame = df_gst['frame'].iloc[i]
            corr_gap = df_gst['corr_gap'].iloc[i]
            prev_cumsum = df_gst['cumsum'].iloc[i-1]

            # mark as a possible frame reset (-7777 error)
            if frame == 0 and corr_gap > 0.60 and corr_gap < 0.608 and prev_cumsum > 5:
                try: 

                    # check that the next gap is also OK
                    next_corr_gap = df_gst['corr_gap'].iloc[i+1]
                    if next_corr_gap > 0.60 and next_corr_gap < 0.608:    
                        # check that the next frame is 1 
                        next_frame = df_gst['frame'].iloc[i+1]
                        if next_frame == 1:
                            df_gst.at[i,'corr_gap_count'] = -7777
                except:
                    pass

        # separate checks if there are large gaps (but smaller than four hours)
        idx_list = df_gst[(df_gst['corr_gap'] >= 53.1344) & (df_gst['corr_gap'] < 14400)].index.tolist()
        for i in idx_list:
            # logging.info(df_gst.iloc[i:i+1].to_string())
            corr_gap = df_gst['corr_gap'].iloc[i]
            last_frame = df_gst['frame'].iloc[i-1]
            new_frame = df_gst['frame'].iloc[i]
            # prev_cumsum = df_gst['cumsum'].iloc[i-1]
            # if the previous cumsum was high, this may be a real gap 
            # note, can only check this once the cumsum has be set, before that,
            # all the gaps will be set as -8888
            # if prev_cumsum > config.cumsum_test: 
            # delta4_mean = delta4_arr[i:valid_idx+1].mean()

            delta4_mean = df_gst['delta4'].iloc[i-1:i+1].mean()
            # logging.info('delta4_mean')            
            # logging.info(delta4_mean)
            if math.isnan(delta4_mean):
                delta4_mean = DELTA*4

            # logging.info(df_gst.to_string())
            # logging.info(delta4_mean)

            # logging.info('corr gap count in the guess')
            # logging.info(df_gst['corr_gap_count'].iloc[i])
            
            single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(last_frame,new_frame,corr_gap,delta4=delta4_mean,frame_correction=0)
                
            # if actual_frame_gap == 1042:
            #     print('weird!!')
            #     logging.info('weird!!!')
            #     single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(last_frame,new_frame,corr_gap,delta4=delta4_mean,frame_correction=-10)
            # if actual_frame_gap == 529:
            #     print('weird!!')
            #     logging.info('weird!!!')
                # single_gap, actual_frame_gap, full_frames_est, time_error, percent_error = frame_diff_positive(last_frame,new_frame,corr_gap,delta4=delta4_mean,frame_correction=-3)
            if single_gap > config.low_single_gap and single_gap < config.high_single_gap:
                df_gst.at[i,'corr_gap_count'] = actual_frame_gap


            # # if there's more than 120 s between them, make a guess 
            # elif corr_gap > 120.0:
            #     # make some kind of guess 
            # 
            #     # logging.info(df_gst.iloc[i-1:i+2].to_string())
            #     frame_correction = round(est_frame_gap) - actual_frame_gap
            #     single_gap2, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(last_frame,new_frame,corr_gap,delta4=delta4_mean,frame_correction=frame_correction)
            #     df_gst.at[i,'corr_gap_count'] = -9999
            #     # logging.info('WARNING: Guessing')  
            #     df_gst.at[i,'guess_count'] = actual_frame_gap  
            #     # logging.info('old single gap = {}, new single gap = {}'.format(single_gap,single_gap2))
            #     # 
            #     # logging.info('Fixed guess')
            #     # logging.info(df_gst.iloc[i-1:i+2].to_string())


        # the first record should be set to what it was before 
        # df_gst.at[0,'corr_gap'] = corr_gap1
        # df_gst.at[0,'frame_change'] = frame_change1
        # df_gst.at[0,'corr_gap_count'] = corr_gap_count1
        # the first record should be reset
        df_gst.at[0,'corr_gap'] = pd.NA
        df_gst.at[0,'frame_change'] = 1
        df_gst.at[0,'corr_gap_count'] = 1

        # logging.info('just before')
        # logging.info(df_gst.to_string())
        gaps_long = (df_gst['corr_gap_count']>1).sum()
        gaps_8888 = (df_gst['corr_gap_count']==-8888).sum()
    
    return df_gst, gaps_long, gaps_8888

def frame_diff_positive(last_frame,new_frame,gap,delta4,frame_correction=0):
    # get the estimate of the number of frames

    


    actual_frame_gap = pd.NA
    full_frames_est = pd.NA
    single_gap = pd.NA
    time_error = pd.NA
    percent_error = pd.NA

    if gap < 0:
        neg_gap = True
        gap *= -1
    else:
        neg_gap = False

    # get an approximate number for the frame gap 
    est_frame_gap = gap/0.6038

    # get an approximate number for the full number of frames
    full_frames_est1 = int(est_frame_gap // 90)

    # if there's a frame correction, add it to new frame 
    if frame_correction == 0:
        corrected_new_frame = new_frame
    else:
        corrected_new_frame = add_or_minus_frame(new_frame,frame_correction)

    # find difference between the last frame and the new frame 
    # (this number will ignore the full frames)
    diff = minus_frame(old_frame=last_frame,new_frame=corrected_new_frame)

    # full_frames_est1 can be wrong if we are very close to 
    # dividing by 90 frames, so we check the esimate minus and 
    # plus one, and then find the value of the single gap closest
    # to the nominal inteval value 
    # actual_frame_gap is the total number of frames in the gap (that 
    # fits with the last_frame and new_frame)

    test_frame_gap_normal = (full_frames_est1)*90 + diff
    test_frame_gap_minus = (full_frames_est1-1)*90 + diff    
    test_frame_gap_plus = (full_frames_est1+1)*90 + diff

    if test_frame_gap_normal > 0: 
        single_gap_normal = gap/test_frame_gap_normal
    else: 
        single_gap_normal = -999999
    if test_frame_gap_minus > 0:
        single_gap_minus = gap/test_frame_gap_minus
    else:
        single_gap_minus = -999999
    if test_frame_gap_plus > 0:
        single_gap_plus = gap/test_frame_gap_plus 
    else:
        single_gap_plus = -999999

    
    # now work out which is closest to delta4
    if abs(delta4 - single_gap_plus) >  abs(delta4 - single_gap_normal):
        if abs(delta4 - single_gap_normal) > abs(delta4 - single_gap_minus):
            type = 'minus'
        else:
            type = 'normal'
    else: 
        if abs(delta4 - single_gap_plus) > abs(delta4 - single_gap_minus):
            type = 'minus'
        else:
            type = 'plus'

    if type == 'minus':
        full_frames_est = full_frames_est1 - 1
        single_gap = single_gap_minus
    elif type == 'normal':
        full_frames_est = full_frames_est1
        single_gap = single_gap_normal
    elif type == 'plus':
        full_frames_est = full_frames_est1 + 1
        single_gap = single_gap_plus

    actual_frame_gap = (full_frames_est)*90 + diff
    
    # actual_frame_gap_plus = (full_frames_est1+1)*90 + diff
    # full_frames_est = full_frames_est1 + 1
    # 
    # single_gap_plus = gap/actual_frame_gap_plus
    # single_gap = single_gap_plus
    # actual_frame_gap = actual_frame_gap_plus
    # 
    # actual_frame_gap_normal = full_frames_est1*90 + diff
    # if actual_frame_gap_normal != 0: 
    #     single_gap_normal = gap/actual_frame_gap_normal
    #     if abs(single_gap_normal) <  abs(single_gap):
    # 
    #         single_gap = single_gap_normal
    #         actual_frame_gap = actual_frame_gap_normal
    #         full_frames_est = full_frames_est1
    # 
    # if full_frames_est1-1 >= 0:
    #     actual_frame_gap_minus = (full_frames_est1-1)*90 + diff
    #     if actual_frame_gap_minus != 0: 
    #         single_gap_minus = gap/actual_frame_gap_minus
    #         if abs(single_gap_minus) <  abs(single_gap):
    #             single_gap = single_gap_minus
    #             actual_frame_gap = actual_frame_gap_minus
    #             full_frames_est = full_frames_est1 - 1


    # logging.info('actual_frame_gap {} ,est_frame_gap {},diff {},full_frames_est {},single_gap {}, gap={} '.format(actual_frame_gap,est_frame_gap,diff,full_frames_est,gap/actual_frame_gap,gap))

    # if single_gap > 0.6009 and single_gap < 0.6071: 
    #     fit = -4444
    # else:
    #     fit = -3333
    if neg_gap:
        gap *= -1

    # diff1 = actual_frame_gap - (full_frames_est*90)
    time_error = (actual_frame_gap * delta4) - gap
    percent_error = abs(100* (gap-(actual_frame_gap * delta4))/gap)
        

    #     if percent_error > 0.02:
    #         if abs(time_error) < 0.003:
    #             fit = -4444
    #         else: 
    #             fit = -3333
    #     else:
    #         fit = -4444
    # 
    # else: 
    #     fit = pd.NA


    return single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap

    #     logging.info('actual_frame_gap {} ,est_frame_gap {},diff {},a {},single_gap {} '.format(actual_frame_gap,est_frame_gap,diff,a,gap/actual_frame_gap))
    #     # one extra check 
    #     if actual_frame_gap == 0:
    #         frame_change = -8888
    #     elif abs((actual_frame_gap - est_frame_gap)/actual_frame_gap) < 0.001:
    #         single_gap = gap/actual_frame_gap
    #         if single_gap < 0.6009 or single_gap > 0.6061: 
    #             frame_change = -8888
    #         else:
    #             frame_change = actual_frame_gap
    #     else:
    #         frame_change = -8888
    # else: 
    #     frame_change = -8888



def frame_diff(last_frame,new_frame,gap):
    # get the estimate of the number of frames
    if gap > 0: 
        # get an approximate number for the frame gap 
        est_frame_gap = gap/0.6038
        # get an approximate number for the full number of frames
        full_frames_est = int(est_frame_gap // 90)

        # find difference between the last frame and the new frame 
        # (this number will ignore the full frames)
        diff = minus_frame(old_frame=last_frame,new_frame=new_frame)

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
            if abs(single_gap_normal - (DELTA*4)) <  abs(single_gap - (DELTA*4)):

                single_gap = single_gap_normal
                actual_frame_gap = actual_frame_gap_normal

        if full_frames_est-1 >= 0:
            actual_frame_gap_minus = (full_frames_est-1)*90 + diff
            if actual_frame_gap_minus != 0: 
                single_gap_minus = gap/actual_frame_gap_minus
                if abs(single_gap_minus - (DELTA*4)) <  abs(single_gap - (DELTA*4)):
                    single_gap = single_gap_minus
                    actual_frame_gap = actual_frame_gap_minus

        # logging.info('actual_frame_gap {} ,est_frame_gap {},diff {},full_frames_est {},single_gap {} '.format(actual_frame_gap,est_frame_gap,diff,full_frames_est,gap/actual_frame_gap))

        if single_gap > 0.6009 and single_gap < 0.6071: 
            frame_change = actual_frame_gap
        else:
            frame_change = -8888
    else: 
        frame_change = -8888

    #     logging.info('actual_frame_gap {} ,est_frame_gap {},diff {},a {},single_gap {} '.format(actual_frame_gap,est_frame_gap,diff,a,gap/actual_frame_gap))
    #     # one extra check 
    #     if actual_frame_gap == 0:
    #         frame_change = -8888
    #     elif abs((actual_frame_gap - est_frame_gap)/actual_frame_gap) < 0.001:
    #         single_gap = gap/actual_frame_gap
    #         if single_gap < 0.6009 or single_gap > 0.6061: 
    #             frame_change = -8888
    #         else:
    #             frame_change = actual_frame_gap
    #     else:
    #         frame_change = -8888
    # else: 
    #     frame_change = -8888

    return frame_change


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