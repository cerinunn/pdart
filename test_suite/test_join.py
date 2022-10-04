#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
Test Suite 
Test clean is designed to test the join code - 
    and is probably better than this code. 
    Probably obsolete? 
:copyright:
    The pdart Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

import inspect
import io
import os
import pickle
import platform
import unittest
import warnings
from copy import deepcopy
import random
import pandas as pd
import shutil
import glob
# from decimal import Decimal, getcontext TODO Remove 
from pandas.testing import assert_frame_equal
import logging, logging.handlers
import datetime

import numpy as np

from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, read
from pdart.extra_plots.plot_timing_divergence import plot_timing

import matplotlib 
import matplotlib.pyplot as plt
# from pdart.csv_check_work_tapes import (
#   all_drop_duplicates,
#   all_split, 
#     # station_fix_timestamps, 
#   station_fix_missing_timestamps, 
#   initial_cleanup, detailed_report,
#   station_remove_damaged_timestamps,calculate_gaps,
#   station_fix_frames,
#   # station_fix_timeskips, 
#   # add_frame, 
#   minus_frame, frame_diff,
#   station_consec_frames,station_simple_fix_timestamp,
#   get_new_ground_station,all_drop_station_duplicates,all_station_duplicates_v2,
#   all_flat_seismogram, to_Int64
# )
import pdart.config as config
# from pdart.test_suite.test_import import basic_timeseries_station as basic_timeseries_station
# from pdart.test_suite.test_import import default_config as default_config
from pdart.csv_join_work_tapes import call_csv_join_work_tapes
# from pdart.csv_check_work_tapes import calculate_gaps, to_Int64, 

# # consecutive stations (the order is important)
# STATIONS = ['S12', 'S15', 'S16', 'S14', 'S17']
# # consecutive frames
# ORIG_FRAMES = list(range(1,61))
# 
# FRAMES = list(range(0,90))

# nominal delta (sampling rate is 6.25 samples/s)
DELTA = 0.1509433962
INTERVAL = '603774us'

# default date format 
DATEFORMAT='%Y-%m-%dT%H:%M:%S.%fZ'

GZIP_FILENAME='test.csv.gz'

'''
A testcase is created by subclassing unittest.TestCase. The three individual 
tests are defined with methods whose names start with the letters test. This 
naming convention informs the test runner about which methods represent tests.

The crux of each test is a call to assertEqual() to check for an expected 
result; assertTrue() or assertFalse() to verify a condition; or assertRaises() 
to verify that a specific exception gets raised. These methods are used instead 
of the assert statement so the test runner can accumulate all test results and 
produce a report.

The setUp() and tearDown() methods allow you to define instructions that will 
be executed before and after each test method.
'''


class ImportTestCase(unittest.TestCase):
    """
    Test suite for csv_import_work_tapes.py.
    """

    def assertDataFrameEqual(self, a, b, msg):
        try:
            pd.testing.assert_frame_equal(a, b, check_dtype=False)
        except AssertionError as e:
            try:
                pd.testing.assert_frame_equal(a, b, check_dtype=True)
            except AssertionError as e:
                raise self.failureException(msg) from e

    def assertSeriesEqual(self, a, b, msg):
        try:
            pd.testing.assert_series_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        # set up to use pandas to check whether data frames are equal
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataFrameEqual)
        # example 
        # self.assertEqual(df, df_copy)  

        self.addTypeEqualityFunc(pd.Series, self.assertSeriesEqual)
        # example 
        # self.assertEqual(df['orig_timestamp'], df_copy['orig_timestamp']) 

    def test_gap_import(self):
        processed_dir='test_data'
        join_dir='test_data_output'
        config.combine_locations=True
        config.clean_spikes=True
        config.view_corrected_traces = True


        ############################
        # Read in a test file and check when it ends 
        call_csv_join_work_tapes(
        processed_dir=processed_dir,
        join_dir=join_dir,
        log_dir=processed_dir,
        wildcard_style='pse',
        year_start=2003,
        year_end=2003,
        day_start=181,
        day_end=181,
        stations=['S14'],
        logging_level=logging.DEBUG,
        )
        # 
        plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime(year=2003, julday=181),stations=['S14'],out_dir='test_data_output')
        
        stream_2003 = read('test_data_output/s14/2003/181/xa.s14..att.2003.181.0.mseed')
        stream_2003.merge()
        end_time_2003 = stream_2003[-1].stats.endtime
        print('End sample time for the comparison file ', end_time_2003)
        print('End timestamp for the comparison file ', UTCDateTime(stream_2003[-1].data[-1]))


        test_cutout_start=pd.to_datetime(datetime.datetime(2003, 6, 30, 0, 30),utc=True)
        test_cutout_end=pd.to_datetime(datetime.datetime(2, 7, 1, 0, 31 ),utc=True)
        
        call_csv_join_work_tapes(
        processed_dir=processed_dir,
        join_dir=join_dir,
        log_dir=processed_dir,
        wildcard_style='pse',
        year_start=1993,
        year_end=1993,
        day_start=181,
        day_end=181,
        # stations=['S12','S14','S15','S16'],
        # stations=['S12'],
        stations=['S14'],
        # stations=['S15'],
        # stations=['S16'],
        logging_level=logging.DEBUG,
        test=True, # these parameters are used for testing ONLY
        test_start=None, # these parameters are used for testing ONLY
        test_end=None, # these parameters are used for testing ONLY
        test_framejump=None, # these parameters are used for testing ONLY     
        test_cutout_start=test_cutout_start,
        test_cutout_end=test_cutout_end,
        )
        # # 
        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime(year=1993, julday=181),stations=['S14'],out_dir='test_data_output')
        # 
        # stream_1993 = read('test_data_output/s14/1993/181/xa.s14..att.1993.181.0.mseed')
        # stream_1993.merge()
        # end_time_1993 = stream_1993[-1].stats.endtime

        # print(UTCDateTime('1973-06-30T00:00:00.00000Z').julday)
        # exit()
        
        # '1973-06-30T00:00:00.00000Z'
        # 

        # # check the basic file (1973)
        # call_csv_join_work_tapes(
        # processed_dir=processed_dir,
        # join_dir=join_dir,
        # log_dir=processed_dir,
        # wildcard_style='pse',
        # year_start=1973,
        # year_end=1973,
        # day_start=181,
        # day_end=181,
        # # stations=['S12','S14','S15','S16'],
        # # stations=['S12'],
        # stations=['S14'],
        # # stations=['S15'],
        # # stations=['S16'],
        # logging_level=logging.DEBUG)
        # 
        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime(year=1973, julday=181),stations=['S14'],out_dir='test_data_output')

        stream_1973 = read('test_data_output/s14/1973/181/xa.s14..att.1973.181.0.mseed')
        stream_1973.merge()
        end_time_1973 = stream_1973[-1].stats.endtime

        # # small gap in the 1983 test data 
        # call_csv_join_work_tapes(
        # processed_dir=processed_dir,
        # join_dir=join_dir,
        # log_dir=processed_dir,
        # wildcard_style='pse',
        # year_start=1983,
        # year_end=1983,
        # day_start=181,
        # day_end=181,
        # # stations=['S12','S14','S15','S16'],
        # # stations=['S12'],
        # stations=['S14'],
        # # stations=['S15'],
        # # stations=['S16'],
        # logging_level=logging.DEBUG)
        # 
        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime(year=1983, julday=181),stations=['S14'],out_dir='test_data_output')

        stream_1983 = read('test_data_output/s14/1983/181/xa.s14..att.1983.181.0.mseed')
        stream_1983.merge()
        end_time_1983 = stream_1983[-1].stats.endtime

        self.assertEqual(end_time_1973.second, end_time_1983.second)
        self.assertEqual(end_time_1973.microsecond, end_time_1983.microsecond)

        



        # # big gap in the 1993 test data 
        # call_csv_join_work_tapes(
        # processed_dir=processed_dir,
        # join_dir=join_dir,
        # log_dir=processed_dir,
        # wildcard_style='pse',
        # year_start=1993,
        # year_end=1993,
        # day_start=181,
        # day_end=181,
        # # stations=['S12','S14','S15','S16'],
        # # stations=['S12'],
        # stations=['S14'],
        # # stations=['S15'],
        # # stations=['S16'],
        # logging_level=logging.DEBUG)
        # 
        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime(year=1993, julday=181),stations=['S14'],out_dir='test_data_output')
        # 
        # stream_1993 = read('test_data_output/s14/1993/181/xa.s14..att.1993.181.0.mseed')
        # stream_1993.merge()
        # end_time_1993 = stream_1993[-1].stats.endtime
        # 
        # self.assertEqual(end_time_1973.second, end_time_1993.second)
        # self.assertEqual(end_time_1973.microsecond, end_time_1993.microsecond)

        # mess up the frame count for a while in the 1993 test data 

        test_start=pd.to_datetime(datetime.datetime(1993, 6, 30, 20, 0),utc=True)
        test_end=pd.to_datetime(datetime.datetime(1993, 7, 1),utc=True)

        test_cutout_start=pd.to_datetime(datetime.datetime(1993, 6, 30, 20, 0),utc=True)
        test_cutout_end=pd.to_datetime(datetime.datetime(1993, 7, 1),utc=True)
        
        
        call_csv_join_work_tapes(
        processed_dir=processed_dir,
        join_dir=join_dir,
        log_dir=processed_dir,
        wildcard_style='pse',
        year_start=1993,
        year_end=1993,
        day_start=181,
        day_end=181,
        # stations=['S12','S14','S15','S16'],
        # stations=['S12'],
        stations=['S14'],
        # stations=['S15'],
        # stations=['S16'],
        logging_level=logging.DEBUG,
        test=True, # these parameters are used for testing ONLY
        test_start=test_start, # these parameters are used for testing ONLY
        test_end=test_end, # these parameters are used for testing ONLY
        test_framejump=40, # these parameters are used for testing ONLY     
        test_cutout_start=None,
        )
        # # 
        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime(year=1993, julday=181),stations=['S14'],out_dir='test_data_output')
        # 
        # stream_1993 = read('test_data_output/s14/1993/181/xa.s14..att.1993.181.0.mseed')
        # stream_1993.merge()
        # end_time_1993 = stream_1993[-1].stats.endtime

        # self.assertEqual(end_time_1973.second, end_time_1993.second)
        # self.assertEqual(end_time_1973.microsecond, end_time_1993.microsecond)



    test=False, # these parameters are used for testing ONLY
    test_start=None, # these parameters are used for testing ONLY
    test_end=None, # these parameters are used for testing ONLY
    test_framejump=40 # these parameters are used for testing ONLY


        # for file in glob.glob('/Users/cnunn/python_packages/pdart/test_suite/test_data/pse.a14*.1973.*.csv'):
        #     print(file)
        #     new_file_basename = file.replace('1973','1983')
        #     # print(new_file_basename)
        #     with open(os.path.join('test_data',new_file_basename), 'w') as out_file:
        #         new_lines = []
        #         with open(file, 'r') as f:
        #             lines = f.readlines()
        #             for i, line in enumerate(lines):
        #                 # print(line)
        #                 # print(type(line))
        #                 line = line.replace('1973','1983')
        #                 new_lines.append(line)
        #                 # print(line)
        #                 # if i > 0:
        #                 #     break
        #                 # 
        # 
        #     # 
        #         out_file.writelines(new_lines)
        # 
        # for file in glob.glob('/Users/cnunn/python_packages/pdart/test_suite/test_data/pse.a14*.1973.*.csv'):
        #     print(file)
        #     new_file_basename = file.replace('1973','1993')
        #     # print(new_file_basename)
        #     with open(os.path.join('test_data',new_file_basename), 'w') as out_file:
        #         new_lines = []
        #         with open(file, 'r') as f:
        #             lines = f.readlines()
        #             for i, line in enumerate(lines):
        #                 # print(line)
        #                 # print(type(line))
        #                 line = line.replace('1973','1993')
        #                 new_lines.append(line)
        #                 # print(line)
        #                 # if i > 0:
        #                 #     break
        #                 # 
        # 
        #     # 
        #         out_file.writelines(new_lines)


        # st = read('/Users/cnunn/python_packages/pdart/test_suite/test_data_output/s14/1973/181/xa.s14..mh1.1973.181.0.mseed')
        # print(st)
        # st.plot()


        
        # default_config() 
        # config.clean_spikes = True
        # 
        # df_gst = basic_timeseries_station()  
        # starttime0 = df_gst.corr_timestamp.iloc[0]
        # starttime0 = pd.Timestamp(starttime0)
        # 
        # ##############################
        # # first view a basic import
        # df_gst = basic_timeseries_station()  
        # df_gst['time_index'] = np.arange(len(df_gst))
        # 
        # stream = stream_import(df_gst,starttime0=starttime0,index0=0)
        # 
        # stream.plot(method='full',size=(1200,600))
        # 
        # ##############################
        # # add a single upward spike and don't clean it
        # df_gst = basic_timeseries_station()  
        # df_gst['time_index'] = np.arange(len(df_gst))
        # 
        # # add a spike 
        # spikes = 0
        # for row_i in range(25,26):
        #     df_gst.at[row_i,'orig_mh2_1'] = 800
        #     spikes += 1
        # 
        # config.clean_spikes = False
        # stream = stream_import(df_gst,starttime0=starttime0,index0=0)
        # 
        # for tr in stream:
        #     if tr.stats.channel == 'MH2' and tr.stats.station == 'S12':
        #         final = tr.data[row_i*4]
        # 
        # self.assertEqual(final, 800)
        
def suite():

    logging.basicConfig(filename='logs/clean_test.log', filemode='w', 
      level=logging.DEBUG,format='%(message)s')

    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(ImportTestCase, 'test'))
    
    return suite


if __name__ == '__main__':

    unittest.main(defaultTest='suite')
