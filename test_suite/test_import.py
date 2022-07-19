#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
Test Suite 
Old file for testing the imports. 

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
# from decimal import Decimal, getcontext TODO Remove 
from pandas.testing import assert_frame_equal
import logging, logging.handlers

import numpy as np

from obspy.core.utcdatetime import UTCDateTime
from pdart.csv_check_work_tapes import (
  all_drop_duplicates,
  all_split, station_fix_bad_timestamps, station_fix_missing_timestamps2, 
  initial_cleanup, detailed_report,
  calculate_gaps,
  station_fix_frames,
  add_or_minus_frame, minus_frame, frame_diff,
  station_simple_fix_timestamp,
  get_new_ground_station,all_drop_station_duplicates,all_station_duplicates_v2,
  all_flat_seismogram, check_compliance, calculate_delta4,
  calculate_cumsum, frame_diff_positive, calculate_frame_change, fix_timestamps
)
import pdart.config as config

# consecutive stations (the order is important)
STATIONS = ['S12', 'S15', 'S16', 'S14', 'S17']
# consecutive frames
ORIG_FRAMES = list(range(1,61))

FRAMES = list(range(0,90))

# nominal delta (sampling rate is 6.25 samples/s)
DELTA = 0.1509433962
INTERVAL = '603774us'

# default date format 
DATEFORMAT='%Y-%m-%dT%H:%M:%S.%fZ'

GZIP_FILENAME='test.csv.gz'

# from obspy import Stream, Trace, UTCDateTime, read, read_inventory
# from obspy.core.inventory import Channel, Inventory, Network, Station
# from obspy.core.compatibility import mock
# from obspy.core.stream import _is_pickle, _read_pickle, _write_pickle
# from obspy.core.util.attribdict import AttribDict
# from obspy.core.util.base import NamedTemporaryFile, _get_entry_points
# from obspy.core.util.obspy_types import ObsPyException
# from obspy.core.util.testing import streams_almost_equal
# from obspy.io.xseed import Parser

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

def basic_timeseries_station(no_records=6):
    df = basic_timeseries(no_records=no_records)

    # make a timeseries with only 'S12' data in it 
    df = initial_cleanup(df)

    df['delta4'] = 0.604    

    # detailed_report(df,'')
    grouped_st = df.groupby('orig_station')
    df_st = grouped_st.get_group('S12')
    df_st = df_st.copy(deep=True)

    # # make a new index
    df_st.reset_index(inplace=True,drop=True)



    return df_st

def basic_timeseries(no_records=3):

    random.seed(a=2)
    # orig_no,orig_timestamp,orig_frame,orig_station,clock_flag,orig_ground_station,orig_mh1_1,orig_mh2_1,orig_mhz_1,orig_mh1_2,orig_mh2_2,orig_mhz_2,orig_mh1_3,orig_mh2_3,orig_mhz_3,orig_mh1_4,orig_mh2_4,orig_mhz_4
    # 1,1976-03-01T00:01:44.696000Z,1,S12,0,9,531,491,481,531,491,481,531,491,481,531,491,481

    orig_timestamp = [UTCDateTime('1976-03-01T00:01:44.696000Z'), 
      UTCDateTime('1976-03-01T00:01:44.696000Z'),
      UTCDateTime('1976-03-01T00:01:44.378000Z'),
      UTCDateTime('1977-03-01T00:01:44.373000Z'),
      UTCDateTime('1976-03-01T00:01:44.671000Z')
    ]

    orig_mh1_1 = [531, 534, 534, 560, 530]
    orig_mh2_1 = [490, 491, 492, 493, 494]
    orig_mhz_1 = [490, 491, 492, 493, 494]

    orig_mh1_2 = [531, 534, 534, 560, 530]
    orig_mh2_2 = [490, 491, 492, 493, 494]
    orig_mhz_2 = [490, 491, 492, 493, 494]

    orig_mh1_3 = [531, 534, 534, 560, 530]
    orig_mh2_3 = [490, 491, 492, 493, 494]
    orig_mhz_3 = [490, 491, 492, 493, 494]

    orig_mh1_4 = [531, 534, 534, 560, 530]
    orig_mh2_4 = [490, 491, 492, 493, 494]
    orig_mhz_4 = [490, 491, 492, 493, 494]
    

# 1,1976-03-01T00:01:44.696000Z,1,S12,0,9,531,491,481,531,491,481,531,491,481,531,491,481
# 1,1976-03-01T00:01:44.638000Z,2,S15,0,9,504,505,491,504,505,491,504,505,491,504,505,491
# 1,1976-03-01T00:01:44.378000Z,3,,0,9,511,512,469,511,512,469,510,511,469,510,511,470
# 1,1976-03-01T00:01:44.373000Z,4,S14,0,9,511,490,509,511,490,509,511,490,508,511,490,508
# 1,1976-03-01T00:01:44.671000Z,5,S17,0,9,0,0,0,514,504,0,0,0,0,0,0,0

    start_no = 150

    for i in range(60*no_records-5):
        # add a random number onto the last random.choice (which has an index of 
        # -5 because there are 5 stations)
        orig_timestamp.append(orig_timestamp[-5] + random.choice([0.603, 0.604, 0.604, 0.604, 0.604]))
    orig_timestamp = [x.strftime(DATEFORMAT) for x in orig_timestamp]

    # logging.info(STATIONS)
    # logging.info(ORIG_FRAMES)

    # make the data frame from orig_timestamp
    df = pd.DataFrame(orig_timestamp, columns=['orig_timestamp'])
    df['orig_timestamp'] = pd.to_datetime(df['orig_timestamp'], format=DATEFORMAT)

    # tile the frames until the end of the record
    df['orig_frame'] = np.tile(ORIG_FRAMES, len(df)//len(ORIG_FRAMES) + 1)[:len(df)]
    # tile the stations until the end of the record
    df['orig_station'] = np.tile(STATIONS, len(df)//len(STATIONS) + 1)[:len(df)]    

    # repeat the origin number for a whole physical record 
    df['orig_no'] = list(np.repeat(np.arange(start_no,start_no+no_records), 60))
    # set the software clock flag to zero 
    df['clock_flag'] = 0 
    df['orig_ground_station'] = 9
    df['bit_synchronizer'] = '00011'
    df['sync'] = '1110001001000011101101'


    frame = [3,10,13,20,23]

    for i in range(60*no_records-5):
        # add a random number onto the last random.choice (which has an index of 
        # -5 because there are 5 stations)
        orig_mh1_1.append(orig_mh1_1[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mh2_1.append(orig_mh2_1[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mhz_1.append(orig_mhz_1[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))

        orig_mh1_2.append(orig_mh1_2[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mh2_2.append(orig_mh2_2[-5] + random.choice([0, 0, -1, 1, 2, -2, 0, 0]))
        orig_mhz_2.append(orig_mhz_2[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))

        orig_mh1_3.append(orig_mh1_3[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mh2_3.append(orig_mh2_3[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mhz_3.append(orig_mhz_3[-5] + random.choice([0, 0, 1, -1, -2, 2, 0, 0]))

        orig_mh1_4.append(orig_mh1_4[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mh2_4.append(orig_mh2_4[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mhz_4.append(orig_mhz_4[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))

        frame.append(add_or_minus_frame(frame[-5],1))

    df['orig_mh1_1'] = orig_mh1_1
    df['orig_mh2_1'] = orig_mh2_1
    df['orig_mhz_1'] = orig_mhz_1

    df['orig_mh1_2'] = orig_mh1_2
    df['orig_mh2_2'] = orig_mh2_2
    df['orig_mhz_2'] = orig_mhz_2

    df['orig_mh1_3'] = orig_mh1_3
    df['orig_mh2_3'] = orig_mh2_3
    df['orig_mhz_3'] = orig_mhz_3

    df['orig_mh1_4'] = orig_mh1_4
    df['orig_mh2_4'] = orig_mh2_4
    df['orig_mhz_4'] = orig_mhz_4

    df['frame'] = frame

    df['corr_ground_station'] = df['orig_ground_station']
    df['corr_timestamp'] = df['orig_timestamp']

    return df

def default_config():
    config.extra_ground_stations = []
    config.station_order = ['S12', 'S15', 'S16', 'S14', 'S17']
    config.cumsum_test = 10
    config.time_test = 5*60
    config.initial=False
    config.last_station=None
    config.last_timestamp=None

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

    def test_isolate(self):
        # isolated tests for chceck_compliance()

        #############################
        # Test calculate_delta4 (basic version)
        df_gst = basic_timeseries_station(40)
        
        df_gst = calculate_cumsum(df_gst)
        df_gst = calculate_delta4(df_gst)
        self.assertAlmostEqual(df_gst['delta4'].iloc[266], 0.603844, 5)
        self.assertAlmostEqual(df_gst['delta4'].iloc[267], 0.603777, 5)

        #############################
        # Test calculate_delta4 when using clock error calculation 
        df_gst = basic_timeseries_station(40)
        df_gst['clock_flag'] = 1
        
        df_gst = calculate_cumsum(df_gst)
        df_gst = calculate_delta4(df_gst)
        self.assertAlmostEqual(df_gst['delta4'].iloc[266], 0.603844, 5)
        self.assertAlmostEqual(df_gst['delta4'].iloc[267], 0.603777, 5)

        logging.info(df_gst.to_string())
        

        return 

        # #############################
        # Set up a -7777 error where the frame is reset from some number
        # to zero 
        
        df_gst = basic_timeseries_station(64)
        
        df_gst['shift_frame'] = df_gst['frame'].shift(-30)
        for row_i in range(147,720):
            df_gst.at[row_i,'frame'] = df_gst['shift_frame'].iloc[row_i]
            # df_gst.at[row_i,'clock_flag'] = 1
        
        # drop the end, because it's not valid, leaving 720 records 
        to_drop = list(range(720,768))
        df_gst.drop(to_drop,inplace=True)

        df_gst  = calculate_cumsum(df_gst)
        df_gst = calculate_delta4(df_gst)

        df_gst = calculate_sample_fit(df_gst,start_idx=None,end_idx=None)
        self.assertEqual(df_gst['actual_frame_gap'].iloc[148], -212.0)
        self.assertEqual(df_gst['actual_frame_gap'].iloc[146], -214.0)
        self.assertTrue(pd.isna(df_gst['sample_fit']).all())

        # XXXX

        return

        ##############################
        # Test double clock jump, without setting the clock flag 
        df_gst = basic_timeseries_station(60)
        
        before_len = len(df_gst)
        
        jumps = 0
        for row_i in range(300,360):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=1.206)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        self.assertEqual(df_gst['corr_gap_count'].iloc[300],-8888)
        self.assertEqual(df_gst['corr_gap_count'].iloc[360],-8888)
        
        df_gst = calculate_cumsum(df_gst)
        self.assertEqual(df_gst['cumsum'].iloc[299],300)
        self.assertEqual(df_gst['cumsum'].iloc[300],0)
        
        df_gst = calculate_delta4(df_gst)
        self.assertAlmostEqual(df_gst['delta4'].iloc[300], 0.603842, 5)
        
        df_gst = calculate_sample_fit(df_gst,start_idx=None,end_idx=None)
        self.assertTrue((df_gst['sample_fit'].iloc[300:359] == -2222).all())
        
        df_gst, rec_adjusted_timestamps = fix_timestamps(df_gst)
        self.assertEqual(rec_adjusted_timestamps,jumps)


        ##############################
        # Check calculate_sample_fit() without the whole record 
        df_gst = basic_timeseries_station(60)
        
        before_len = len(df_gst)
        
        jumps = 0
        for row_i in range(60,62):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=3.7)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        self.assertEqual(df_gst['corr_gap_count'].iloc[60],-8888)
        self.assertEqual(df_gst['corr_gap_count'].iloc[61],1)
        self.assertEqual(df_gst['corr_gap_count'].iloc[62],-8888)
        
        df_gst = calculate_cumsum(df_gst)
        self.assertEqual(df_gst['cumsum'].iloc[59],60)
        self.assertEqual(df_gst['cumsum'].iloc[60],0)
        
        df_gst = calculate_delta4(df_gst)
        self.assertAlmostEqual(df_gst['delta4'].iloc[25], 0.603844, 5)
        
        # test a small section with calculate_sample_fit()
        # before it's fixed, its -2222
        df_gst = calculate_sample_fit(df_gst,start_idx=60,end_idx=61)
        self.assertTrue((df_gst['sample_fit'].iloc[60:61] == -2222).all())
        

        df_gst, rec_adjusted_timestamps = fix_timestamps(df_gst)
        self.assertEqual(rec_adjusted_timestamps,jumps)
        
        # test a small section with calculate_sample_fit()
        # fixed in fix_timestamps, so recheck
        df_gst = calculate_sample_fit(df_gst,start_idx=60,end_idx=61)

        # the two affected records should be null
        self.assertTrue(pd.isna(df_gst['sample_fit'].iloc[60]))
        self.assertTrue(pd.isna(df_gst['sample_fit'].iloc[61]))

        # the next record should also be null
        self.assertTrue(pd.isna(df_gst['sample_fit'].iloc[61]))
        self.assertEqual(df_gst['actual_frame_gap'].iloc[60],-300.0)


        ##############################
        # Test a clock jump with 1.2 second jump that's a few records long
        df_gst = basic_timeseries_station(60)
        
        before_len = len(df_gst)
        
        jumps = 0
        
        for row_i in range(25,29):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=1.2)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1

        # broken from record 25 to 28 
        # so begin should be 24 because that is good? 
        # and end should be 29 because that is also good (it does have the -8888 
        # error but that is OK)

        
        for row_i in range(60,61):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=3.7)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1

        # begin should be 59 and end should be 61
        # broken at record 60
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        self.assertEqual(df_gst['corr_gap_count'].iloc[25],-8888)
        self.assertEqual(df_gst['corr_gap_count'].iloc[29],-8888)

        df_gst = calculate_cumsum(df_gst)
        self.assertEqual(df_gst['cumsum'].iloc[24],25)
        self.assertEqual(df_gst['cumsum'].iloc[25],0)

        df_gst = calculate_delta4(df_gst)
        self.assertAlmostEqual(df_gst['delta4'].iloc[25], 0.603844, 5)

        df_gst = calculate_sample_fit(df_gst,start_idx=None,end_idx=None)
        self.assertTrue((df_gst['sample_fit'].iloc[25:29] == -2222).all())
        self.assertTrue((df_gst['sample_fit'].iloc[60:61] == -2222).all())

        df_gst, rec_adjusted_timestamps = fix_timestamps(df_gst)
        self.assertEqual(rec_adjusted_timestamps,jumps)


    def test_compliance(self):
        # For testing check_compliance()
        # see also isolate(), for testing each of the individual methods
        # 
        # ##############################
        # # Test a clock jump with 1.2 second jump that's a few records long
        # # Test check_compliance()
        # df_gst = basic_timeseries_station()
        # # add a clock jump 
        # jumps = 0
        # for row_i in range(25,29):
        #     df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=1.2)
        #     df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
        #     df_gst.at[row_i,'clock_flag'] = 1
        #     jumps += 1
        # 
        # for row_i in range(60,61):
        #     df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=3.7)
        #     df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
        #     df_gst.at[row_i,'clock_flag'] = 1
        # single = 1
        # 
        # df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        # 
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # 
        # self.assertEqual(rec_adjusted_timestamps,jumps)
        # self.assertEqual(rec_fixed_simple_timestamp,single)
        # # valid_idx = 40, so actual_frame_gap = 60 - 40 = 20
        # self.assertEqual(df_gst['actual_frame_gap'].iloc[60], 20.0)
        # 
        # 
        # ##############################
        # # Test very small clock jump, without setting the clock flag 
        # df_gst = basic_timeseries_station(60)
        # 
        # jumps = 0
        # for row_i in range(300,360):
        #     df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=0.102)
        #     df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
        #     jumps += 1
        # 
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # 
        # # too small to trigger a change
        # self.assertTrue(pd.isna(df_gst['sample_fit']).all())
        # self.assertEqual(rec_fixed_frame,0)
        # self.assertEqual(rec_fixed_simple_timestamp,0)
        # self.assertEqual(rec_adjusted_timestamps,0)
        # self.assertEqual(rec_deleted_non_uniq,0)
        # 
        # ##############################
        # # Test bigger clock jump, without setting the clock flag 
        # df_gst = basic_timeseries_station(60)
        # 
        # before_len = len(df_gst)
        # 
        # jumps = 0
        # for row_i in range(300,360):
        #     df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=0.803)
        #     df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
        #     jumps += 1
        # 
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # # 
        # self.assertTrue(pd.isna(df_gst['sample_fit']).all())
        # self.assertEqual(rec_fixed_frame,0)
        # self.assertEqual(rec_fixed_simple_timestamp,0)
        # self.assertEqual(rec_adjusted_timestamps,jumps)
        # self.assertEqual(rec_deleted_non_uniq,0)
        # 
        # ##############################
        # # Test double clock jump, setting the clock flag 
        # df_gst = basic_timeseries_station(60)
        # 
        # before_len = len(df_gst)
        # 
        # jumps = 0
        # for row_i in range(300,360):
        #     df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=1.206)
        #     df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
        #     df_gst.at[row_i,'clock_flag'] = 1
        #     jumps += 1
        # 
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # 
        # after_len = len(df_gst)
        # 
        # # clock jump should be fixed 
        # self.assertTrue(pd.isna(df_gst['sample_fit']).all())
        # 
        # idx_list = (df_gst[(df_gst['corr_gap_count'] == -8888)]).index.tolist()
        # self.assertTrue(len(idx_list) == 0)
        # 
        # self.assertEqual(before_len,after_len)
        # 
        # # valid_idx = 461, so actual_frame_gap = 300 - 461 = -161
        # self.assertEqual(df_gst['actual_frame_gap'].iloc[300], -161.)
        # 
        # 
        # # return
        # 
        # ##############################
        # # Test really big jump
        # df_gst = basic_timeseries_station(60)
        # 
        # before_len = len(df_gst)
        # 
        # jumps = 0
        # for row_i in range(300,310):
        #     df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=1440000)
        #     df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
        #     df_gst.at[row_i,'clock_flag'] = 1
        #     jumps += 1
        # 
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # # 
        # idx_list = (df_gst[(df_gst['corr_gap'] < 0)]).index.tolist()
        # self.assertTrue(len(idx_list) == 0)
        # idx_list = (df_gst[(df_gst['sample_fit'] == -3333) | (df_gst['sample_fit'] == -2222)]).index.tolist()
        # self.assertTrue(len(idx_list) == 0)
        # 
        # idx_list = (df_gst[(df_gst['corr_gap_count'] == -8888)]).index.tolist()
        # self.assertTrue(len(idx_list) == 0)
        # 
        # after_len = len(df_gst)
        # self.assertEqual(before_len,after_len)
        # 
        # # return
        # # 
        # #############################
        # # Frame with the wrong number
        # df_gst = basic_timeseries_station(3)
        # 
        # df_gst.at[7,'frame'] = 27
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # 
        # self.assertEqual(rec_fixed_frame,1)
        # self.assertEqual(rec_fixed_simple_timestamp,0)
        # self.assertEqual(rec_adjusted_timestamps,0)
        # self.assertEqual(rec_deleted_non_uniq,0)
        # 
        # 
        # #############################
        # # Test repeated frame 
        # df_gst = basic_timeseries_station(60)
        # 
        # x1=df_gst.loc[[7],:]
        # 
        # df_gst = pd.concat([df_gst,x1]).sort_index()
        # df_gst = pd.concat([df_gst,x1]).sort_index()
        # df_gst.reset_index(inplace=True,drop=True)
        # repeated_frame = 2
        # before_len = len(df_gst)
        # 
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # after_len = len(df_gst)
        # 
        # self.assertEqual(before_len-repeated_frame,after_len)
        # self.assertEqual(rec_fixed_frame,0)
        # self.assertEqual(rec_fixed_simple_timestamp,0)
        # self.assertEqual(rec_adjusted_timestamps,0)
        # self.assertEqual(rec_deleted_non_uniq,repeated_frame)
        # 
        # # return
        # # 
        # # 
        # #############################
        # # Test clock jumps at the end 
        # df_gst = basic_timeseries_station(60)
        # 
        # # print(len(df_gst))
        # 
        # jumps = 0
        # for row_i in range(700,720):
        #     df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=.52)
        #     df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
        #     df_gst.at[row_i,'clock_flag'] = 1
        #     jumps += 1
        # 
        # print('This should raise SEVERE error, because it deletes records at the end.')
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # 
        # 
        # self.assertEqual(len(df_dropped),jumps)
        # self.assertEqual(len(df_gst),720-jumps)
        # self.assertEqual(before_len-repeated_frame,after_len)
        # 
        # self.assertEqual(rec_fixed_frame,0)
        # self.assertEqual(rec_fixed_simple_timestamp,0)
        # self.assertEqual(rec_adjusted_timestamps,0)
        # self.assertEqual(rec_deleted_non_uniq,0)
        # self.assertEqual(rec_deleted_timestamps,jumps)
        # 
        # #############################
        # # Test clock jumps at the beginning
        # df_gst = basic_timeseries_station(60)
        # 
        # print(len(df_gst))
        # 
        # jumps = 0
        # for row_i in range(0,20):
        #     df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=.52)
        #     df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
        #     df_gst.at[row_i,'clock_flag'] = 1
        #     jumps += 1
        # 
        # print('This should raise SEVERE error, because it deletes records at the beginning.')    
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # 
        # 
        # self.assertEqual(len(df_dropped),jumps)
        # self.assertEqual(len(df_gst),720-jumps)
        # 
        # self.assertEqual(rec_fixed_frame,0)
        # self.assertEqual(rec_fixed_simple_timestamp,0)
        # self.assertEqual(rec_adjusted_timestamps,0)
        # self.assertEqual(rec_deleted_non_uniq,0)
        # self.assertEqual(rec_deleted_timestamps,jumps)
        
        # return 
        # #############################
        # # Test clock jumps at the beginning
        # df_gst = basic_timeseries_station(3)
        # 
        # # Rubbish at the beginning of the record 
        # df_gst = basic_timeseries_station()
        # df_gst.at[0,'frame'] = 7
        # df_gst.at[1,'frame'] = 88
        # df_gst.at[2,'frame'] = 7
        # deleted_timestamps = 1
        # deleted_non_uniq = 2
        # total_len = len(df_gst) - deleted_timestamps - deleted_non_uniq
        # 
        # print('This should raise SEVERE error, because it deletes records at the beginning.')       
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # 
        # # 0 0 0 2 1
        # # print(len(df_dropped))
        # # logging.info(df_dropped.to_string())
        # 
        # self.assertEqual(rec_fixed_frame,0)
        # self.assertEqual(rec_fixed_simple_timestamp,0)
        # self.assertEqual(rec_adjusted_timestamps,0)
        # self.assertEqual(rec_deleted_non_uniq,deleted_non_uniq)
        # self.assertEqual(rec_deleted_timestamps,deleted_timestamps)
        # 
        # # 3 fewer records 
        # self.assertEqual(len(df_gst),total_len)
        # # 3 dropped records 
        # self.assertEqual(len(df_dropped),deleted_timestamps + deleted_non_uniq)


        # ##############################
        # # Test big jump thats real 
        # df_gst = basic_timeseries_station(600)
        # before_len = len(df_gst)
        # 
        # df_gst.drop(df_gst.index[list(range(100,7000))],inplace=True)
        # df_gst.reset_index(inplace=True,drop=True)
        # # left with only 300 records 
        # 
        # # adds the missing records back in 
        # df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        # print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        # rec_missing1 = df_gst['orig_no'].isna().sum()
        # # print(rec_missing1)
        # after_len = len(df_gst)
        # # print(len(df_gst))
        # 
        # self.assertEqual(rec_fixed_frame,0)
        # self.assertEqual(rec_fixed_simple_timestamp,0)
        # self.assertEqual(rec_adjusted_timestamps,0)
        # self.assertEqual(rec_deleted_non_uniq,0)
        # self.assertEqual(rec_deleted_timestamps,0)
        # self.assertEqual(before_len,after_len)
        # self.assertEqual(rec_missing1,7000-100)

        # #############################
        # Set up a -7777 error where the frame is reset from some number
        # to zero 
        df_gst = basic_timeseries_station(64)
        
        df_gst['shift_frame'] = df_gst['frame'].shift(-30)
        for row_i in range(147,720):
            df_gst.at[row_i,'frame'] = df_gst['shift_frame'].iloc[row_i]
            # df_gst.at[row_i,'clock_flag'] = 1
        
        # drop the end, because it's not valid, leaving 720 records 
        to_drop = list(range(720,768))
        df_gst.drop(to_drop,inplace=True)

        before_len = len(df_gst)
        
        # df_gst = calculate_frame_change(df_gst)
        # df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        after_len = len(df_gst)

        self.assertEqual(rec_fixed_frame,0)
        self.assertEqual(rec_fixed_simple_timestamp,0)
        self.assertEqual(rec_adjusted_timestamps,0)
        self.assertEqual(rec_deleted_non_uniq,0)
        self.assertEqual(rec_deleted_timestamps,0)
        self.assertEqual(before_len,after_len)    


        # #############################
        # Set up a two -7777 errors where the frame is reset from some number
        # to zero (the valid idx is 384, and the errors occur before)
        df_gst = basic_timeseries_station(64)
        
        # df_gst['shift_frame'] = df_gst['frame'].shift(-30)
        
        for row_i in range(147,768):
            
            df_gst.at[row_i,'frame'] = add_or_minus_frame(df_gst['frame'].iloc[row_i],30)
            # df_gst.at[row_i,'clock_flag'] = 1

        for row_i in range(207,768):
            df_gst.at[row_i,'frame'] = add_or_minus_frame(df_gst['frame'].iloc[row_i],30)


        # df_gst['shift_frame'] = df_gst['frame'].shift(-30)
        # for row_i in range(207,720):
        #     df_gst.at[row_i,'frame'] = df_gst['shift_frame'].iloc[row_i]
        
        # drop the end, because it's not valid, leaving 720 records 
        # to_drop = list(range(720,768))
        # df_gst.drop(to_drop,inplace=True)

        before_len = len(df_gst)

        # df_gst = calculate_frame_change(df_gst)
        # df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        # df_gst = calculate_cumsum(df_gst)
        # df_gst = calculate_delta4(df_gst)
        # df_gst = calculate_sample_fit(df_gst,start_idx=None,end_idx=None)



        # self.assertEqual(df_gst['actual_frame_gap'].iloc[148], -212.0)
        # self.assertEqual(df_gst['actual_frame_gap'].iloc[146], -214.0)
        # self.assertTrue(pd.isna(df_gst['sample_fit']).all())

        df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        after_len = len(df_gst)

        # logging.info(df_gst.to_string())
        
        self.assertEqual(rec_fixed_frame,0)
        self.assertEqual(rec_fixed_simple_timestamp,0)
        self.assertEqual(rec_adjusted_timestamps,0)
        self.assertEqual(rec_deleted_non_uniq,0)
        self.assertEqual(rec_deleted_timestamps,0)
        self.assertEqual(before_len,after_len)    


        # #############################
        # Set up a two -7777 errors where the frame is reset from some number
        # to zero (the valid idx is )
        df_gst = basic_timeseries_station(64)
        
        # df_gst['shift_frame'] = df_gst['frame'].shift(-30)
        
        for row_i in range(417,768):
            df_gst.at[row_i,'frame'] = add_or_minus_frame(df_gst['frame'].iloc[row_i],30)
            # df_gst.at[row_i,'clock_flag'] = 1

        for row_i in range(477,768):
            df_gst.at[row_i,'frame'] = add_or_minus_frame(df_gst['frame'].iloc[row_i],30)


        # df_gst['shift_frame'] = df_gst['frame'].shift(-30)
        # for row_i in range(207,720):
        #     df_gst.at[row_i,'frame'] = df_gst['shift_frame'].iloc[row_i]
        
        # drop the end, because it's not valid, leaving 720 records 
        # to_drop = list(range(720,768))
        # df_gst.drop(to_drop,inplace=True)


        # df_gst = calculate_frame_change(df_gst)
        # df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        before_len = len(df_gst)
        
        # df_gst = calculate_cumsum(df_gst)
        # df_gst = calculate_delta4(df_gst)
        # df_gst = calculate_sample_fit(df_gst,start_idx=None,end_idx=None)
        # logging.info(df_gst.to_string())
        # return



        # self.assertEqual(df_gst['actual_frame_gap'].iloc[148], -212.0)
        # self.assertEqual(df_gst['actual_frame_gap'].iloc[146], -214.0)
        # self.assertTrue(pd.isna(df_gst['sample_fit']).all())
        
        df_gst, df_dropped, rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps = check_compliance(df_gst,'test.gzip')
        print(rec_fixed_frame, rec_fixed_simple_timestamp, rec_adjusted_timestamps, rec_deleted_non_uniq, rec_deleted_timestamps)
        after_len = len(df_gst)
        
        
        self.assertEqual(rec_fixed_frame,0)
        self.assertEqual(rec_fixed_simple_timestamp,0)
        self.assertEqual(rec_adjusted_timestamps,0)
        self.assertEqual(rec_deleted_non_uniq,0)
        self.assertEqual(rec_deleted_timestamps,0)
        self.assertEqual(before_len,after_len)    


        # XXXX



        print('End of compliance')


    def test_station_fix_bad_timestamps(self):

        # ##############################
        # Test frame_diff_positive() 

        # really basic case 
        # 33 is the valid frame, 37 is the new one 
        single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(
            last_frame=33,new_frame=37,gap=2.415,delta4=0.603811,frame_correction=0)
        print('single_gap={}, actual_frame_gap={}, full_frames_est={}, time_error={}, percent_error={}'.format(single_gap, actual_frame_gap, full_frames_est, time_error, percent_error))
        self.assertEqual(actual_frame_gap,4)
        self.assertAlmostEqual(time_error,0.000244,6)
        self.assertAlmostEqual(percent_error,0.010104,6)
        self.assertAlmostEqual(single_gap,2.415/4)
        
        # reverse everything (-1 x a negative gap)
        single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(
            last_frame=29,new_frame=33,gap=-1*-2.415,delta4=0.603811,frame_correction=0)
        print('single_gap={}, actual_frame_gap={}, full_frames_est={}, time_error={}, percent_error={}'.format(single_gap, actual_frame_gap, full_frames_est, time_error, percent_error))
        self.assertEqual(actual_frame_gap,4)
        self.assertAlmostEqual(time_error,0.000244,6)
        self.assertAlmostEqual(percent_error,0.010104,6)
        self.assertAlmostEqual(single_gap,2.415/4)

        # # check some addition
        # new = add_or_minus_frame(89,1)
        # self.assertEqual(new, 0)
        # 
        # new = add_or_minus_frame(89,90)
        # self.assertEqual(new, 89)
        # 
        # new = add_or_minus_frame(89,89)
        # self.assertEqual(new, 88)
        # 
        # new = add_or_minus_frame(89,90*3+1)
        # self.assertEqual(new, 0)
        # 
        # new = add_or_minus_frame(89,90*3-1)
        # self.assertEqual(new, 88)
        # 
        # 
        # # check a zero 
        # new = add_or_minus_frame(89,0)
        # self.assertEqual(new, 89)
        # 
        # 
        # new = add_or_minus_frame(89,-1)
        # self.assertEqual(new, 88)
        # # 
        # new = add_or_minus_frame(89,-90)
        # self.assertEqual(new, 89)
        # # 
        # new = add_or_minus_frame(89,-90*3)
        # self.assertEqual(new, 89)
        # # # 
        # new = add_or_minus_frame(89,(-90*3)+1)
        # self.assertEqual(new, 0)
        # 
        # new = add_or_minus_frame(89,(-90*3)-1)
        # self.assertEqual(new, 88)
        # 
        # new = add_or_minus_frame(89,-89)
        # self.assertEqual(new, 0)
        # 
        # new = add_or_minus_frame(89,-91)
        # self.assertEqual(new, 88)


        single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(
            2, 33, -127.406*-1, 0.603855)
        # logging.info('single_gap={}, actual_frame_gap={}, full_frames_est={}, time_error={}, percent_error={}'.format(single_gap, actual_frame_gap, full_frames_est, time_error, percent_error))
        # print('single_gap={}, actual_frame_gap={}, full_frames_est={}, time_error={}, percent_error={}'.format(single_gap, actual_frame_gap, full_frames_est, time_error, percent_error))
        # .single_gap=0.6038199052132702, actual_frame_gap=211, full_frames_est=2, time_error=0.007405000000005657, percent_error=0.005812128157234084

        # single_gap=0.6038199052132702, actual_frame_gap=211, full_frames_est=2, time_error=0.007405000000005657, percent_error=0.005812128157234084
        self.assertEqual(actual_frame_gap,211)
        self.assertAlmostEqual(time_error,0.007405,6)
        self.assertAlmostEqual(percent_error,0.005812,6)
        self.assertAlmostEqual(single_gap,127.406/211,6)


        # without a correction, this will show the SAME actual_frame_gap,
        # but with a large error 
        single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(
            2, 33, (-127.406+0.6038*3)*-1, 0.603855)
        self.assertEqual(actual_frame_gap,211)
        self.assertAlmostEqual(time_error,0.007405+(0.6038*3),6)
        self.assertAlmostEqual(single_gap,(-127.406+0.6038*3)*-1/211,6)

        # case with a -7777 error
        # the gap will remain the same, but the new frame will be different 
        # positive gaps
        # valid frame is 33
        single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(
            last_frame=33,new_frame=36,gap=2.415,delta4=0.603811,frame_correction=1)
        print('single_gap={}, actual_frame_gap={}, full_frames_est={}, time_error={}, percent_error={}'.format(single_gap, actual_frame_gap, full_frames_est, time_error, percent_error))
        self.assertEqual(actual_frame_gap,4)
        self.assertAlmostEqual(time_error,0.000244,6)
        self.assertAlmostEqual(percent_error,0.010104,6)
        self.assertAlmostEqual(single_gap,2.415/4)



        # negative gaps
        # valid frame is 33
        # reverse everything (-1 x a negative gap)
        single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(
            last_frame=30,new_frame=33,gap=-1*-2.415,delta4=0.603811,frame_correction=1)
        print('single_gap={}, actual_frame_gap={}, full_frames_est={}, time_error={}, percent_error={}'.format(single_gap, actual_frame_gap, full_frames_est, time_error, percent_error))
        self.assertEqual(actual_frame_gap,4)
        self.assertAlmostEqual(time_error,0.000244,6)
        self.assertAlmostEqual(percent_error,0.010104,6)
        self.assertAlmostEqual(single_gap,2.415/4)

        # some real examples (this one is from before the error)
        # valid frame is 33
        # reverse everything (-1 x a negative gap)
        single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(
            last_frame=1,new_frame=33,gap=-1*-128.010,delta4=0.603855,frame_correction=0)
        print('single_gap={}, actual_frame_gap={}, full_frames_est={}, time_error={}, percent_error={}'.format(single_gap, actual_frame_gap, full_frames_est, time_error, percent_error))
        self.assertEqual(actual_frame_gap,212)
        # approximate time_error / percent_error, because we're not 
        # the same delta4
        self.assertTrue(time_error < 0.01)
        self.assertTrue(percent_error < 0.01)
        self.assertAlmostEqual(single_gap,128.010/212)


        # some real examples (this one is from after the error)
        # valid frame is 33
        # reverse everything (-1 x a negative gap)
        single_gap, actual_frame_gap, full_frames_est, time_error, percent_error, est_frame_gap = frame_diff_positive(
            last_frame=59,new_frame=33,gap=-1*-129.218,delta4=0.603855,frame_correction=-30)
        print('single_gap={}, actual_frame_gap={}, full_frames_est={}, time_error={}, percent_error={}'.format(single_gap, actual_frame_gap, full_frames_est, time_error, percent_error))
        self.assertEqual(actual_frame_gap,214)
        # approximate time_error / percent_error, because we're not 
        # the same delta4
        self.assertTrue(time_error < 0.01)
        self.assertTrue(percent_error < 0.01)
        self.assertAlmostEqual(single_gap,129.218/214)


    def test_consec_frame(self):
        default_config()
        
        ###################
        # Test the basic consecutive frame
        df_gst = basic_timeseries_station()
        frame1 = df_gst['frame'].iloc[0]

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        valid_frame = df_gst['frame'].iloc[-1]

        df_gst = station_consec_frames(df_gst,len(df_gst)-1,valid_frame)
        frame_tile1 = df_gst['frame_tile'].iloc[0]

        # logging.info('df_gst {}'.format(df_gst.to_string()))

        self.assertEqual(frame1, frame_tile1)

        ###################
        # Test frames with a tile zero error (backwards)
        df_gst = basic_timeseries_station(no_records=12)
        frame1 = df_gst['frame'].iloc[0]

        for i, row_i in enumerate(range(82,len(df_gst))):
            df_gst.at[row_i,'frame'] = i
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        # find the 'good' records
        df_gst['good'] = np.where((df_gst['corr_gap_count'] == 1 ), True, False)
        df_gst['good'] = np.where((df_gst['clock_flag'] == 0 ), df_gst['good'] , False)
        df_gst['cumsum'] = df_gst['good'] * (df_gst['good'].groupby((df_gst['good'] != df_gst['good'].shift()).cumsum()).cumcount() + 1)

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        # logging.info('df_gst after\n{}'.format(df_gst.to_string()))

        valid_frame = df_gst['frame'].iloc[-1]

        df_gst = station_consec_frames(df_gst,len(df_gst)-1,valid_frame)
        frame_tile1 = df_gst['frame_tile'].iloc[0]

        # logging.info('df_gst \n{}'.format(df_gst.to_string()))
        
        self.assertEqual(frame1, frame_tile1)


        ###################
        # Test frames with a tile zero error (forwards)
        df_gst = basic_timeseries_station(no_records=12)

        for i, row_i in enumerate(range(82,len(df_gst))):
            df_gst.at[row_i,'frame'] = i
        frame1 = df_gst['frame'].iloc[0]
        frame_last = df_gst['frame'].iloc[-1]
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        # find the 'good' records
        df_gst['good'] = np.where((df_gst['corr_gap_count'] == 1 ), True, False)
        df_gst['good'] = np.where((df_gst['clock_flag'] == 0 ), df_gst['good'] , False)
        df_gst['cumsum'] = df_gst['good'] * (df_gst['good'].groupby((df_gst['good'] != df_gst['good'].shift()).cumsum()).cumcount() + 1)

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        # logging.info('df_gst after\n{}'.format(df_gst.to_string()))

        valid_frame = df_gst['frame'].iloc[0]

        df_gst = station_consec_frames(df_gst,0,valid_frame)
        frame_tile1 = df_gst['frame_tile'].iloc[0]
        frame_tile_last = df_gst['frame_tile'].iloc[-1]

        self.assertEqual(frame1, frame_tile1)
        self.assertEqual(frame_last, frame_tile_last)

        ###################
        # Test frames with two tile zero errors (forwards)
        df_gst = basic_timeseries_station(no_records=18)

        for i, row_i in enumerate(range(82,len(df_gst))):
            df_gst.at[row_i,'frame'] = i
        for i, row_i in enumerate(range(168,len(df_gst))):
            df_gst.at[row_i,'frame'] = i
        frame1 = df_gst['frame'].iloc[0]
        frame_last = df_gst['frame'].iloc[-1]
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)



        # find the 'good' records
        df_gst['good'] = np.where((df_gst['corr_gap_count'] == 1 ), True, False)
        df_gst['good'] = np.where((df_gst['clock_flag'] == 0 ), df_gst['good'] , False)
        df_gst['cumsum'] = df_gst['good'] * (df_gst['good'].groupby((df_gst['good'] != df_gst['good'].shift()).cumsum()).cumcount() + 1)

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        # logging.info('df_gst after\n{}'.format(df_gst.to_string()))

        valid_frame = df_gst['frame'].iloc[0]

        df_gst = station_consec_frames(df_gst,0,valid_frame)
        frame_tile1 = df_gst['frame_tile'].iloc[0]
        frame_tile_last = df_gst['frame_tile'].iloc[-1]

        # logging.info('df_gst \n{}'.format(df_gst.to_string()))
        # logging.info('df_gst \n{}'.format(df_gst.to_string()))
        
        self.assertEqual(frame1, frame_tile1)
        self.assertEqual(frame_last, frame_tile_last)

        ###################
        # Test frames with two tile zero error (backwards)
        df_gst = basic_timeseries_station(no_records=18)


        for i, row_i in enumerate(range(82,len(df_gst))):
            df_gst.at[row_i,'frame'] = i
        for i, row_i in enumerate(range(168,len(df_gst))):
            df_gst.at[row_i,'frame'] = i
        frame1 = df_gst['frame'].iloc[0]     
        frame_last = df_gst['frame'].iloc[-1]   

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        # find the 'good' records
        df_gst['good'] = np.where((df_gst['corr_gap_count'] == 1 ), True, False)
        df_gst['good'] = np.where((df_gst['clock_flag'] == 0 ), df_gst['good'] , False)
        df_gst['cumsum'] = df_gst['good'] * (df_gst['good'].groupby((df_gst['good'] != df_gst['good'].shift()).cumsum()).cumcount() + 1)

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        # logging.info('df_gst after\n{}'.format(df_gst.to_string()))

        valid_frame = df_gst['frame'].iloc[-1]

        df_gst = station_consec_frames(df_gst,len(df_gst)-1,valid_frame)
        frame_tile1 = df_gst['frame_tile'].iloc[0]
        frame_tile_last = df_gst['frame_tile'].iloc[-1]

        logging.info('df_gst \n{}'.format(df_gst.to_string()))
        
        self.assertEqual(frame1, frame_tile1)
        self.assertEqual(frame_last, frame_tile_last)

    
    # Tested
    def test_frame_stuff(self):
        default_config()

        #############################
        # repeated frame 
        df_gst = basic_timeseries_station()
        
        x1=df_gst.loc[[7],:]
        
        # x1.A+=2
        # x1.B+=7
        df_gst = pd.concat([df_gst,x1]).sort_index()
        df_gst = pd.concat([df_gst,x1]).sort_index()
        # reindex 
        df_gst.reset_index(inplace=True,drop=True)
        repeated_frame = 2

        # add a clock jump 
        jumps = 0
        
        for row_i in range(8,9):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=-0.001)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  

        logging.info('Before xx \n{}'.format(df_gst.to_string()))
        
        df_gst, rec_fixed_frame, rec_repeated_frame = station_fix_frames(df_gst)
        
        self.assertEqual(rec_fixed_frame, 0)
        self.assertEqual(rec_repeated_frame, repeated_frame)

    def test_breaks(self):
        default_config()

        ###################
        # Test what a real break looks like 
        df_gst = basic_timeseries_station(20)

        # delete lots of records records 
        lst = list(range(25,125))
        df_gst.drop(lst,inplace=True)
        # reindex 
        df_gst.reset_index(inplace=True,drop=True)
        missing_total = lst[-1] - lst[0] + 2
        jumps = missing_total - 1

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  

        df_gst = station_fix_missing_timestamps2(df_gst,GZIP_FILENAME)

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 

        rec_missing_total = df_gst['orig_no'].isna().sum()

        self.assertEqual(rec_missing_total, jumps)

        frame_change = frame_diff(last_frame=27,new_frame=38,gap=60.988)
        self.assertEqual(frame_change, missing_total)

        new = add_or_minus_frame(4,2)
        self.assertEqual(new, 6)

        new = add_or_minus_frame(4,90)
        self.assertEqual(new, 4)

        new = add_or_minus_frame(4,360)
        self.assertEqual(new, 4)

        new = add_or_minus_frame(89,1)
        self.assertEqual(new, 0)

        new = add_or_minus_frame(89,2)
        self.assertEqual(new, 1)

        new = add_or_minus_frame(30,70)
        self.assertEqual(new, 10)

        diff = minus_frame(old_frame=30,new_frame=10)
        self.assertEqual(diff, 70)

        diff = minus_frame(old_frame=10,new_frame=30)
        self.assertEqual(diff, 20)

        frame_change = frame_diff(last_frame=10,new_frame=11,gap=1*0.6038)  
        self.assertEqual(frame_change, 1)

        frame_change = frame_diff(last_frame=10,new_frame=30,gap=20*0.6038)
        self.assertEqual(frame_change, 20)

        frame_change = frame_diff(last_frame=30,new_frame=10,gap=70*0.6038)
        self.assertEqual(frame_change, 70)
        
        print('End of break')

    def test_frame_change(self):
        default_config()

        ##############################
        # Test a very long frame change

        frame_no = (1*90+1)
        frame_change = frame_diff(last_frame=0,new_frame=1,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)
        
        frame_no = (1*90+0)
        frame_change = frame_diff(last_frame=0,new_frame=0,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)
        
        frame_no = (20*90+1)
        frame_change = frame_diff(last_frame=0,new_frame=1,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)
        
        frame_no = (20*90+0)
        frame_change = frame_diff(last_frame=0,new_frame=0,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)
        
        frame_no = (20*90-1)
        frame_change = frame_diff(last_frame=0,new_frame=89,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)
        
        frame_no = (200*90+0)
        frame_change = frame_diff(last_frame=0,new_frame=0,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)
        
        frame_no = (2000*90+0)
        frame_change = frame_diff(last_frame=0,new_frame=0,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)

        frame_no = (2000*90+3)
        frame_change = frame_diff(last_frame=0,new_frame=3,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)

        frame_no = (20*90+3)
        frame_change = frame_diff(last_frame=0,new_frame=3,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)

        frame_no = (2000*90+3)
        frame_change = frame_diff(last_frame=0,new_frame=3,gap=frame_no*0.6038)
        self.assertEqual(frame_change, frame_no)

        # the gap needs to be very accurate to work now 
        frame_no = (20000*90+3)
        print('gap=',frame_no*DELTA*4, frame_no*DELTA*4/3600 )
        frame_change = frame_diff(last_frame=0,new_frame=3,gap=frame_no*DELTA*4)
        self.assertEqual(frame_change, frame_no)

    def test_ground_station(self):
        default_config()
        config.extra_ground_stations = [600,700,701]
        new_ground_station = get_new_ground_station(ground_station1=6)
        self.assertEqual(new_ground_station, 601)

    def test_damaged_timestamps(self):
        default_config()

        ##############################
        # One end clock jump, so the data are removed at the end
        df_gst = basic_timeseries_station()
        # add a clock jump 
        jumps = 0
        
        for row_i in range(69,72):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst.reset_index(inplace=True,drop=True)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  
        
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        logging.info('After\n{}'.format(df_gst.to_string()))  
        self.assertEqual(rec_adjusted_timestamps, 0)
        self.assertEqual(rec_deleted_timestamps, jumps)   
        
        
        
        ##############################
        # One end clock jump, penultimate record
        df_gst = basic_timeseries_station()
        # add a clock jump 
        jumps = 0
        
        for row_i in range(70,72):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst.reset_index(inplace=True,drop=True)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 
        
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        
        self.assertEqual(rec_adjusted_timestamps, 0)
        self.assertEqual(rec_deleted_timestamps, jumps)  
        
        
        
        ##############################
        # One end clock jump, last record
        df_gst = basic_timeseries_station()
        # add a clock jump 
        jumps = 0
        
        for row_i in range(71,72):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst.reset_index(inplace=True,drop=True)
        
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 
        
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME) 
        
        self.assertEqual(rec_adjusted_timestamps, 0)
        self.assertEqual(rec_deleted_timestamps, jumps)   



        ##############################
        # Begining clock jump, so the data are updated at the beginning
        df_gst = basic_timeseries_station()
        # add a clock jump 
        jumps = 0
        
        for row_i in range(0,2):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst.reset_index(inplace=True,drop=True)

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        # logging.info('Before\n{}'.format(df_gst.to_string()))    
        
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)

        self.assertEqual(rec_adjusted_timestamps, 0)
        self.assertEqual(rec_deleted_timestamps, jumps)    

        print('End of damaged_timestamps')


    def test_damaged_frames(self):
        default_config()

        ##############################
        # one damaged frame
        df_gst = basic_timeseries_station()
        
        df_gst.at[8,'frame'] = 7
        fixed_frame = 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  
        
        df_gst, rec_fixed_frame, rec_repeated_frame = station_fix_frames(df_gst)
        
        self.assertEqual(rec_fixed_frame, fixed_frame)
        self.assertEqual(rec_repeated_frame, 0)
        
        ##############################
        # one damaged frame
        df_gst = basic_timeseries_station()
        
        df_gst.at[8,'frame'] = 7
        df_gst.at[9,'frame'] = 14
        damaged_frame = 2
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  
        # logging.info('Before\n{}'.format(df_gst.to_string()))
        df_gst, rec_fixed_frame, rec_repeated_frame = station_fix_frames(df_gst)
        print(rec_fixed_frame, rec_repeated_frame)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 
        # logging.info('Middle\n{}'.format(df_gst.to_string()))

        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        # logging.info('After\n{}'.format(df_gst.to_string()))
        
        # 2 are damaged, and 3 are removed in station_fix_frames
        self.assertEqual(rec_fixed_frame, 0)
        self.assertEqual(rec_repeated_frame, damaged_frame+1)
        self.assertEqual(rec_adjusted_timestamps, 0)
        self.assertEqual(rec_deleted_timestamps, 0)
        
        
        ##############################
        # one frame equal to previous frame
        df_gst = basic_timeseries_station()
        
        df_gst.at[5,'frame'] = 7
        fixed_frame = 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)      
        
        df_gst, rec_fixed_frame, rec_repeated_frame = station_fix_frames(df_gst)
        
        self.assertEqual(rec_repeated_frame, 0)
        self.assertEqual(rec_fixed_frame, fixed_frame)
        self.assertEqual(rec_repeated_frame, 0)
        
        ##############################
        # repeated frame 
        df_gst = basic_timeseries_station()
        
        x1=df_gst.loc[[7],:]
        # x1.A+=2
        # x1.B+=7
        df_gst = pd.concat([df_gst,x1]).sort_index()
        df_gst = pd.concat([df_gst,x1]).sort_index()
        # reindex 
        df_gst.reset_index(inplace=True,drop=True)
        repeated_frame = 2
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  
        
        df_gst, rec_fixed_frame, rec_repeated_frame = station_fix_frames(df_gst)
        
        self.assertEqual(rec_fixed_frame, 0)
        self.assertEqual(rec_repeated_frame, repeated_frame)
        
        ##############################    
        # 90 missing records, checking that it doesn't flag as an error 
        df_gst = basic_timeseries_station(20)
        
        # delete 90 records 
        lst = list(range(25,114))
        df_gst.drop(lst,inplace=True)
        # reindex 
        df_gst.reset_index(inplace=True,drop=True)
        missing_total = lst[-1] - lst[0] + 2
        jumps = missing_total - 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  
        
        # there is no frame error 
        df_gst, rec_fixed_frame, rec_repeated_frame = station_fix_frames(df_gst)
        
        self.assertEqual(rec_fixed_frame, 0)
        self.assertEqual(rec_repeated_frame, 0)
        
        
        ##############################    
        # frames equal to previous frame
        df_gst = basic_timeseries_station(12)
        
        fixed_frame = 0
        df_gst.at[5,'frame'] = df_gst.at[4,'frame']
        fixed_frame += 1
        df_gst.at[70,'frame'] = df_gst['frame'].iloc[71]
        fixed_frame += 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  
        
        # there is no frame error 
        df_gst, rec_fixed_frame, rec_repeated_frame = station_fix_frames(df_gst)
        
        self.assertEqual(rec_fixed_frame, fixed_frame)
        self.assertEqual(rec_repeated_frame, 0)
        print('End of damaged_frames()')

    def test_missing(self):
        default_config()

        ##############################
        # test inserting blank records into the gaps
        df_gst = basic_timeseries_station()

        # delete some records 
        df_gst.drop([7,8,12],inplace=True)
        # reindex 
        df_gst.reset_index(inplace=True,drop=True)
        missing_total = 3

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  

        df_gst = station_fix_missing_timestamps2(df_gst,GZIP_FILENAME)

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 

        rec_missing_total = df_gst['orig_no'].isna().sum()

        self.assertEqual(rec_missing_total, missing_total)

        logging.info('after\n{}'.format(df_gst.to_string()))


    


    def test_begin_end(self):
        default_config()

        ####################
        # first record has clock jump (delete first 2 records)
        df_gst = basic_timeseries_station(20)
        
        jumps = 0
        # clock jump 
        for row_i in range(0,1):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        # now we want it to delete the first two records
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        self.assertEqual(rec_deleted_timestamps, jumps)
        
        ####################
        # second record has clock jump (adjust the second and third records)
        df_gst = basic_timeseries_station(20)
        
        jumps = 0
        # clock jump 
        for row_i in range(1,2):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        # now we want it to adjust the second and third records
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        # logging.info('After\n{}'.format(df_gst.to_string()))
        
        self.assertEqual(rec_adjusted_timestamps, jumps)
        
        ####################
        # third record has clock jump (adjust the third and fourth records)
        df_gst = basic_timeseries_station(20)
        
        jumps = 0
        # clock jump 
        for row_i in range(3,4):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)

        logging.info('Before\n{}'.format(df_gst.to_string()))
        
        # now we want it to adjust the second and third records
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        logging.info('After\n{}'.format(df_gst.to_string()))
        
        # because the cummulative sum isn't high enough, it uses the
        # first record - this means that there is an extra record lost
        self.assertEqual(rec_adjusted_timestamps, jumps+1)
        # Checked

        ####################
        # third record has clock jump (adjust the third and fourth records)
        df_gst = basic_timeseries_station(20)
        
        # print(len(df_gst))
        
        jumps = 0
        # clock jump 
        for row_i in range(3,4):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        # now we want it to adjust the second and third records
        df_gst, rec_fixed_simple_timestamp = station_simple_fix_timestamp(df_gst)
        
        logging.info('After\n{}'.format(df_gst.to_string()))
        
        # this will fix just one
        self.assertEqual(rec_fixed_simple_timestamp, jumps)


        ####################
        # test station_simple_fix_timestamp with a small jump 
        df_gst = basic_timeseries_station(20)
        
        # print(len(df_gst))
        
        jumps = 0
        # clock jump 
        for row_i in range(210,211):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=0.003)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        # delete a record before this error 
        to_drop = [200]
        df_gst.drop(to_drop,inplace=True)
        # jumps -= len(to_drop)
        
        df_gst.reset_index(inplace=True,drop=True)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        logging.info('Before\n{}'.format(df_gst.to_string()))
        
        # # now we want it to adjust the second and third records
        # df_gst, rec_fixed_simple_timestamp = station_simple_fix_timestamp(df_gst)
        # 
        
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        
        print(rec_adjusted_timestamps, rec_deleted_timestamps)
        
        logging.info('After\n{}'.format(df_gst.to_string()))
        
        # this will fix just one
        # self.assertEqual(rec_fixed_simple_timestamp, jumps)


        ####################
        # test station_simple_fix_timestamp with a small jump 
        df_gst = basic_timeseries_station(20)

        
        print(len(df_gst))
        
        jumps = 0
        # clock jump 
        # error in first section
        for row_i in range(10,11):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=0.003)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1

        # error in second section
        for row_i in range(150,151):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=0.003)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        # delete a record in the middle
        to_drop = [140]
        df_gst.drop(to_drop,inplace=True)
        # jumps -= len(to_drop)
        
        df_gst.reset_index(inplace=True,drop=True)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        logging.info('Before\n{}'.format(df_gst.to_string()))
        
        # # now we want it to adjust the second and third records
        # df_gst, rec_fixed_simple_timestamp = station_simple_fix_timestamp(df_gst)
        # 
        
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        
        print(rec_adjusted_timestamps, rec_deleted_timestamps)
        
        logging.info('After\n{}'.format(df_gst.to_string()))
        
        # this will fix just one
        self.assertEqual(len(df_gst), 240)
        


        ####################
        # really big one
        df_gst = basic_timeseries_station(20)

        
        print(len(df_gst))
        
        jumps = 0
        # clock jump 
        # error in first section
        for row_i in range(10,11):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=0.003)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1

        # error in second section
        for row_i in range(150,151):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=0.003)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1

        # error in second section
        for row_i in range(160,161):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=0.003)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        for row_i in range(161,152):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.003)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        for row_i in range(162,163):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=1.003)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        for row_i in range(163,164):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=0.003)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        # delete a record in the middle
        to_drop = [140]
        df_gst.drop(to_drop,inplace=True)
        # jumps -= len(to_drop)
        
        df_gst.reset_index(inplace=True,drop=True)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)
        
        logging.info('Before\n{}'.format(df_gst.to_string()))
        
        # # now we want it to adjust the second and third records
        # df_gst, rec_fixed_simple_timestamp = station_simple_fix_timestamp(df_gst)
        # 

        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        
        print(rec_adjusted_timestamps, rec_deleted_timestamps)
        
        logging.info('After\n{}'.format(df_gst.to_string()))
        
        # this will fix just one
        self.assertEqual(len(df_gst), 240)

        print('End of begin_end x')


    def test_nasty(self):
        default_config()

        ####################
        # a frame jump with a lot of missing records in the middle 
        # station_fix_frames
        # station_fix_timeskips
        # station_fix_missing_timestamps

        df_gst = basic_timeseries_station(20)

        jumps = 0
        # test a clock jump with a couple of missing values 
        # add a clock jump 

        for row_i in range(25,45):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1
        
        for row_i in range(65,80):
            df_gst.at[row_i,'orig_timestamp'] = df_gst['orig_timestamp'].iloc[row_i] + pd.Timedelta(seconds=4.02)
            df_gst.at[row_i,'corr_timestamp'] = df_gst['orig_timestamp'].iloc[row_i]
            df_gst.at[row_i,'clock_flag'] = 1
            jumps += 1

        repeated_frame = 0
        # repeated frame 
        x1=df_gst.loc[[160],:]
        df_gst = pd.concat([df_gst,x1]).sort_index()
        # reindex 
        # 158, 159
        df_gst.reset_index(inplace=True,drop=True)
        repeated_frame += 1

        x1=df_gst.loc[[31],:]
        df_gst = pd.concat([df_gst,x1]).sort_index()
        df_gst = pd.concat([df_gst,x1]).sort_index()
        # reindex 
        df_gst.reset_index(inplace=True,drop=True)
        repeated_frame += 2

        fixed_frame = 0
        # frames equal to previous frame
        df_gst.at[5,'frame'] = df_gst.at[4,'frame']
        fixed_frame += 1
        df_gst.at[70,'frame'] = df_gst['frame'].iloc[71]
        fixed_frame += 1


        # delete lots of records in the middle
        lst = []
        lst = [6, 9]
        lst.extend(list(range(36,39)))
        lst.extend([237])
        df_gst.drop(lst,inplace=True)
        # reindex 
        df_gst.reset_index(inplace=True,drop=True)
        missing_total = len(lst)
        actualjumps = jumps - missing_total

        # since record 6 is deleted, record 5 cannot be fixed 
        missing_total += 1
        fixed_frame -= 1
        # record 5 is both a repeated frame and an out of sync frame, 
        # but it can't be fixed because of the missing record near it
        repeated_frame += 1

        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 

        # df_gst, df_dropped, rec_damaged_outside_tolerance, rec_fixed_frame, rec_repeated_frame, rec_adjusted_timestamps, rec_deleted_timestamps, rec_fixed_simple_timestamp = station_fix_timestamps(df_gst,GZIP_FILENAME)   

        # self.assertEqual(rec_fixed_frame, fixed_frame)
        # self.assertEqual(rec_repeated_frame, repeated_frame)    
        # rec_missing_total = df_gst['orig_no'].isna().sum()
        # self.assertEqual(rec_missing_total, missing_total)

        # THIS ISN"T REALLY TESTING ANYTHING NOW!

        print('End of nasty')

    def test_read_in(self):
        default_config()


        # working with some nasty data 


        df = pd.read_csv('/Users/cnunn/lunar_data/PDART_CSV_WORK_TAPES_FULL/wtn.13.15.csv.gz', dtype=str)

        df = initial_cleanup(df)
        
        start_i = (df[(df['orig_timestamp'] == '1977-02-27T12:49:06.698000Z') & 
          (df['orig_station'] == 'S12') & 
          (df['corr_ground_station'] == 7)].index[0])
        end_i = (df[(df['orig_timestamp'] == '1977-02-27T12:49:09.221000Z') & 
          (df['orig_station'] == 'S14') & 
          (df['corr_ground_station'] == 7)].index[0])
        
        # 86,"1977-02-27T12:49:06.698000Z",1,"S12",1,7,182,344,117,802,740,913,11,604,861,998,445,532,"10001",31,"0101101111111001110010"
        # 89,"1977-02-27T12:49:09.221000Z",60,"S14",0,7,508,519,449,511,520,449,514,521,449,518,522,450,"00011",80,"1110001001000011101101"
        
        # # these 2 are examples 
        # 307754       88 1977-02-27 12:48:54.128000+00:00          15          S16           0                    7         495         519         448         494         519         448         494         518         448         493         518         448            00011     55  1110001001000011101101 1977-02-27 12:48:54.128000+00:00       NaN             1               1                    7       0  False  False
        # 307755       88 1977-02-27 12:48:54.128000+00:00          16          S14           0      
        
        df = df.iloc[start_i:end_i+1].copy()
        df.reset_index(inplace=True)
        logging.info('head \n{}'.format(df.to_string()))
        logging.info(len(df))
        
        df, rec_station_duplicates = all_drop_station_duplicates(df,'')
    
        
    def test_row_duplicatons(self):
        default_config()

        ####################
        # 4 duplications 1 row apart 
        df = basic_timeseries(no_records=100)
        cols = list(df.columns) 
        # copy 4 rows from another station
        copied = 0
        for row_i in range(72,92,5):
            current_station = df['orig_station'].iloc[row_i]
            df.at[row_i,cols] = df[cols].iloc[row_i-1]
            df.at[row_i,'orig_station'] = current_station
            copied += 1
        
        # copy 1 row from the SAME station (not an error)
        for row_i in range(105,106):
            current_timestamp = df['orig_timestamp'].iloc[row_i]
            current_frame = df['frame'].iloc[row_i]
            df.at[row_i,cols] = df[cols].iloc[row_i-5]
            df.at[row_i,'orig_timestamp'] = current_timestamp
            df.at[row_i,'corr_timestamp'] = current_timestamp
            df.at[row_i,'frame'] = current_frame
        
        df, rec_station_duplicates = all_station_duplicates_v2(df,'test_import.py')
        self.assertEqual(rec_station_duplicates, copied) 
        
        
        ####################
        # 4 duplications 2 rows apart 
        df = basic_timeseries(no_records=100)
        cols = list(df.columns) 
        # copy 4 rows from another station
        copied = 0
        for row_i in range(72,92,5):
            current_station = df['orig_station'].iloc[row_i]
            df.at[row_i,cols] = df[cols].iloc[row_i-2]
            df.at[row_i,'orig_station'] = current_station
            copied += 1
        
        df, rec_station_duplicates = all_station_duplicates_v2(df,'test_import.py')
        self.assertEqual(rec_station_duplicates, copied) 
        
        ####################
        # 4 duplications 3 rows apart 
        df = basic_timeseries(no_records=100)
        cols = list(df.columns) 
        # copy 4 rows from another station
        copied = 0
        for row_i in range(72,92,5):
            current_station = df['orig_station'].iloc[row_i]
            df.at[row_i,cols] = df[cols].iloc[row_i-3]
            df.at[row_i,'orig_station'] = current_station
            copied += 1
        
        df, rec_station_duplicates = all_station_duplicates_v2(df,'test_import.py')
        self.assertEqual(rec_station_duplicates, copied) 
        
        ####################
        # 4 duplications 4 rows apart 
        df = basic_timeseries(no_records=100)
        cols = list(df.columns) 
        # copy 4 rows from another station
        copied = 0
        for row_i in range(72,92,5):
            current_station = df['orig_station'].iloc[row_i]
            df.at[row_i,cols] = df[cols].iloc[row_i-4]
            df.at[row_i,'orig_station'] = current_station
            copied += 1
        
        df, rec_station_duplicates = all_station_duplicates_v2(df,'test_import.py')
        self.assertEqual(rec_station_duplicates, copied) 
        
        ####################
        # 10 duplications 1 rows apart, only 2 stations
        df = basic_timeseries(no_records=100)
        df = df[(df.orig_station == 'S12') | (df.orig_station == 'S14')]
        df.reset_index(inplace=True)
        cols = list(df.columns) 
        # copy 4 rows from another station
        copied = 0
        for row_i in range(72,92,2):
            current_station = df['orig_station'].iloc[row_i]
            df.at[row_i,cols] = df[cols].iloc[row_i-1]
            df.at[row_i,'orig_station'] = current_station
            copied += 1
        
        config.station_order = ['S12','S14']
        df, rec_station_duplicates = all_station_duplicates_v2(df,'test_import.py')
        self.assertEqual(rec_station_duplicates, copied) 
        # reset the station order 
        config.station_order = ['S12', 'S15', 'S16', 'S14', 'S17']


        ####################
        # row duplications for the same station
        df = basic_timeseries(no_records=3)
        cols = list(df.columns) 
        # copy 4 rows from the same station
        copied = 0
        for i, row_i in enumerate(range(72,92,5)):
            print(row_i)
            current_frame = df['frame'].iloc[row_i]
            df.at[row_i,cols] = df[cols].iloc[row_i-5*(i+1)]
            df.at[row_i,'frame'] = current_frame
            copied += 1

        # logging.info('head \n{}'.format(df.to_string()))
        all_flat_seismogram(df,'test_import.py')
        # self.assertEqual(rec_station_duplicates, copied) 

        print('End of row_duplications')


    # def test_duplications_in_stations(self):
    # old version
    #     df = basic_timeseries(no_records=100)
    #     cols = list(df.columns) 
    #     # copy 4 rows from another station
    #     copied = 0
    #     for row_i in range(72,92,5):
    #         current_station = df['orig_station'].iloc[row_i]
    #         df.at[row_i,cols] = df[cols].iloc[row_i-1]
    #         df.at[row_i,'orig_station'] = current_station
    #         copied += 1
    # 
    #     # copy 1 row from the SAME station (not an error)
    #     for row_i in range(105,106):
    #         current_timestamp = df['orig_timestamp'].iloc[row_i]
    #         current_frame = df['frame'].iloc[row_i]
    #         df.at[row_i,cols] = df[cols].iloc[row_i-5]
    #         df.at[row_i,'orig_timestamp'] = current_timestamp
    #         df.at[row_i,'corr_timestamp'] = current_timestamp
    #         df.at[row_i,'frame'] = current_frame
    # 
    #     logging.info('begin \n{}'.format(df.iloc[58:106].to_string()))
    #     df, rec_station_duplicates = all_drop_station_duplicates(df,'')
    #     logging.info('end \n{}'.format(df.iloc[58:106].to_string()))
    #     self.assertEqual(rec_station_duplicates, copied)    


        # 
        # df['found'] = 0
        # from datetime import datetime, timedelta



        # for a in range(0,10):
        #     # df['Draws'] = soc_iter('Arsenal', df['HomeTeam'].values, df['AwayTeam'].values,df['FTR'.values])
        # df['found'] = search_copied_station(df['corr_timestamp'].values)



        # (180000, 4)
        # print(array.shape)
        # # very basic
        # # t4 : 0.001164
        # print(len(df))

        # old way 
        # t1 = datetime.now()
        # for i in range(0,len(df) -1):
        #     corr_timestamp1 = df['corr_timestamp'].iloc[i]
        #     corr_timestamp_low = corr_timestamp1 - pd.Timedelta(seconds=5)
        #     corr_timestamp_high = corr_timestamp1 + pd.Timedelta(seconds=5)
        #     frame1 = df['frame'].iloc[i]
        #     station1 = df['orig_station'].iloc[i]
        # t2 = datetime.now()
        # logging.info('old : {}'.format((t2-t1).total_seconds()))
        # print('old ', (t2-t1).total_seconds())
        # old  2.536819

        # new way

        # df['corr_timestamp_low'] = df['corr_timestamp'] - pd.Timedelta(seconds=5)
        # df['corr_timestamp_high'] = df['corr_timestamp'] + pd.Timedelta(seconds=5)
        # array = df[['corr_timestamp_low', 'corr_timestamp_high', 'frame', 'orig_station','corr_timestamp', 'orig_mh1_1',  'orig_mh2_1',  'orig_mhz_1',  'orig_mh1_2',  'orig_mh2_2',  'orig_mhz_2',  'orig_mh1_3',  'orig_mh2_3',  'orig_mhz_3',  'orig_mh1_4',  'orig_mh2_4',  'orig_mhz_4']].to_numpy()
        # 
        # 
        # t1 = datetime.now()
        # df.reset_index(inplace=True)
        # df = df.rename(columns = {'index':'old_index'})
        # df.set_index(['corr_timestamp'], drop=False, append=False, inplace=True, verify_integrity=False)
        # df.sort_index()
        # logging.info(df.iloc[58:92].to_string())

        # logging.info(df.head().to_string())
        # found  = 0
        # for i in range(0,len(df) -1):
        # # for i in range(68,69):
        #     corr_timestamp_low = array[i,0]
        #     corr_timestamp_high = array[i,1]
        #     frame = array[i,2]
        #     orig_station = array[i,3]
        #     corr_timestamp = array[i,4]
        #     orig_mh1_1 = array[i,5]
        #     orig_mh2_1 = array[i,6]  
        #     orig_mhz_1 = array[i,7]  
        #     orig_mh1_2 = array[i,8] 
        #     orig_mh2_2 = array[i,9] 
        #     orig_mhz_2 = array[i,10]
        #     orig_mh1_3 = array[i,11] 
        #     orig_mh2_3 = array[i,12] 
        #     orig_mhz_3 = array[i,13] 
        #     orig_mh1_4 = array[i,14] 
        #     orig_mh2_4 = array[i,15]  
        #     orig_mhz_4 = array[i,16]
        # 
        #     # logging.info('{} {} {} {}'.format(corr_timestamp_low,corr_timestamp_high,frame,orig_station))
        #     # lst = df[(df.index > corr_timestamp_low) & (df.index < corr_timestamp_high) & (df.frame == frame) & (df.orig_station != orig_station)].old_index.tolist()
        #     lst = df[(df.corr_timestamp > corr_timestamp_low) & (df.corr_timestamp < corr_timestamp_high) & (df.frame == frame) & (df.orig_station != orig_station) & 
        #         (df.orig_mh1_1 == orig_mh1_1) &
        #         (df.orig_mh2_1 == orig_mh2_1) &
        #         (df.orig_mhz_1 == orig_mhz_1) &
        #         (df.orig_mh1_2 == orig_mh1_2) &
        #         (df.orig_mh2_2 == orig_mh2_2) &
        #         (df.orig_mhz_2 == orig_mhz_2) &
        #         (df.orig_mh1_3 == orig_mh1_3) &
        #         (df.orig_mh2_3 == orig_mh2_3) &
        #         (df.orig_mhz_3 == orig_mhz_3) &
        #         (df.orig_mh1_4 == orig_mh1_4) &
        #         (df.orig_mh2_4 == orig_mh2_4) &
        #         (df.orig_mhz_4 == orig_mhz_4)
        #         ].index.tolist()
        #     # lst = df[(df['frame'] == frame)].index.tolist()
        #     # print(lst)
        #     # lst = df[(df.corr_timestamp > corr_timestamp_low)].index.tolist()
        # 
        #     if len(lst) > 0:
        #         logging.info(df.iloc[lst].to_string())
        #         # logging.info('found')
        #         # logging.info(df.iloc[lst].to_string())
        #         # print(lst)
        #         found += 1
        # # logging.info('found {}'.format(found))
        # print(found)
        # 
        # 
        # t2 = datetime.now()
        # logging.info('new : {}'.format((t2-t1).total_seconds()))
        # print('new ', (t2-t1).total_seconds())


        # new  0.23747
        # df['Draws'] = search_copied_station(df['corr_timestamp_low'].values,df['corr_timestamp_high'].values,df['frame'].values,df['orig_station'].values)






# def search_copied_station(corr_timestamp_low,corr_timestamp_high,frame,orig_station):
#     # df['Draws'] = 0
#     df.loc[('orig_station' == orig_station), 'found'] = 1



# XXXX
    def test_find_good_records(self):

        config.cumsum_test = 4

        ####################
        # simple test, include everything
        lst = list(range(1,8))
        
        low, high = find_good_records(lst)
        lists = make_good_lists(lst,low,high)
        bad_lists = make_bad_lists(lst,low,high)

        self.assertEqual(low, [0])  
        self.assertEqual(high, [6])  
        self.assertEqual([lst], lists)  
        self.assertEqual([], bad_lists)  
        
        ####################
        # two valid lists
        lst1 = list(range(1,8))
        lst2 = list(range(0,8))
        lst = []
        lst.extend(lst1)
        lst.extend(lst2)
        
        low, high = find_good_records(lst)
        lists = make_good_lists(lst,low,high)
        bad_lists = make_bad_lists(lst,low,high)
        
        self.assertEqual(low, [0,7])  
        self.assertEqual(high, [6,14])  
        self.assertEqual(lst1, lists[0])  
        self.assertEqual(lst2, lists[1])  
        self.assertEqual([], bad_lists)  


        
        ####################
        # two valid lists with a gap between them 
        lst1 = list(range(0,8))
        lstx = [0,0,0,1,0,1,2]
        lst2 = list(range(0,8))
        lst = []
        lst.extend(lst1)
        lst.extend(lstx)
        lst.extend(lst2)
        
        low, high = find_good_records(lst)
        print(low)
        print(high)
        
        
        lists = make_good_lists(lst,low,high)
        bad_lists = make_bad_lists(lst,low,high)
        
        self.assertEqual(low, [0,15])  
        self.assertEqual(high, [7,22])  
        self.assertEqual(lst1, lists[0])  
        self.assertEqual(lst2, lists[1])  
        self.assertEqual(lstx, bad_lists[0])  

        # XXXX
        ####################
        # nothing valid
        lstx = [0,0,0,1,0,1,2]
        lst = []
        lst.extend(lstx)
        
        low, high = find_good_records(lst)
        lists = make_good_lists(lst,low,high)
        bad_lists = make_bad_lists(lst,low,high)
        
        # everything should be empty 
        self.assertEqual(low, [])  
        self.assertEqual(high, [])  
        self.assertEqual(lists, [])  
        # the bad list should be full
        self.assertEqual(lstx, bad_lists[0])  


        ####################
        # list at the end is invalid
        lst1 = list(range(1,8))
        lstx = [0,0,0,1,0,1,2]
        lst = []
        lst.extend(lst1)
        lst.extend(lstx)

        low, high = find_good_records(lst)
        lists = make_good_lists(lst,low,high)
        bad_lists = make_bad_lists(lst,low,high)

        self.assertEqual(low, [0])  
        self.assertEqual(high, [6])  
        self.assertEqual([lst1], lists)   
        self.assertEqual(lstx, bad_lists[0])  



        ####################
        # list at the beginning is invalid
        lstx = [0,0,0,1,0,1,2]
        lst1 = list(range(0,8))
        lst = []
        lst.extend(lstx)
        lst.extend(lst1)

        low, high = find_good_records(lst)
        lists = make_good_lists(lst,low,high)
        bad_lists = make_bad_lists(lst,low,high)

        self.assertEqual(low, [7])  
        self.assertEqual(high, [14])  
        self.assertEqual([lst1], lists) 
        self.assertEqual(lstx, bad_lists[0])  

        ####################
        # three valid lists
        lst1 = list(range(1,8))
        lst2 = list(range(0,8))
        lst3 = list(range(0,8))
        lst = []
        lst.extend(lst1)
        lst.extend(lst2)
        lst.extend(lst3)

        low, high = find_good_records(lst)
        lists = make_good_lists(lst,low,high)
        bad_lists = make_bad_lists(lst,low,high)
        
        self.assertEqual(lst1, lists[0])  
        self.assertEqual(lst2, lists[1])   
        self.assertEqual(lst3, lists[2])     
        # nothing is invalid
        self.assertEqual([], bad_lists)  

    def test_read_in_local(self):
        default_config()


        # working with some nasty data 


        df = pd.read_csv('/Users/cnunn/lunar_data/PDART_CSV_WORK_TAPES_FULL/wtn.13.15.csv.gz', dtype=str)

        df = initial_cleanup_local(df)


def make_good_lists(lst,low,high):
    lists = []
    for l, h in zip(low,high):
        l1 = lst[l:h+1]
        lists.append(l1)

    return lists 

def make_bad_lists(lst,low,high):
    lists = []
    low_bad = low.copy()
    high_bad = high.copy()

    # print('low_bad ', low_bad)
    # print('high_bad ', high_bad)
    # 
    # if the first good record is greater than zero, then the first records 
    # contain errors

    high_bad.insert(0,-1)
    low_bad.append(len(lst))

    # print('low_bad ', low_bad)
    # print('high_bad ', high_bad)

    # add any records between good records to the bad list
    for h, l in zip(high_bad,low_bad):
        # print(h, l)
        l1 = lst[h+1:l]
        # print(l1)
        if len(l1) > 0:
            lists.append(l1)

    return lists 

def find_good_records(lst):

    # search backwards for first value greater than cumsum_test
    end_idx = len(lst) - 1
    high = []
    low = []

    while True:

        idx = None
        try:
            # idx_backward = list[end_idx::-1].index(>config.cumsum_test)
            
            # print('end_idx ', end_idx)
            # print(lst)
            idx_backward = next(x[0] for x in enumerate(lst[end_idx::-1]) if x[1] >= config.cumsum_test)
            # print('idx_backward ', idx_backward)
            idx = end_idx - idx_backward
            # print('idx ', idx)
            high.append(idx)
            # print('value ', lst[idx])
        except StopIteration:
            # didn't find any high values 
            pass

        low1 = None
        # if a value larger than the cumsum test was found, now find the first zero 
        if idx is not None:
            end_idx = idx
            low1 = 0
            try:
                # idx_backward = list[end_idx::-1].index(>config.cumsum_test)
                idx_backward = next(x[0] for x in enumerate(lst[end_idx::-1]) if x[1] == 0)
                # print('idx_backward ', idx_backward)
                low1 = end_idx - idx_backward
            except StopIteration:
                # didn't find a zero (this can be OK when it is the first record)
                # print('did not find the low')
                pass
            low.append(low1)

        # print(' low 1 at the end ', low1)
        if low1 is not None:
            end_idx = low1
        else:
            break

        if end_idx == 0:
            break

    low.reverse()
    high.reverse()

    return low, high

def suite():

    logging.basicConfig(filename='logs/import_test.log', filemode='w', 
      level=logging.DEBUG,format='%(message)s')

    # copy these to the right place in the code to view the output: 
    # df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 
    # logging.info('Before\n{}'.format(df_gst.to_string()))
    # logging.info('After\n{}'.format(df_gst.to_string()))

    suite = unittest.TestSuite()


    # run all the tests 
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test'))

    # run a single test
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_isolate'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_compliance'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_station_fix_bad_timestamps'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_consec_frame'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_frame_stuff'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_breaks'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_frame_change'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_ground_station'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_damaged_timestamps'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_damaged_frames'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_missing'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_begin_end'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_nasty'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_read_in'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_row_duplicatons'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_find_good_records'))

    suite.addTest(unittest.makeSuite(ImportTestCase, 'test_read_in_local'))



    return suite

def initial_cleanup_local(df):

    # df['orig_timestamp'] = df['orig_timestamp'].astype('datetime64[ns, UTC]')
    # 
    # print('a')

    # # # # main tapes extracted without the S 
    # print('b')
    df['orig_station'].replace('12', 'S12',inplace=True)
    print('a')
    df['orig_station'].replace('14', 'S14',inplace=True)
    # print('c')
    df['orig_station'].replace('15', 'S15',inplace=True)
    # print('d')
    df['orig_station'].replace('16', 'S16',inplace=True)
    # print('e')    
    df['orig_station'] = df['orig_station'].astype('string')

    # if 'bit_synchronizer' in df.columns:
    #     df['bit_synchronizer'] = df['bit_synchronizer'].astype('string')
    #     print('f')
    # else:
    #     df['bit_synchronizer'] = pd.NA
    #     print('g')
    # df['sync'] = df['sync'].astype('string')
    # print('h')

    return df

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
