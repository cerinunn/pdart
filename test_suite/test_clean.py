#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
Test Suite 
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
from obspy.core import Stream, read



# Qt5Agg seems to work best on Mac - try 'TkAgg' if that works for you
# put this after the other imports, otherwise it can be overridden
import matplotlib  
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
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
from pdart.test_suite.test_import import basic_timeseries_station as basic_timeseries_station
from pdart.test_suite.test_import import default_config as default_config
from pdart.csv_join_work_tapes import stream_import, initial_cleanup, merge_channel_stream, despike3, loose_frame_diff, read_file, process_list
from pdart.csv_check_work_tapes import calculate_gaps, to_Int64, add_or_minus_frame
from pdart.extra_plots.plot_timing_divergence import plot_timing
from pdart.snippet import relative_timing_trace

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

def local_initial_cleanup(df):

    # df['orig_timestamp'] = df['orig_timestamp'].astype('datetime64[ns, UTC]')

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
    df['frame'] = to_Int64(df['frame'])

    # replace empty ground station and orig station with values
    df['orig_ground_station'].fillna(-1,inplace=True)
    df['orig_station'].fillna('S0',inplace=True)
    df['orig_frame'].fillna(-99,inplace=True)
    df['frame'].fillna(-99,inplace=True)

    df['corr_timestamp'] = df['corr_timestamp'].astype('datetime64[ns, UTC]')
    df['corr_ground_station'] = to_Int64(df['corr_ground_station'])

    df['time_index'] = 0

    # add new columns for timestamps
    # df['corr_timestamp'] = df['orig_timestamp']   
    # df['corr_gap'] = np.NaN
    # df['frame_change'] = 1
    # df['corr_gap_count'] = 1
    # df['corr_ground_station'] = df['orig_ground_station']
    # df['cumsum'] = 0
    # 
    # df['begin'] = False
    # df['end'] = False

    return df

def overlapping_timeseries():
    ##############################
    df_gst = basic_timeseries_station()  
    starttime0 = df_gst.corr_timestamp.iloc[0]
    starttime0 = pd.Timestamp(starttime0) 
    df_gst['time_index'] = np.arange(len(df_gst))
    df_gst2 = df_gst.copy(deep=True)

    # make the first one a bit shorter
    df_gst.drop(df_gst.index[30:],inplace=True)

    # make the second one a bit shorter, with some overlap
    df_gst2.drop(df_gst2.index[:20],inplace=True)
    df_gst2.reset_index(inplace=True,drop=True)
    df_gst2.corr_timestamp = df_gst2.corr_timestamp + pd.Timedelta(seconds=0.802)
    df_gst2.orig_timestamp = df_gst2.orig_timestamp + pd.Timedelta(seconds=0.802)
    df_gst2.corr_ground_station = 7

    stream1 = stream_import(df_gst,starttime0=starttime0,index0=0) 
    stream2 = stream_import(df_gst2,starttime0=starttime0,index0=0) 

    return(stream1, stream2, starttime0)

def basic_timeseries_stream():
    df_gst = basic_timeseries_station()
    starttime0 = df_gst.corr_timestamp.iloc[0]
    starttime0 = pd.Timestamp(starttime0)  
    df_gst['time_index'] = np.arange(len(df_gst)) 
    stream = stream_import(df_gst,starttime0=starttime0,index0=0) 
    return stream

def overlapping_timeseries_five():
    ##############################
    df_gst = basic_timeseries_station()  
    starttime0 = df_gst.corr_timestamp.iloc[0]
    starttime0 = pd.Timestamp(starttime0) 
    df_gst['time_index'] = np.arange(len(df_gst))
    df_gst2 = df_gst.copy(deep=True)
    df_gst3 = df_gst.copy(deep=True)
    df_gst4 = df_gst.copy(deep=True)
    df_gst5 = df_gst.copy(deep=True)

    df_gst.corr_ground_station = 1

    # first timeseries - rows 0:29
    df_gst.drop(df_gst.index[30:],inplace=True)


    # second timeseries - rows 20:39 (some overlap with the first)
    df_gst2.drop(df_gst2.index[:20],inplace=True)
    df_gst2.drop(df_gst2.index[40-20:],inplace=True)
    df_gst2.reset_index(inplace=True,drop=True)
    # df_gst2.corr_timestamp = df_gst2.corr_timestamp + pd.Timedelta(seconds=0.802)
    # df_gst2.orig_timestamp = df_gst2.orig_timestamp + pd.Timedelta(seconds=0.802)
    df_gst2.corr_ground_station = 2

    # third timeseries - rows 50:69 (no overlap with the second)
    df_gst3.drop(df_gst3.index[:50],inplace=True)
    df_gst3.drop(df_gst3.index[70-50:],inplace=True)
    df_gst3.reset_index(inplace=True,drop=True)
    # df_gst3.corr_timestamp = df_gst3.corr_timestamp + pd.Timedelta(seconds-0.702)
    # df_gst3.orig_timestamp = df_gst3.orig_timestamp + pd.Timedelta(seconds-0.702)
    df_gst3.corr_ground_station = 3

    # fourth timeseries - rows 60:64 (contained in the third)
    df_gst4.drop(df_gst4.index[:60],inplace=True)
    df_gst4.drop(df_gst4.index[65-60:],inplace=True)
    df_gst4.reset_index(inplace=True,drop=True)
    # df_gst4.corr_timestamp = df_gst4.corr_timestamp + pd.Timedelta(seconds=0.804)
    # df_gst4.orig_timestamp = df_gst4.orig_timestamp + pd.Timedelta(seconds=0.804)
    df_gst4.corr_ground_station = 4

    # fifth timeseries - rows 62:71 (overlap with third and fourth)
    df_gst5.drop(df_gst5.index[:62],inplace=True)
    logging.info(df_gst5.to_string())
    df_gst5.reset_index(inplace=True,drop=True)
    # df_gst5.corr_timestamp = df_gst2.corr_timestamp + pd.Timedelta(seconds=-0.803)
    # df_gst5.orig_timestamp = df_gst2.orig_timestamp + pd.Timedelta(seconds=-0.803)
    df_gst5.corr_ground_station = 5
    

    stream1 = stream_import(df_gst,starttime0=starttime0,index0=0) 
    stream2 = stream_import(df_gst2,starttime0=starttime0,index0=0) 
    stream3 = stream_import(df_gst3,starttime0=starttime0,index0=0) 
    stream4 = stream_import(df_gst4,starttime0=starttime0,index0=0) 
    stream5 = stream_import(df_gst5,starttime0=starttime0,index0=0) 

    return stream1, stream2,stream3, stream4, stream5, starttime0

def basic_timeseries_station_good_timing(no_records=6):

    df_gst = basic_timeseries_station(no_records)  

    orig_timestamps = []
    orig_timestamps.append(df_gst['orig_timestamp'].iloc[0])
    
    for i in range(1, len(df_gst)):
        orig_timestamps.append(orig_timestamps[-1] + pd.Timedelta(seconds=0.604))
    orig_timestamps = [x.strftime(DATEFORMAT) for x in orig_timestamps]
    df_gst['orig_timestamp'] = orig_timestamps
    df_gst['orig_timestamp'] = pd.to_datetime(orig_timestamps)
    df_gst['corr_timestamp'] = df_gst['orig_timestamp']

    x = [x/180 for x in range(0,4801*4)]
    sine = 20* np.sin(x)

    st = read()
    st = st.select(channel='EHZ')

    mhz1 = sine[0:4800*4:4]
    mhz2 = sine[1:4800*4:4]
    mhz3 = sine[2:4800*4:4]
    mhz4 = sine[3:4800*4:4]

    df_gst['orig_mhz_1'] = mhz1
    df_gst['orig_mhz_2'] = mhz2
    df_gst['orig_mhz_3'] = mhz3
    df_gst['orig_mhz_4'] = mhz4

    # print(len(mhz4))
    # print(mhz1[0:4])
    # # mhz1 = [x.strftime(DATEFORMAT) for x in orig_timestamps]
    # # print(len(st[0].data)*4)
    # print(len(df_gst))
    # exit() 
    # pd.to_datetime(df_gst['corr_timestamp'])
    # pd.to_datetime(df_gst['orig_timestamp'])

    # print(df_gst.tail().to_string())
    
    # print(df_gst['orig_timestamp'].dtypes)
    
    # logging.info(df_gst.to_string())

    df_gst['delta4'] = 0.604

    return df_gst

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


    def test_combine1(self):
        default_config() 
        config.clean_spikes = True

        # frame_change = loose_frame_diff(21,84,2699.741)
        # self.assertEqual(frame_change,4473)

        # +1 change
        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(last_frame=44,new_frame=45,gap=1)
        self.assertEqual(actual_frame_gap,1)
        self.assertEqual(frame_gap_correct,True)    
        
        # +1 change
        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(last_frame=44,new_frame=45,gap=-0.3)
        self.assertEqual(actual_frame_gap,1)
        self.assertEqual(frame_gap_correct,True)  

        # frame later, time just a bit earlier
        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(44,45,-0.5)
        self.assertEqual(actual_frame_gap,1)
        self.assertEqual(frame_gap_correct,False)

        # +2 change
        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(last_frame=44,new_frame=46,gap=0.7)
        self.assertEqual(actual_frame_gap,2)
        self.assertEqual(frame_gap_correct,True)

        # +2 change
        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(last_frame=44,new_frame=46,gap=-0.7)
        self.assertEqual(actual_frame_gap,2)
        self.assertEqual(frame_gap_correct,False)


        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(44,45,-0.1)
        self.assertEqual(actual_frame_gap,1)
        self.assertEqual(frame_gap_correct,True)

        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(44,45,0)
        self.assertEqual(actual_frame_gap,1)
        self.assertEqual(frame_gap_correct,True)
        
        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(44,45,0.1)
        self.assertEqual(actual_frame_gap,1)
        self.assertEqual(frame_gap_correct,True)

        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(last_frame=27,new_frame=38,gap=60.988)
        self.assertEqual(actual_frame_gap, 101)
        self.assertEqual(frame_gap_correct,True)

        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(last_frame=10,new_frame=11,gap=1*0.6038)  
        self.assertEqual(actual_frame_gap, 1)
        self.assertEqual(frame_gap_correct,True)

        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(last_frame=10,new_frame=30,gap=20*0.6038)
        self.assertEqual(actual_frame_gap, 20)
        self.assertEqual(frame_gap_correct,True)

        frame_gap_correct, actual_frame_gap, absolute_error, percent_error = loose_frame_diff(last_frame=30,new_frame=10,gap=70*0.6038)
        self.assertEqual(actual_frame_gap, 70)
        self.assertEqual(frame_gap_correct,True)



        print('End of test_combine1')

    def test_despike3(self):

        default_config() 
        config.clean_spikes = True
    
        TEST_SPIKE_V_SMALL = 5
        TEST_SPIKE_SMALL = 10
        TEST_SPIKE_HIGH = 30

        ##############################
        # test the spike types - large
        ##############################
        # Complicated spike with a base above the mean  
        stream = basic_timeseries_stream()   
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        
        # set a mean
        for row_i in range(0,288):
            channel_stream1[0].data[row_i] = 512
        # spike
        for row_i in range(144,145):
            channel_stream1[0].data[row_i] = 512+4+TEST_SPIKE_HIGH+1
        # base to the left 
        for row_i in range(140,144):
            channel_stream1[0].data[row_i] = 512+4
        # base to the right 
        for row_i in range(145,148):
            channel_stream1[0].data[row_i] = 512-4
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[90]))
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[170]))
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        
        # set a mean
        for row_i in range(0,288):
            channel_stream1[0].data[row_i] = 512
        # spike
        for row_i in range(144,145):
            channel_stream1[0].data[row_i] = 512-(TEST_SPIKE_HIGH+1)
        # base to the left 
        for row_i in range(140,144):
            channel_stream1[0].data[row_i] = 512
        # base to the right 
        for row_i in range(145,148):
            channel_stream1[0].data[row_i] = 512
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[90]))
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[170]))
        
        
        ##############################
        # test the spike types - very small
        stream = basic_timeseries_stream()
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        
        # set a mean
        for row_i in range(0,288):
            channel_stream1[0].data[row_i] = 512
        # spike
        for row_i in range(144,145):
            channel_stream1[0].data[row_i] = 512-(TEST_SPIKE_V_SMALL+1)
        # base to the left 
        for row_i in range(140,144):
            channel_stream1[0].data[row_i] = 512
        # base to the right 
        for row_i in range(145,148):
            channel_stream1[0].data[row_i] = 512
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[90]))
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[170]))
        
        
        ##############################
        # test the spike types - small    
        stream = basic_timeseries_stream()
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        
        # set a mean
        for row_i in range(0,288):
            channel_stream1[0].data[row_i] = 512
        # spike
        for row_i in range(144,145):
            channel_stream1[0].data[row_i] = 512+2+TEST_SPIKE_SMALL+1
        # base to the left 
        for row_i in range(140,144):
            channel_stream1[0].data[row_i] = 512+2
        # base to the right 
        for row_i in range(145,148):
            channel_stream1[0].data[row_i] = 512-2
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[90]))
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[170]))
        
        ##############################
        # Test removing spike before a masked value
        stream = basic_timeseries_stream()
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        # set a mean   
        for i in range(0,len(channel_stream1[0].data)):
            channel_stream1[0].data[i] = 512
        # allow the data to be masked
        channel_stream1[0].data = np.ma.masked_array(data=channel_stream1[0].data,
                     mask=False)
        
        print('this is before null:')
        # spike index 142
        channel_stream1[0].data[142] = 32
        # mask index 143
        channel_stream1[0].data[143] = np.ma.masked
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[141]))        
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[142]))
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[143]))
        
        ##############################
        # Test removing spike after a masked value
        stream = basic_timeseries_stream()
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        # set a mean   
        for i in range(0,len(channel_stream1[0].data)):
            channel_stream1[0].data[i] = 512
        # allow the data to be masked
        channel_stream1[0].data = np.ma.masked_array(data=channel_stream1[0].data,
                     mask=False)
        
        print('this is after null:')
        # mask index 140
        channel_stream1[0].data[141] = np.ma.masked
        # spike index 141
        channel_stream1[0].data[142] = 32
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))

        print('after null end')
        
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[141]))
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[142]))
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[143]))
        
        ##############################
        # Spike on last record of the trace
        stream = basic_timeseries_stream()
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        # set a mean   
        for i in range(0,len(channel_stream1[0].data)):
            channel_stream1[0].data[i] = 512
        # allow the data to be masked
        channel_stream1[0].data = np.ma.masked_array(data=channel_stream1[0].data,
                     mask=False)
        
        channel_stream1[0].data[-1] = 32
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[-1]))
        
        ##############################
        # Spike on first record of the trace
        stream = basic_timeseries_stream()
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        # set a mean   
        for i in range(0,len(channel_stream1[0].data)):
            channel_stream1[0].data[i] = 512
        # allow the data to be masked
        channel_stream1[0].data = np.ma.masked_array(data=channel_stream1[0].data,
                     mask=False)
        
        channel_stream1[0].data[0] = 32
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[0]))
        
        
        ##############################
        # Test removing spike 2 after masked value 
        stream = basic_timeseries_stream()
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        # set a mean   
        for i in range(0,len(channel_stream1[0].data)):
            channel_stream1[0].data[i] = 512
        # allow the data to be masked
        channel_stream1[0].data = np.ma.masked_array(data=channel_stream1[0].data,
                     mask=False)
        
        print('removing spike 2 after masked value ')
        # mask index 140
        channel_stream1[0].data[140] = np.ma.masked
        # normal on 141
        # spike index 142
        channel_stream1[0].data[142] = 32
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[140]))
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[141]))
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[142]))
        
        ##############################
        # Test removing spike 2 before masked value 
        stream = basic_timeseries_stream()
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        # set a mean   
        for i in range(0,len(channel_stream1[0].data)):
            channel_stream1[0].data[i] = 512
        # allow the data to be masked
        channel_stream1[0].data = np.ma.masked_array(data=channel_stream1[0].data,
                     mask=False)
        
        print('removing spike 2 before masked value ')
        
        # spike index 142
        channel_stream1[0].data[142] = 32
        # normal on 143
        # mask index 144
        channel_stream1[0].data[144] = np.ma.masked
        
        # mask index 144
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[142]))
        self.assertTrue(~np.ma.is_masked(channel_stream1[0].data[143]))
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[144]))
        
        
        
        ##############################
        # Test removing spike after a masked value
        stream = basic_timeseries_stream()
        
        # add an error to the MH1 channel on stream1
        channel_stream1 = stream.select(channel='MH1')
        
        # set a mean   
        for i in range(0,len(channel_stream1[0].data)):
            channel_stream1[0].data[i] = 512
        # allow the data to be masked
        channel_stream1[0].data = np.ma.masked_array(data=channel_stream1[0].data,
                     mask=False)
        

        # -- 512 32 --
        print('another spike near null situation')
        # mask index 140
        channel_stream1[0].data[140] = np.ma.masked
        channel_stream1[0].data[142] = 32
        channel_stream1[0].data[143] = np.ma.masked
        
        channel_stream2 = channel_stream1.copy()
        
        for tr in channel_stream1:
            tr.stats.station = 'despiked'
            rec_spikes_found = despike3(tr)
        
        channel_stream1 += channel_stream2
        # channel_stream1.plot(method='full',size=(1200,600))
        
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[140]))
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[142]))
        self.assertTrue(np.ma.is_masked(channel_stream1[0].data[143]))


        print('End despike3 ')

        

    def test_merge(self):
        default_config() 
        config.clean_spikes = True

        ##############################
        # for this test, ignore the spikes and just test the merge 
        config.clean_spikes = False

        # get two dataframes with overlapping timeseries 
        stream1, stream2, starttime0 = overlapping_timeseries()
        plot_stream = Stream()

        channel_stream1 = stream1.select(channel='MH1')
        # add gap on first trace
        for row_i in range(85,86):
            if np.ma.is_masked(channel_stream1[0].data) == False:    
                channel_stream1[0].data = np.ma.masked_array(channel_stream1[0].data, mask=False)
            channel_stream1[0].data[row_i] = np.ma.masked
        # add error on first trace
        for row_i in range(93,96):
            channel_stream1[0].data[row_i] = 800

        channel_stream2 = stream2.select(channel='MH1')
        # add gap on second trace
        for row_i in range(10,11):
            if np.ma.is_masked(channel_stream2[0].data) == False:    
                channel_stream2[0].data = np.ma.masked_array(channel_stream2[0].data, mask=False)
            channel_stream2[0].data[row_i] = np.ma.masked
        # add error on second trace
        for row_i in range(20,22):
            channel_stream2[0].data[row_i] = 200

        channel_stream = channel_stream1 + channel_stream2

        stream1a = channel_stream1.copy()
        for tr in stream1a:
            tr.stats.location = ''
            # tr.stats.channel = tr.stats.channel + '_orig2b'
        stream2a = channel_stream2.copy()
        for tr in stream2a:
            tr.stats.location = ''
            # tr.stats.channel = tr.stats.channel + '_orig2a'
    
        plot_stream += stream1a
        plot_stream += stream2a    
        
        config.view_corrected_traces = True
        new_channel_stream = merge_channel_stream(channel_stream,delta=DELTA)

        plot_stream += new_channel_stream    
        
        stream = stream1 + stream2 
        for tr in stream:
            tr.stats.location = ''
        
        # add an error to the MH1 channel on stream1

        # add a spike 
        for row_i in range(90,91):
            channel_stream1[0].data[row_i] = 800
        
        # add a gap to the MH1 channel on stream1
        for row_i in range(95,96):
            if np.ma.is_masked(channel_stream1[0].data) == False:    
                channel_stream1[0].data = np.ma.masked_array(channel_stream1[0].data, mask=False)
            channel_stream1[0].data[row_i] = np.ma.masked
        
        # add an error to the MH1 channel on stream2
        channel_stream2 = stream2.select(channel='MH1')
        # add a spike 
        for row_i in range(20,21):
            channel_stream2[0].data[row_i] = 200
        

        for tr in channel_stream2:
            rec_spikes_found = despike3(tr)

        # add a gap to the MH1 channel on stream2
        for row_i in range(25,26):
            if np.ma.is_masked(channel_stream2[0].data) == False:    
                channel_stream2[0].data = np.ma.masked_array(channel_stream2[0].data, mask=False)
            channel_stream2[0].data[row_i] = np.ma.masked
        

        
        plot_stream = channel_stream1 + channel_stream2
        channel_stream = stream.select(channel='MH1')

        plot_stream += channel_stream2
        

        self.assertTrue(~np.ma.is_masked(new_channel_stream[0].data[85]))
        self.assertTrue(~np.ma.is_masked(new_channel_stream[0].data[10+80]))


        self.assertTrue(~np.ma.is_masked(new_channel_stream[0].data[95]))

        # replaced from the good trace 
        self.assertTrue(~(new_channel_stream[0].data[20+80] == 200))
        # TODO THIS IS FAILING
        self.assertTrue(np.ma.is_masked(new_channel_stream[0].data[85]))
        
        plot_stream.plot(method='full',size=(1200,600),automerge=True)

        config.clean_spikes = True




    def test_joins(self):
#         from obspy.core import read
#         st = read()
#         st = st.select(channel='EHZ')
#         st.decimate(10)
# # BW.RJOB..EHZ | 2009-08-24T00:20:03.000000Z - 2009-08-24T00:20:32.990000Z | 100.0 Hz, 3000 samples
#         st.cutout(starttime=UTCDateTime('2009-08-24T00:20:15.000000Z'),endtime=UTCDateTime('2009-08-24T00:20:20.000000Z'))
#         # st.split()
#         print(st)
#         st[1].stats.starttime = st[1].stats.starttime + 0.04
#         print(st)
#         st.merge()
#         print(st)        
# 
#         # st.plot()
# 
#         return 

        default_config() 
        config.clean_spikes = True

        df_gst = basic_timeseries_station_good_timing(400)  
        

        end_sample_time = None
        # #############################
        # # Basic process list - only one file 
        # df_list = []
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst)    
        # 
        # stream = process_list(df_list,'test_data_output')
        # # return 
        # # 
        # # 
        # end_sample_time = stream.select(channel='ATT')[-1].stats.endtime
        # print(end_sample_time)

        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')

        # run it like this so we don't need to rerun it all the time 
        if end_sample_time is None: 
            end_sample_time = UTCDateTime('1976-03-01T00:50:01.962264Z')

        # ###############################
        # # Process list - gap between files - but no errors
        # df_gst1 = df_gst.copy(deep=True)
        # df_gst2 = df_gst.copy(deep=True)
        # 
        # # delete some records 
        # df_gst1.drop(range(1500,len(df_gst1)),inplace=True)
        # df_gst1.reset_index(inplace=True,drop=True)
        # 
        # # delete some records 
        # df_gst2.drop(range(0,2000),inplace=True)
        # df_gst2.reset_index(inplace=True,drop=True)
        # 
        # df_list = []
        # 
        # 
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst1)
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst2)
        # 
        # 
        # stream = process_list(df_list,'test_data_output')
        # end_sample_time_with_gap = stream.select(channel='ATT')[-1].stats.endtime
        # print('end_sample_time_with_gap ', end_sample_time_with_gap)
        # self.assertEqual(end_sample_time,end_sample_time_with_gap)
        # 
        # 
        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')
        # 
        # 
        # ###############################
        # # Process list - gap between files - frame is wrong, but time was fine 
        # df_gst1 = df_gst.copy(deep=True)
        # df_gst2 = df_gst.copy(deep=True)
        # 
        # 
        # # delete some records 
        # df_gst1.drop(range(1500,len(df_gst1)),inplace=True)
        # df_gst1.reset_index(inplace=True,drop=True)
        # 
        # # delete some records 
        # df_gst2.drop(range(0,2000),inplace=True)
        # df_gst2.reset_index(inplace=True,drop=True)
        # 
        # test_framejump = 40
        # for idx in df_gst2.index.tolist():
        #     frame = df_gst2['frame'].iloc[idx]
        #     df_gst2.at[idx,'frame'] = add_or_minus_frame(frame,test_framejump)
        # 
        # # print('sample time ', df_gt2.corr_timestamp.iloc[0])
        # 
        # df_list = []
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst1)
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst2)
        # 
        # print(df_list)
        # 
        # stream = process_list(df_list,'test_data_output')
        # # stream = stream.split()
        # # 
        # # print(stream)
        # # return
        # 
        # # stream.select(channel='ATT').plot()
        # # return
        # 
        # end_sample_time_with_framejump = stream.select(channel='ATT')[-1].stats.endtime
        # print('end_sample_time_with_framejump ', end_sample_time_with_framejump)
        # self.assertEqual(end_sample_time,end_sample_time_with_framejump)
        # 
        # print(stream)
        # 
        # 
        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')
        # 
        # 
        # 
        # ###############################
        # # Process list - gap between files - time is wrong
        # df_gst1 = df_gst.copy(deep=True)
        # df_gst2 = df_gst.copy(deep=True)
        # 
        # # delete some records 
        # df_gst1.drop(range(1500,len(df_gst1)),inplace=True)
        # df_gst1.reset_index(inplace=True,drop=True)
        # 
        # # delete some records 
        # df_gst2.drop(range(0,2000),inplace=True)
        # df_gst2.reset_index(inplace=True,drop=True)
        # 
        # # the second set of records have the wrong time
        # for idx in range(0,len(df_gst2)):
        #     corr_timestamp1 = df_gst2['corr_timestamp'].iloc[idx]
        #     # df_gst2.at[idx,'corr_timestamp'] = corr_timestamp1 + pd.Timedelta(seconds=5.433962+0.1)
        #     df_gst2.at[idx,'corr_timestamp'] = corr_timestamp1 + pd.Timedelta(seconds=0.604*8)
        # df_gst2.orig_timestamp = df_gst2.corr_timestamp
        # 
        # df_list = []
        # 
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst1)
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst2)
        # 
        # stream = process_list(df_list,'test_data_output')
        # end_sample_time_with_gap = stream.select(channel='ATT')[-1].stats.endtime
        # print('end_sample_time_with_gap ', end_sample_time_with_gap)
        # 
        # # the gap should be a straight line, and there should be a severe error in the log 
        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')
        # 
        # ###############################
        # # Process list - small negative gap between frames
        # df_gst1 = df_gst.copy(deep=True)
        # df_gst2 = df_gst.copy(deep=True)
        # 
        # # delete some records 
        # df_gst1.drop(range(1500,len(df_gst1)),inplace=True)
        # df_gst1.reset_index(inplace=True,drop=True)
        # 
        # # delete some records - no gap 
        # df_gst2.drop(range(0,1500),inplace=True)
        # df_gst2.reset_index(inplace=True,drop=True)
        # 
        # # the second set of records have the wrong time
        # for idx in range(0,len(df_gst2)):
        #     corr_timestamp1 = df_gst2['corr_timestamp'].iloc[idx]
        #     # df_gst2.at[idx,'corr_timestamp'] = corr_timestamp1 + pd.Timedelta(seconds=5.433962+0.1)
        #     df_gst2.at[idx,'corr_timestamp'] = corr_timestamp1 + pd.Timedelta(seconds=-0.4)
        # df_gst2.orig_timestamp = df_gst2.corr_timestamp
        # 
        # df_list = []
        # 
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst1)
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst2)
        # 
        # stream = process_list(df_list,'test_data_output')
        # 
        # # the gap should have a small kink - and no error message in the log 
        # # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')
        # 
        # # st = read('test_data_output/s12/1976/061/xa.s12.*.afr.1976.061.*.mseed')
        # # st.plot()
        # 
        # end_sample_time_with_negative = stream.select(channel='ATT')[-1].stats.endtime
        # print('end_sample_time_with_negative ', end_sample_time_with_negative)
        # self.assertEqual(end_sample_time,end_sample_time_with_negative)
        # 
        # ###############################
        # # Process list - three sets, with second set with the wrong time 
        # df_gst1 = df_gst.copy(deep=True)
        # df_gst2 = df_gst.copy(deep=True)
        # df_gst3 = df_gst.copy(deep=True)
        # 
        # # delete some records 
        # df_gst1.drop(range(1500,len(df_gst1)),inplace=True)
        # print(df_gst1.index[0], df_gst1.index[-1])
        # df_gst1.reset_index(inplace=True,drop=True)
        # 
        # # delete some records - 
        # df_gst2.drop(range(0,2000),inplace=True)
        # df_gst2.drop(range(2500,len(df_gst)),inplace=True)
        # print(df_gst2.index[0], df_gst2.index[-1])
        # df_gst2.reset_index(inplace=True,drop=True)
        # 
        # # delete some records - 
        # df_gst3.drop(range(0,3000),inplace=True)
        # print(df_gst3.index[0], df_gst3.index[-1])
        # df_gst3.reset_index(inplace=True,drop=True)
        # 
        # 
        # # the second set of records have the wrong time
        # for idx in range(0,len(df_gst2)):
        #     corr_timestamp1 = df_gst2['corr_timestamp'].iloc[idx]
        #     df_gst2.at[idx,'corr_timestamp'] = corr_timestamp1 + pd.Timedelta(seconds=5.433962+0.1)
        # df_gst2.orig_timestamp = df_gst2.corr_timestamp
        # 
        # df_list = []
        # 
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst1)
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst2)
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst3)
        # 
        # stream = process_list(df_list,'test_data_output')
        # 
        # # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')
        # 
        # # st = read('test_data_output/s12/1976/061/xa.s12.*.afr.1976.061.*.mseed')
        # # st.plot()
        # 
        # end_sample_time_with_third_ok = stream.select(channel='ATT')[-1].stats.endtime
        # print('end_sample_time_with_third_ok ', end_sample_time_with_third_ok)
        # self.assertEqual(end_sample_time,end_sample_time_with_third_ok)
        # 
        # ###############################
        # # Process list - two sets, with second set overlaps the first 
        # df_gst1 = df_gst.copy(deep=True)
        # df_gst2 = df_gst.copy(deep=True)
        # 
        # # delete some records 
        # df_gst1.drop(range(2000,len(df_gst1)),inplace=True)
        # df_gst1.reset_index(inplace=True,drop=True)
        # 
        # # delete some records - but with an overlap
        # df_gst2.drop(range(0,1500),inplace=True)
        # df_gst2.reset_index(inplace=True,drop=True)
        # 
        # logging.info(df_gst.head().to_string())
        # logging.info(df_gst.dtypes)
        # 
        # 
        # # the second set of records have a slight change of time (but within the tolerance)
        # for idx in range(0,len(df_gst2)):
        #     corr_timestamp1 = df_gst2['corr_timestamp'].iloc[idx]
        #     df_gst2.at[idx,'corr_timestamp'] = corr_timestamp1 + pd.Timedelta(seconds=0.1)
        # df_gst2.orig_timestamp = df_gst2.corr_timestamp
        # 
        # df_list = []
        # 
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst1)
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst2)
        # 
        # stream = process_list(df_list,'test_data_output')
        # 
        # # the record plots with a slight kink 
        # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')
        # 
        # # # XXXX    
        # st = read('test_data_output/s12/1976/061/xa.s12.*.mhz.1976.061.*.mseed')
        # st.plot()
        # 
        # 
        # end_sample_time_with_overap = stream.select(channel='ATT')[-1].stats.endtime
        # print('end_sample_time_with_overap ', end_sample_time_with_overap)
        # self.assertEqual(end_sample_time,end_sample_time_with_overap)
        # 
        # ###############################
        # # Process list - two sets, with second set overlaps the first,
        # # with some gaps on left and right 
        # df_gst1 = df_gst.copy(deep=True)
        # df_gst2 = df_gst.copy(deep=True)
        # 
        # # delete some records 
        # df_gst1.drop(range(2000,len(df_gst1)),inplace=True)
        # df_gst1.reset_index(inplace=True,drop=True)
        # 
        # 
        # # delete some records - but with an overlap
        # df_gst2.drop(range(0,1500),inplace=True)
        # df_gst2.reset_index(inplace=True,drop=True)
        # 
        # 
        # 
        # 
        # # the first set have a few gaps
        # for idx in range(1600,1700):
        #     df_gst1.at[idx,'orig_mhz_1'] = pd.NA
        #     df_gst1.at[idx,'orig_mhz_2'] = pd.NA
        #     df_gst1.at[idx,'orig_mhz_3'] = pd.NA
        #     df_gst1.at[idx,'orig_mhz_4'] = pd.NA
        # 
        # # logging.info(df_gst1.to_string())
        # # return
        # 
        # # logging.info()
        # 
        # # the second set have a few gaps 
        # for idx in range(300,400):
        #     df_gst2.at[idx,'orig_mhz_1'] = pd.NA
        #     df_gst2.at[idx,'orig_mhz_2'] = pd.NA
        #     df_gst2.at[idx,'orig_mhz_3'] = pd.NA
        #     df_gst2.at[idx,'orig_mhz_4'] = pd.NA
        # 
        # # the second set of records have a slight change of time (but within the tolerance)
        # for idx in range(0,len(df_gst2)):
        #     corr_timestamp1 = df_gst2['corr_timestamp'].iloc[idx]
        #     df_gst2.at[idx,'corr_timestamp'] = corr_timestamp1 + pd.Timedelta(seconds=0.1)
        # df_gst2.orig_timestamp = df_gst2.corr_timestamp
        # 
        # df_list = []
        # 
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst1)
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst2)
        # 
        # stream = process_list(df_list,'test_data_output')
        # 
        # # the record plots with a slight kink 
        # # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')
        # 
        # # # XXXX    
        # # should look continuous
        # st = read('test_data_output/s12/1976/061/xa.s12.*.mhz.1976.061.*.mseed')
        # st.plot()
        # 
        # 
        # end_sample_time_with_overap = stream.select(channel='ATT')[-1].stats.endtime
        # print('end_sample_time_with_overap ', end_sample_time_with_overap)
        # self.assertEqual(end_sample_time,end_sample_time_with_overap)

        ###############################
        # Process list - two sets, with second set overlaps the first - but with the wrong time 
        # Doesn't try to merge the trace
        df_gst1 = df_gst.copy(deep=True)
        df_gst2 = df_gst.copy(deep=True)
        
        # delete some records 
        df_gst1.drop(range(2000,len(df_gst1)),inplace=True)
        df_gst1.reset_index(inplace=True,drop=True)
        
        # delete some records - but with an overlap
        df_gst2.drop(range(0,1500),inplace=True)
        df_gst2.reset_index(inplace=True,drop=True)
        
        
        # the second set of records have the wrong time 
        for idx in range(0,len(df_gst2)):
            corr_timestamp1 = df_gst2['corr_timestamp'].iloc[idx]
            df_gst2.at[idx,'corr_timestamp'] = corr_timestamp1 + pd.Timedelta(seconds=5.1)
        df_gst2.orig_timestamp = df_gst2.corr_timestamp
        
        df_list = []
        
        read_file('',df_list,join_dir='.',log_dir='.',df=df_gst1)
        read_file('',df_list,join_dir='.',log_dir='.',df=df_gst2)
        
        config.view_corrected_traces = False
        stream = process_list(df_list,'test_data_output')
        
        # the record plots with a slight kink 
        plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')
        
        st = read('test_data_output/s12/1976/061/xa.s12.*.mhz.1976.061.*.mseed')
        st.plot(method='full')
        


        # ###############################
        # # Process list - two sets, with second set overlaps the first - but with the wrong time, but the frames will be close enough to overlap on the wrong time
        # df_gst1 = df_gst.copy(deep=True)
        # df_gst2 = df_gst.copy(deep=True)
        # 
        # # delete some records 
        # df_gst1.drop(range(2000,len(df_gst1)),inplace=True)
        # df_gst1.reset_index(inplace=True,drop=True)
        # 
        # # delete some records - but with an overlap
        # df_gst2.drop(range(0,1500),inplace=True)
        # df_gst2.reset_index(inplace=True,drop=True)
        # 
        # 
        # # the second set of records have the wrong time 
        # for idx in range(0,len(df_gst2)):
        #     corr_timestamp1 = df_gst2['corr_timestamp'].iloc[idx]
        #     df_gst2.at[idx,'corr_timestamp'] = corr_timestamp1 + pd.Timedelta(seconds=54.36)
        # df_gst2.orig_timestamp = df_gst2.corr_timestamp
        # 
        # df_list = []
        # 
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst1)
        # read_file('',df_list,join_dir='.',log_dir='.',df=df_gst2)
        # 
        # config.view_corrected_traces = True
        # config.combine_locations = True
        # 
        # stream = process_list(df_list,'test_data_output')
        # 
        # # plot_timing(top_level_dir='test_data_output',start_time=UTCDateTime('1976-03-01'),stations=['S12'],out_dir='test_data_output')
        # 
        # st = read('test_data_output/s12/1976/061/xa.s12.*.mhz.1976.061.*.mseed')
        # st.plot()
        # # return
        # 
        # end_sample_time_with_overap = stream.select(channel='ATT')[-1].stats.endtime
        # print('end_sample_time_with_overap ', end_sample_time_with_overap)
        # # self.assertEqual(end_sample_time,end_sample_time_with_overap)

    def test_combined_import(self):
        default_config() 
        config.clean_spikes = True

        df_gst = basic_timeseries_station()  
        starttime0 = df_gst.corr_timestamp.iloc[0]
        starttime0 = pd.Timestamp(starttime0)

        ##############################
        # first view a combined import (test stream_import, and check that 
        # the overlaps are correct)
        df_gst = basic_timeseries_station()  
        df_gst['time_index'] = np.arange(len(df_gst))
        
        stream = stream_import(df_gst,starttime0=starttime0,index0=0)
        for tr in stream:
            tr.stats.location = ''
        # stream.plot(method='full',size=(1200,600))

        stream1, stream2, stream3, stream4, stream5, starttime0 = overlapping_timeseries_five()
        stream_overlap = stream1 + stream2 + stream3 + stream4 + stream5

        for tr in stream_overlap:
            tr.stats.location = ''
        # stream_overlap.plot(method='full',size=(1200,600))

        stream_combined = stream + stream_overlap
        stream_combined.plot(method='full',size=(1200,600))


    def test_combined_merge(self):
        default_config() 
        config.clean_spikes = True

        df_gst = basic_timeseries_station()  
        starttime0 = df_gst.corr_timestamp.iloc[0]
        starttime0 = pd.Timestamp(starttime0)  
        
        ##############################
        # Test the merge code (note that this isn't a method, 
        # it's just a copy of the logic in merge_channel_stream())
        data = {'data_x': [501, 1000, pd.NA, 1000, pd.NA, 501], 'data_y': [1000, 501, 500, pd.NA, pd.NA, 501] }
        df_result = pd.DataFrame.from_dict(data)
        df_result['data_x'] = to_Int64(df_result['data_x'])
        df_result['data_y'] = to_Int64(df_result['data_y'])
        
        logging.info('Input')
        logging.info(df_result.to_string())
        
        # fake the mean for this test 
        mean_x = 500
        
        # if the first stream has nulls, replace with non nulls from the second stream
        df_result.loc[df_result['data_x'].isnull(),'data_x'] = df_result['data_y']
        
        count = np.count_nonzero( ~(df_result['data_x'].isna()) & ~(df_result['data_y'].isna()) & (df_result['data_y'] !=  df_result['data_x']))
        
        df_result['data'] = np.where( (~(df_result['data_x'].isna()) & ~(df_result['data_y'].isna()) & ((abs(df_result['data_y'] - mean_x)) < abs(df_result['data_x'] - mean_x))), 
            df_result['data_y'], df_result['data_x'])
        
        logging.info('Output')
        logging.info(df_result.to_string())
        
        # Input
        #    data_x  data_y
        # 0     501    1000
        # 1    1000     501
        # 2    <NA>     500
        # 3    1000    <NA>
        # 4    <NA>    <NA>
        # 5    501     501
        
        # Output
        # data
        # 0     501    1000   501
        # 1    1000     501   501
        # 2     500     500   500
        # 3    1000    <NA>  1000
        # 4    <NA>    <NA>  <NA>
        # 5    501     501   501
        
        # smallest is on the left
        self.assertEqual(df_result['data'].iloc[0],501)
        # smallest is on the right
        self.assertEqual(df_result['data'].iloc[1],501)
        # left was null, so replace by right
        self.assertEqual(df_result['data'].iloc[2],500)
        # only value is from the left
        self.assertEqual(df_result['data'].iloc[3],1000)
        # both values null 
        self.assertTrue(isinstance(df_result['data'].iloc[4], pd._libs.missing.NAType))
        
        
        ##############################
        # get a basic stream for using later 
        df_gst = basic_timeseries_station()  
        df_gst['time_index'] = np.arange(len(df_gst))
        
        basic_stream = stream_import(df_gst,starttime0=starttime0,index0=0)
        for tr in basic_stream:
            tr.stats.location = ''
        
        ##############################
        # Test merging - no errors 
        # by data from an overlapping trace 
        
        # get some streams with overlapping timeseries, and add some errors 
        # (see overlapping_timeseries_five() method for which rows make 
        # up which streams)
        stream1, stream2, stream3, stream4, stream5, starttime0 = overlapping_timeseries_five()
        
        stream_overlap = stream1 + stream2 + stream3 + stream4 + stream5
        overlap_channel_stream = stream_overlap.select(channel='MH2').copy()
        
        merged_channel_stream = merge_channel_stream(overlap_channel_stream.select(channel='MH2'),delta=DELTA)
        
        # plot if required 
        combined_stream = Stream()
        combined_stream += merged_channel_stream
        # combined_stream += basic_stream.select(channel='MH2')
        combined_stream += overlap_channel_stream
        # combined_stream.plot(method='full',size=(1200,600))
        
        end_seconds_merged = merged_channel_stream.select(channel='MH2')[-1].stats.endtime.second
        end_seconds_stream5 = stream5.select(channel='MH2')[-1].stats.endtime.second
        end_seconds_basic = basic_stream.select(channel='MH2')[-1].stats.endtime.second
        self.assertEqual(end_seconds_merged,end_seconds_stream5)
        self.assertEqual(end_seconds_merged,end_seconds_basic)
        
        # check that the values are the same
        l_new = list(merged_channel_stream.select(channel='MH2')[0].data[100:120])
        l_orig = list(basic_stream.select(channel='MH2')[0].data[100:120])
        for a, b in zip(l_new, l_orig):
            self.assertEqual(a,b)
        
        ##############################
        # Test merging - one trace has errors, but is replaced 
        # by data from an overlapping trace 
        
        # get some streams with overlapping timeseries, and add some errors 
        # (see overlapping_timeseries_five() method for which rows make 
        # up which streams)
        stream1, stream2, stream3, stream4, stream5, starttime0 = overlapping_timeseries_five()
        
        # Add some errors 
        for tr in stream1.select(channel='MH2'):
            # add some errors (this one won't be fixed)
            for row_i in range(10,20):
                tr.data[row_i] = 100
        
        for tr in stream1.select(channel='MH2'):
            # add some errors
            for row_i in range(118,120):
                tr.data[row_i] = 700
        
        for tr in stream2.select(channel='MH2'):
            # add some errors
            for row_i in range(20,22):
                tr.data[row_i] = 10
        
        stream_overlap = stream1 + stream2 + stream3 + stream4 + stream5
        overlap_channel_stream = stream_overlap.select(channel='MH2').copy()
        
        merged_channel_stream = merge_channel_stream(overlap_channel_stream.select(channel='MH2'),delta=DELTA)
        
        # plot if required 
        combined_stream = Stream()
        combined_stream += merged_channel_stream
        # combined_stream += basic_stream.select(channel='MH2')
        combined_stream += overlap_channel_stream
        # combined_stream.plot(method='full',size=(1200,600))
        
        end_seconds_merged = merged_channel_stream.select(channel='MH2')[-1].stats.endtime.second
        end_seconds_stream5 = stream5.select(channel='MH2')[-1].stats.endtime.second
        end_seconds_basic = basic_stream.select(channel='MH2')[-1].stats.endtime.second
        self.assertEqual(end_seconds_merged,end_seconds_stream5)
        self.assertEqual(end_seconds_merged,end_seconds_basic)
        
        # check that the merged stream1 has been updated with the correct values 
        # from stream2
        l_new = list(merged_channel_stream.select(channel='MH2')[0].data[110:120])
        l_orig = list(basic_stream.select(channel='MH2')[0].data[110:120])
        for a, b in zip(l_new, l_orig):
            self.assertEqual(a,b)



        ##############################
        # Test removing spike after a masked value
        basic_stream1 = basic_stream.copy()
        
        for tr in basic_stream1.select(channel='MH2'):
            tr.data = np.ma.masked_array(data=tr.data,
                         mask=False)
            # Add a masked value followed by a spike 
            tr.data[9] = np.ma.masked
            tr.data[10] = 100
        
        basic_stream2 = basic_stream1.copy()
        for tr in basic_stream2.select(channel='MH2'):
            # despike it
            rec_spikes_found = despike3(tr)
            tr.stats.station = tr.stats.station + '_despike'
        
        combined_stream = Stream()
        # with spike
        combined_stream += basic_stream1
        # despiked
        combined_stream += basic_stream2
        # combined_stream.select(channel='MH2').plot(method='full',size=(1200,600))
        
        self.assertTrue(np.ma.is_masked(basic_stream2.select(channel='MH2')[0].data[9]))
        self.assertTrue(np.ma.is_masked(basic_stream2.select(channel='MH2')[0].data[10]))




    def test_spike(self):
        default_config() 
        config.clean_spikes = True

        df_gst = basic_timeseries_station()  
        starttime0 = df_gst.corr_timestamp.iloc[0]
        starttime0 = pd.Timestamp(starttime0)

        ##############################
        # first view a basic import
        df_gst = basic_timeseries_station()  
        df_gst['time_index'] = np.arange(len(df_gst))
        
        stream = stream_import(df_gst,starttime0=starttime0,index0=0)
        
        stream.plot(method='full',size=(1200,600))

        ##############################
        # add a single upward spike and don't clean it
        df_gst = basic_timeseries_station()  
        df_gst['time_index'] = np.arange(len(df_gst))
        
        # add a spike 
        spikes = 0
        for row_i in range(25,26):
            df_gst.at[row_i,'orig_mh2_1'] = 800
            spikes += 1
        
        config.clean_spikes = False
        stream = stream_import(df_gst,starttime0=starttime0,index0=0)
        
        for tr in stream:
            if tr.stats.channel == 'MH2' and tr.stats.station == 'S12':
                final = tr.data[row_i*4]
        
        self.assertEqual(final, 800)

        ##############################
        # add a single upward spike and clean it
        df_gst = basic_timeseries_station()  
        df_gst['time_index'] = np.arange(len(df_gst))
        
        # add a spike 
        spikes = 0
        for row_i in range(25,26):
            df_gst.at[row_i,'orig_mh2_1'] = 800
            spikes += 1
        
        config.clean_spikes = True
        stream = stream_import(df_gst,starttime0=starttime0,index0=0)
        
        for tr in stream:
            if tr.stats.channel == 'MH2' and tr.stats.station == 'S12':
                final = tr.data[row_i*4]
        
        # TODO Test is failing 
        self.assertTrue(isinstance(final, np.ma.core.MaskedArray))

        ##############################
        # add a single downward spike
        df_gst = basic_timeseries_station()  
        df_gst['time_index'] = np.arange(len(df_gst))
        
        # add a spike 
        spikes = 0
        for row_i in range(25,26):
            df_gst.at[row_i,'orig_mh2_1'] = 100
            spikes += 1
        
        config.clean_spikes = True
        stream = stream_import(df_gst,starttime0=starttime0,index0=0)
        
        for tr in stream:
            if tr.stats.channel == 'MH2' and tr.stats.station == 'S12':
                final = tr.data[row_i*4]
        
        self.assertTrue(isinstance(final, np.ma.core.MaskedArray))
 
        ##############################
        # add a single downward spike to an array with a few gaps
        df_gst = basic_timeseries_station()  
        
        # add a spike 
        spikes = 0
        for row_i in range(25,26):
            df_gst.at[row_i,'orig_mh1_1'] = 100
            spikes += 1
        
        # delete some records 
        to_drop = [28, 31, 32]
        df_gst.drop(to_drop,inplace=True)
        
        df_gst.reset_index(inplace=True,drop=True)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  
        
        # Skipped timestamps should be added 
        df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
        
        df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 
        df_gst = local_initial_cleanup(df_gst)
        
        df_gst['time_index'] = np.arange(len(df_gst))
        
        config.clean_spikes = True
        stream = stream_import(df_gst,starttime0=starttime0,index0=0)
        # stream.plot(method='full',size=(1200,600))
        
        for tr in stream:
            if tr.stats.channel == 'MH1' and tr.stats.station == 'S12':
                final = tr.data[row_i*4]
        
        self.assertTrue(isinstance(final, np.ma.core.MaskedArray))

# This is not a test, but instead uses the code to remove the spikes 
# TODO - move to extra_plots
def plot_spikes_removal():
    ##############################
    # add a single downward spike to an array with a few gaps
    df_gst = basic_timeseries_station()  
    starttime0 = df_gst.corr_timestamp.iloc[0]
    starttime0 = pd.Timestamp(starttime0) 
    
    # add a spike 
    spikes = 0
    for row_i in range(25,26):
        df_gst.at[row_i,'orig_mh1_1'] = 100
        spikes += 1
    
    # delete some records 
    # to_drop = [28, 31, 32]
    # df_gst.drop(to_drop,inplace=True)
    
    # df_gst.reset_index(inplace=True,drop=True)
    
    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  
    
    # Skipped timestamps should be added 
    # df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
    
    # df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 
    df_gst = local_initial_cleanup(df_gst)
    
    df_gst['time_index'] = np.arange(len(df_gst))
    

    stream = stream_import(df_gst,starttime0=starttime0,index0=0)
    stream_before = stream.select(channel='MH1')
    # replace with a sine wave 
    for tr in stream_before:
        tr.data = (511.2+1.5*np.sin(tr.times())).astype(int)
        # mask the array
        tr.data = np.ma.masked_array(data=tr.data,
                     mask=False)
        tr.data = tr.data[0:201]

    # add a small spike 
    # spikes = 0
    tr.data[24] = 498
        # spikes += 1

    # add a small spike before a gap 
    tr.data[49] = 520
    tr.data[50] = np.ma.masked

    # add a spike after a gap 
    tr.data[75] = np.ma.masked
    tr.data[76] = 498

    # add an up and down spike 
    tr.data[150] = 498 
    tr.data[151] = 520   

    # add a double spike 
    tr.data[175] = 514 
    tr.data[176] = 520   


    # print(tr.data)
    stream_despike = stream_before.copy()
    for tr in stream_despike:
        tr.stats.location = 'despike'

    config.clean_spikes = True
    # if required, clean the spikes 
    if config.clean_spikes:
        for tr in stream_despike:
            # print('trace ', len(tr.data), tr.id)
            despike3(tr)

    total_stream = stream_before + stream_despike

    fig = total_stream.plot(method='full',size=(1200,600),type='relative',handle=True,show=False)
    axes = plt.gcf().get_axes()
    for i, ax in enumerate(axes):
        ax.tick_params(axis='x', labelsize=14 )
        ax.tick_params(axis='y', labelsize=14 )
        if i == 1:
            ax.set_xlabel(xlabel='[s]',size=14)    
        ax.set_ylabel(ylabel='Digital Units',size=14) 
    plt.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.10)
    fig.savefig('spike_cleaning.png')
    plt.show()


        
def suite():



    logging.basicConfig(filename='logs/clean_test.log', filemode='w', 
      level=logging.DEBUG,format='%(message)s')

    suite = unittest.TestSuite()
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test'))

    # add tests individually 
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_combine1'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_despike3'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_merge'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_combined_import'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_combined_merge'))
    # suite.addTest(unittest.makeSuite(ImportTestCase, 'test_spike'))
    suite.addTest(unittest.makeSuite(ImportTestCase, 'test_joins'))



    return suite


if __name__ == '__main__':

    st = read('/Users/cnunn/lunar_data/PDART_CONTINUOUS_MAIN_TAPES/s14/1973/181/xa.s14..shz.1973.181.0.mseed')
    print(st)
    exit()

    # plot_spikes_removal()
    unittest.main(defaultTest='suite')


