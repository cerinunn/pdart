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
from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, Trace, Stats, read

import numpy as np

from obspy.core.utcdatetime import UTCDateTime
from pdart.csv_import_work_tapes import (
  filename_from_trace
)


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

    # def assertDataFrameEqual(self, a, b, msg):
    #     try:
    #         pd.testing.assert_frame_equal(a, b, check_dtype=False)
    #     except AssertionError as e:
    #         try:
    #             pd.testing.assert_frame_equal(a, b, check_dtype=True)
    #         except AssertionError as e:
    #             raise self.failureException(msg) from e
    # 
    # def assertSeriesEqual(self, a, b, msg):
    #     try:
    #         pd.testing.assert_series_equal(a, b)
    #     except AssertionError as e:
    #         raise self.failureException(msg) from e
    # 
    # def setUp(self):
    #     # set up to use pandas to check whether data frames are equal
    #     self.addTypeEqualityFunc(pd.DataFrame, self.assertDataFrameEqual)
    #     # example 
    #     # self.assertEqual(df, df_copy)  
    # 
    #     self.addTypeEqualityFunc(pd.Series, self.assertSeriesEqual)
    #     # example 
    #     # self.assertEqual(df['orig_timestamp'], df_copy['orig_timestamp'])   

    def test_spikes(self):
        


    # def test_save(self):
    #     directory = 'test_data/'
    #     filename = 'xa.s14.06.mhz.1976.066.1.0.mseed'
    #     filepath = os.path.join(directory,filename)
    # 
    #     # # set up the data 
    # 
    #     # stream = read(filepath)
    #     # trace = read(filepath)[0]
    #     # print(trace.stats.starttime)
    #     # 
    #     # stream1 = stream.copy()
    #     # stream1[0].stats.starttime = UTCDateTime('1976-03-06T20:41:06.339000Z') + 3600
    #     # filename1 = 'xa.s14.06.mhz.1976.066.2.0.mseed'
    #     # stream1.write(os.path.join(directory,filename1))
    #     # 
    #     # stream2 = stream.copy()
    #     # stream2[0].stats.starttime = UTCDateTime('1976-03-06T20:41:06.339000Z') + 3600*2
    #     # filename2 = 'xa.s14.06.mhz.1976.066.3.0.mseed'
    #     # stream2.write(os.path.join(directory,filename2))
    #     # 
    #     # stream = read(os.path.join(directory,'xa.s14.06.mhz.1976.066*'))
    #     # print(stream)
    # 
    #     # 1. same => overwrite
    #     # 2. different => increment
    #     # 3. make sure numbers are sequential 
    # 
    # 
    #     trace = read(filepath)[0]
    #     current_filename= filename_from_trace(trace,directory='test_data/',lower=True)
    #     # should be equal
    #     self.assertEqual(current_filename,filename)
    # 
    #     filename1 = 'xa.s14.06.mhz.1976.066.2.0.mseed'
    #     trace1 = read(filepath)[0]
    #     trace1.stats.starttime = UTCDateTime('1976-03-06T20:41:06.339000Z') + 3600
    #     current_filename1= filename_from_trace(trace1,directory='test_data/',lower=True)
    #     # should be equal
    #     self.assertEqual(current_filename1,filename1)
    # 
    #     filename2 = 'xa.s14.06.mhz.1976.066.3.0.mseed'
    #     trace2 = read(filepath)[0]
    #     trace2.stats.starttime = UTCDateTime('1976-03-06T20:41:06.339000Z') + 3600*2
    #     current_filename2= filename_from_trace(trace2,directory='test_data/',lower=True)
    #     # should be equal
    #     self.assertEqual(current_filename2,filename2)
    # 
    #     filename3 = 'xa.s14.06.mhz.1976.066.4.0.mseed'
    #     trace4 = read(filepath)[0]
    #     trace4.stats.starttime = UTCDateTime('1976-03-06T20:41:06.339000Z') + 3600*3
    #     current_filename4= filename_from_trace(trace4,directory='test_data/',lower=True)
    #     # should be equal (makes a new file_no bigger than the other 3)
    #     self.assertEqual(current_filename2,filename2)

def suite():

    logging.basicConfig(filename='logs/import_test.log', filemode='w', 
      level=logging.INFO,format='%(message)s')

    # copy these to the right place in the code to view the output: 
    # df_st1, gaps_long, gaps_8888 = calculate_gaps(df_st1) 
    # logging.info('Before\n{}'.format(df_st1.to_string()))
    # logging.info('After\n{}'.format(df_st1.to_string()))

    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(ImportTestCase, 'test'))
    return suite


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
