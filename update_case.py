#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

:copyright:
    The pdart Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)

Update the locations of streams with flat response files from 00 to 01.
There is no problem with rerunning the method - it will read in files 
with either location and write out files with the correct one. 



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

from pdart.view import find_filename_date_lower, find_dir_lower


import natsort # if required, install with pip install natsort

DELTA = 0.1509433962

from pdart.initial_checks import _est_frame
# from pdart.save_24_hours import update_starttimes
import pdart.config as config
# from pdart.csv_import_work_tapes import find_output_dir, make_output_dir, make_filelist
import matplotlib.pyplot as plt

flat_times = [
['S12', UTCDateTime('1974-10-16T14:02:36.073000Z') , UTCDateTime('1975-04-09T15:31:03.702000Z')],
['S12', UTCDateTime('1975-06-28T13:48:23.124000Z') , UTCDateTime('1977-03-27T15:41:06.247000Z')],
['S14', UTCDateTime('1976-09-18T08:24:35.026000Z') , UTCDateTime('1976-11-17T15:34:34.524000Z')],
['S14', UTCDateTime('1976-05-27T08:24:35.026000Z') , UTCDateTime('1976-11-17T15:34:34.524000Z')],
['S15', UTCDateTime('1971-10-24T20:58:47.248000Z') , UTCDateTime('1971-11-08T00:34:39.747000Z')],
['S15', UTCDateTime('1975-06-28T14:36:33.034000Z') , UTCDateTime('1977-03-27T15:24:05.361000Z')],
['S16', UTCDateTime('1972-05-13T14:08:03.157000Z') , UTCDateTime('1972-05-14T14:47:08.185000Z')],
['S16', UTCDateTime('1975-06-29T02:46:45.610000Z') , UTCDateTime('1977-03-26T14:52:05.483000Z')],
]

# flat_times = [
# # test only!
# ['S12', UTCDateTime('1970-10-16T14:02:36.073000Z') , UTCDateTime('1970-04-09T15:31:03.702000Z')],
# # ['S12', UTCDateTime('1975-06-28T13:48:23.124000Z') , UTCDateTime('1977-03-27T15:41:06.247000Z')],
# ]

# flat_times = [
# ['S12', UTCDateTime('1974-10-16T14:02:36.073000Z') , UTCDateTime('1975-04-09T15:31:03.702000Z')],
# # ['S12', UTCDateTime('1975-06-28T13:48:23.124000Z') , UTCDateTime('1977-03-27T15:41:06.247000Z')],
# # ['S14', UTCDateTime('1976-09-18T08:24:35.026000Z') , UTCDateTime('1976-11-17T15:34:34.524000Z')],
# # ['S14', UTCDateTime('1976-05-27T08:24:35.026000Z') , UTCDateTime('1976-11-17T15:34:34.524000Z')],
# # ['S15', UTCDateTime('1971-10-24T20:58:47.248000Z') , UTCDateTime('1971-11-08T00:34:39.747000Z')],
# # ['S15', UTCDateTime('1975-06-28T14:36:33.034000Z') , UTCDateTime('1977-03-27T15:24:05.361000Z')],
# # ['S16', UTCDateTime('1972-05-13T14:08:03.157000Z') , UTCDateTime('1972-05-14T14:47:08.185000Z')],
# # ['S16', UTCDateTime('1975-06-29T02:46:45.610000Z') , UTCDateTime('1977-03-26T14:52:05.483000Z')],
# ]

PEAKED = '00'
FLAT = '01'

from datetime import date, timedelta

def daterange(start_date, end_date):
    start_date = start_date.date
    end_date = end_date.date
    for n in range(int((end_date - start_date).days)):
        # print(start_date + timedelta(n))
        out_date = UTCDateTime(start_date + timedelta(n)) 
        # print(start_date + timedelta(n), out_date)
        yield out_date



def update_case(
    top_level_dir='.',
    # processed_dir='.',
    log_dir='.',
    # filenames=None,
    # # logging_level=logging.DEBUG
    logging_level=logging.INFO,
    # single_station=None,
    # single_ground_station=None
    ):

    log_filename = 'update_case.log'
    log_filename = os.path.join(log_dir,log_filename)
    logging.basicConfig(filename=log_filename, filemode='w', 
      level=logging_level,format='%(message)s')
    print('log file ', log_filename)

    '''
    '''
    # st = read('/Users/cnunn/lunar_data/tmp_PDART/s12/1974/290/xa.s12.*.mhz.1974.290.0.mseed')
    # 
    # for tr in st:
    #     print(tr.stats)
    # print('')
    # # print(st)
    # # 
    # # 
    # # 
    # st = read('/Users/cnunn/lunar_data/tmp_PDART/s12/1974/290/xa.s12.*.mh1.1974.290.0.mseed')
    # for tr in st:
    #     print(tr.stats)
    # # # st = read('/Users/cnunn/lunar_data/tmp_PDART/s12/1974/290/xa.s12.01.mh2.1974.290.0.mseed')
    # # # print(st)
    # # 
    # exit()

    # this code ignores what was there before - so it needs to retrieve both 00 and 01 locations


    wildcard_files = os.path.join(top_level_dir,'*/*/*/xa.*.*.*.*.*.0.mseed')
    print(wildcard_files)

    for fullpath_old  in glob.glob(wildcard_files):
        file_old = os.path.basename(fullpath_old)
        dir = os.path.dirname(fullpath_old)
        
        
        file_new = file_old.upper()
        fullpath_new = os.path.join(dir,file_new)
        # print('renaming ', fullpath_old, fullpath_new)      
        os.rename(fullpath_old,fullpath_new)

if __name__ == "__main__":
    top_level_dir = '/Users/cnunn/lunar_data/PDART'
    log_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'
    update_case(top_level_dir=top_level_dir,log_dir=log_dir)
