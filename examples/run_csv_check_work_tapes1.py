#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run initial checks on csv files - first set of work tapes 

Make checked csv files (will drop some damaged data)
out of the raw csv files extracted from binary. 

Remember to run final_csv_check_main_tapes_all.py AFTER this file, 
which uses slightly different paramters for some of the files which 
didn't have good data recovery.

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from pdart.csv_check_work_tapes import call_csv_check_work_tapes
# from pdart.view import plot_from_stream
from obspy.core.stream import read

import logging
# logging.handlers


import pandas as pd
import numpy as np

def run_csv_check_work_tapes():
    checked_dir_full='/Users/cnunn/lunar_data/PDART_CSV/work_tapes'
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'


    filenames=None

    filenames_all = [
'wtn.1.1.csv.gz', # works?
'wtn.1.10.csv.gz', # Very large clock jumps 
'wtn.1.11.csv.gz',
'wtn.1.12.csv.gz',#1
'wtn.1.13.csv.gz',#2
'wtn.1.14.csv.gz',#3
'wtn.1.15.csv.gz',#this one doesn't work so well - very slow 
'wtn.1.16.csv.gz',
'wtn.1.17.csv.gz',
'wtn.1.18.csv.gz',
'wtn.1.19.csv.gz',
'wtn.1.2.csv.gz',
'wtn.1.20.csv.gz',
'wtn.1.21.csv.gz',
'wtn.1.22.csv.gz',
'wtn.1.23.csv.gz',
'wtn.1.24.csv.gz',
'wtn.1.25.csv.gz',
'wtn.1.26.csv.gz',
'wtn.1.27.csv.gz',
'wtn.1.28.csv.gz',
'wtn.1.29.csv.gz',
'wtn.1.3.csv.gz',
'wtn.1.30.csv.gz',
'wtn.1.31.csv.gz',
'wtn.1.32.csv.gz',
'wtn.1.33.csv.gz',
'wtn.1.34.csv.gz',#checked
'wtn.1.35.csv.gz',
'wtn.1.36.csv.gz',
'wtn.1.37.csv.gz',
'wtn.1.38.csv.gz',
'wtn.1.39.csv.gz',
'wtn.1.4.csv.gz',
'wtn.1.40.csv.gz',
'wtn.1.41.csv.gz',
'wtn.1.42.csv.gz',
'wtn.1.43.csv.gz',
'wtn.1.44.csv.gz',
'wtn.1.45.csv.gz',
'wtn.1.46.csv.gz',
'wtn.1.47.csv.gz',
'wtn.1.48.csv.gz',
'wtn.1.49.csv.gz',
'wtn.1.5.csv.gz',
'wtn.1.50.csv.gz',
'wtn.1.51.csv.gz',
'wtn.1.52.csv.gz',
'wtn.1.53.csv.gz',
'wtn.1.54.csv.gz',
'wtn.1.55.csv.gz',
'wtn.1.6.csv.gz',
'wtn.1.7.csv.gz',
'wtn.1.8.csv.gz',
'wtn.1.9.csv.gz',
'wtn.10.1.csv.gz', #This one has duplicates
'wtn.10.10.csv.gz',
'wtn.10.11.csv.gz',
'wtn.10.12.csv.gz',
'wtn.10.13.csv.gz',
'wtn.10.14.csv.gz',
'wtn.10.15.csv.gz',
'wtn.10.16.csv.gz',
'wtn.10.17.csv.gz',
'wtn.10.18.csv.gz',# a short repeat of frame from earlier in the day, exception handled 
'wtn.10.19.csv.gz',
'wtn.10.2.csv.gz',
'wtn.10.20.csv.gz',
'wtn.10.21.csv.gz',
'wtn.10.22.csv.gz',
'wtn.10.23.csv.gz', # checked, ground station switches a lot, but otherwise fine
'wtn.10.24.csv.gz',
'wtn.10.25.csv.gz',
'wtn.10.26.csv.gz',
'wtn.10.27.csv.gz', # FAILS with exception, need to fix
'wtn.10.28.csv.gz',
'wtn.10.29.csv.gz', # for considerable part of this record, S12 is copied to S14, so S14 is wrong, then same problem with S12 to S15. FAILS with exception, need to fix
'wtn.10.3.csv.gz',
'wtn.10.30.csv.gz',
'wtn.10.31.csv.gz',
'wtn.10.32.csv.gz',
'wtn.10.33.csv.gz',
'wtn.10.34.csv.gz',
'wtn.10.35.csv.gz',
'wtn.10.36.csv.gz',
'wtn.10.37.csv.gz',
'wtn.10.38.csv.gz',
'wtn.10.39.csv.gz',
'wtn.10.4.csv.gz',
'wtn.10.40.csv.gz',
'wtn.10.41.csv.gz',
'wtn.10.42.csv.gz',
'wtn.10.43.csv.gz',
'wtn.10.44.csv.gz',
'wtn.10.45.csv.gz',
'wtn.10.46.csv.gz',
'wtn.10.47.csv.gz',
'wtn.10.48.csv.gz',
'wtn.10.49.csv.gz', #clock flag example
'wtn.10.5.csv.gz',
'wtn.10.50.csv.gz',
'wtn.10.51.csv.gz',
'wtn.10.52.csv.gz',
'wtn.10.53.csv.gz',
'wtn.10.54.csv.gz',#Errors corrected
'wtn.10.55.csv.gz',
'wtn.10.6.csv.gz',
'wtn.10.7.csv.gz',
'wtn.10.8.csv.gz',
'wtn.10.9.csv.gz',
'wtn.11.1.csv.gz',
'wtn.11.10.csv.gz',
'wtn.11.11.csv.gz',
'wtn.11.12.csv.gz',
'wtn.11.13.csv.gz',#Errors corrected, problem with clock flag where it jumps 2, fixed with timestamp_adjust#Problem with single_station='S15',single_ground_station=7, unable to corect 18:17:25.113000. Ground station 9 problem fixed 
'wtn.11.14.csv.gz',
'wtn.11.15.csv.gz',
'wtn.11.16.csv.gz',
'wtn.11.17.csv.gz',#Flat on S12?, but otherwise OK
'wtn.11.18.csv.gz',
'wtn.11.19.csv.gz',
'wtn.11.2.csv.gz',
'wtn.11.20.csv.gz',
'wtn.11.21.csv.gz',
'wtn.11.22.csv.gz',
'wtn.11.23.csv.gz',
'wtn.11.24.csv.gz',
'wtn.11.25.csv.gz',
'wtn.11.26.csv.gz',
'wtn.11.27.csv.gz',
'wtn.11.28.csv.gz',
'wtn.11.29.csv.gz',
'wtn.11.3.csv.gz',
'wtn.11.30.csv.gz',
'wtn.11.31.csv.gz',
'wtn.11.32.csv.gz',
'wtn.11.33.csv.gz',
'wtn.11.34.csv.gz',
'wtn.11.35.csv.gz',
'wtn.11.36.csv.gz',
'wtn.11.37.csv.gz',
'wtn.11.38.csv.gz',
'wtn.11.39.csv.gz',
'wtn.11.4.csv.gz',
'wtn.11.40.csv.gz',
'wtn.11.41.csv.gz',
'wtn.11.42.csv.gz',
'wtn.11.43.csv.gz',
'wtn.11.44.csv.gz',
'wtn.11.45.csv.gz',
'wtn.11.46.csv.gz',
'wtn.11.47.csv.gz',
'wtn.11.48.csv.gz',
'wtn.11.49.csv.gz',
'wtn.11.5.csv.gz',
'wtn.11.50.csv.gz',
'wtn.11.51.csv.gz',
'wtn.11.52.csv.gz',
'wtn.11.53.csv.gz',
'wtn.11.54.csv.gz',
'wtn.11.55.csv.gz',
'wtn.11.6.csv.gz',
'wtn.11.7.csv.gz',#TODO, exception to correct
'wtn.11.8.csv.gz',
'wtn.11.9.csv.gz',
'wtn.12.1.csv.gz',#TODO, exception to correct
'wtn.12.10.csv.gz',
'wtn.12.11.csv.gz',
'wtn.12.12.csv.gz',
'wtn.12.13.csv.gz',
'wtn.12.14.csv.gz',
'wtn.12.15.csv.gz',
'wtn.12.16.csv.gz',
'wtn.12.17.csv.gz',
'wtn.12.18.csv.gz',
'wtn.12.19.csv.gz',
'wtn.12.2.csv.gz',
'wtn.12.20.csv.gz',
'wtn.12.21.csv.gz',
'wtn.12.22.csv.gz',
'wtn.12.23.csv.gz',
'wtn.12.24.csv.gz',
'wtn.12.25.csv.gz',
'wtn.12.26.csv.gz',
'wtn.12.27.csv.gz',
'wtn.12.28.csv.gz',
'wtn.12.29.csv.gz',
'wtn.12.3.csv.gz',
'wtn.12.30.csv.gz',
'wtn.12.31.csv.gz',
'wtn.12.32.csv.gz',
'wtn.12.33.csv.gz',
'wtn.12.34.csv.gz',
'wtn.12.35.csv.gz',
'wtn.12.36.csv.gz',
'wtn.12.37.csv.gz',
'wtn.12.38.csv.gz',
'wtn.12.39.csv.gz',
'wtn.12.4.csv.gz',
'wtn.12.40.csv.gz',
'wtn.12.41.csv.gz',
'wtn.12.42.csv.gz',
'wtn.12.43.csv.gz',
'wtn.12.44.csv.gz',
'wtn.12.45.csv.gz',
'wtn.12.46.csv.gz',
'wtn.12.47.csv.gz',
'wtn.12.48.csv.gz',
'wtn.12.49.csv.gz',
'wtn.12.5.csv.gz',
'wtn.12.50.csv.gz',
'wtn.12.51.csv.gz',
'wtn.12.52.csv.gz',
'wtn.12.53.csv.gz',
'wtn.12.54.csv.gz',
'wtn.12.55.csv.gz',
'wtn.12.6.csv.gz',
'wtn.12.7.csv.gz',
'wtn.12.8.csv.gz',
'wtn.12.9.csv.gz',
'wtn.13.1.csv.gz',#some problems with S14, fixed as far as possible
'wtn.13.10.csv.gz',
'wtn.13.11.csv.gz',
'wtn.13.12.csv.gz',
'wtn.13.13.csv.gz',
'wtn.13.14.csv.gz',
'wtn.13.15.csv.gz',# duplicate section where S14 is copied from S16, also section where the timestamp is the same on two stations, which is fine, the other errors in the log look fine
'wtn.13.16.csv.gz',
'wtn.13.17.csv.gz',
'wtn.13.18.csv.gz',
'wtn.13.19.csv.gz',
'wtn.13.2.csv.gz',
'wtn.13.20.csv.gz',
'wtn.13.21.csv.gz',
'wtn.13.22.csv.gz',
'wtn.13.23.csv.gz',
'wtn.13.24.csv.gz',
'wtn.13.25.csv.gz',
'wtn.13.26.csv.gz',
'wtn.13.27.csv.gz',
'wtn.13.28.csv.gz',
'wtn.13.29.csv.gz',
'wtn.13.3.csv.gz',## duplicate section where S14 is copied from S16
'wtn.13.30.csv.gz',
'wtn.13.31.csv.gz',
'wtn.13.32.csv.gz',
'wtn.13.33.csv.gz',
'wtn.13.34.csv.gz',
'wtn.13.35.csv.gz',
'wtn.13.36.csv.gz',
'wtn.13.37.csv.gz',
'wtn.13.38.csv.gz',
'wtn.13.39.csv.gz',
'wtn.13.4.csv.gz',
'wtn.13.40.csv.gz',
'wtn.13.41.csv.gz',
'wtn.13.42.csv.gz',
'wtn.13.43.csv.gz',
'wtn.13.44.csv.gz',
'wtn.13.45.csv.gz',
'wtn.13.46.csv.gz',
'wtn.13.47.csv.gz',
'wtn.13.48.csv.gz',
'wtn.13.49.csv.gz',
'wtn.13.5.csv.gz',
'wtn.13.50.csv.gz',
'wtn.13.51.csv.gz',
'wtn.13.52.csv.gz',
'wtn.13.53.csv.gz',
'wtn.13.54.csv.gz',
'wtn.13.55.csv.gz',
'wtn.13.56.csv.gz',
'wtn.13.57.csv.gz',
'wtn.13.6.csv.gz',
'wtn.13.7.csv.gz',
'wtn.13.8.csv.gz',
'wtn.13.9.csv.gz',
'wtn.14.1.csv.gz',
'wtn.14.10.csv.gz',
'wtn.14.11.csv.gz',
'wtn.14.12.csv.gz',
'wtn.14.13.csv.gz',
'wtn.14.14.csv.gz',
'wtn.14.15.csv.gz',
'wtn.14.16.csv.gz',
'wtn.14.17.csv.gz',
'wtn.14.18.csv.gz',
'wtn.14.19.csv.gz',
'wtn.14.2.csv.gz',
'wtn.14.20.csv.gz',
'wtn.14.21.csv.gz',
'wtn.14.22.csv.gz',
'wtn.14.23.csv.gz',
'wtn.14.24.csv.gz',
'wtn.14.25.csv.gz',
'wtn.14.26.csv.gz',
'wtn.14.27.csv.gz',
'wtn.14.28.csv.gz',
'wtn.14.29.csv.gz',
'wtn.14.3.csv.gz',
'wtn.14.30.csv.gz',
'wtn.14.31.csv.gz',
'wtn.14.32.csv.gz',
'wtn.14.33.csv.gz',
'wtn.14.34.csv.gz',
'wtn.14.35.csv.gz',
'wtn.14.36.csv.gz',
'wtn.14.37.csv.gz',
'wtn.14.38.csv.gz',
'wtn.14.39.csv.gz',
'wtn.14.4.csv.gz',
'wtn.14.40.csv.gz',
'wtn.14.41.csv.gz',
'wtn.14.42.csv.gz',
'wtn.14.43.csv.gz',
'wtn.14.44.csv.gz',
'wtn.14.45.csv.gz',
'wtn.14.46.csv.gz',
'wtn.14.47.csv.gz',
'wtn.14.48.csv.gz',
'wtn.14.49.csv.gz',
'wtn.14.5.csv.gz',
'wtn.14.50.csv.gz',
'wtn.14.51.csv.gz',
'wtn.14.52.csv.gz',
'wtn.14.53.csv.gz',
'wtn.14.54.csv.gz',
'wtn.14.55.csv.gz',
'wtn.14.6.csv.gz',
'wtn.14.7.csv.gz',
'wtn.14.8.csv.gz',
'wtn.14.9.csv.gz',
'wtn.15.1.csv.gz',
'wtn.15.10.csv.gz',
'wtn.15.11.csv.gz',
'wtn.15.12.csv.gz',
'wtn.15.13.csv.gz',
'wtn.15.14.csv.gz',
'wtn.15.15.csv.gz',
'wtn.15.16.csv.gz',
'wtn.15.17.csv.gz',
'wtn.15.18.csv.gz',
'wtn.15.19.csv.gz',
'wtn.15.2.csv.gz',
'wtn.15.20.csv.gz',
'wtn.15.21.csv.gz',
'wtn.15.22.csv.gz',
'wtn.15.23.csv.gz',
'wtn.15.24.csv.gz',
'wtn.15.25.csv.gz',
'wtn.15.26.csv.gz',
'wtn.15.27.csv.gz',
'wtn.15.28.csv.gz',
'wtn.15.29.csv.gz',
'wtn.15.3.csv.gz',
'wtn.15.30.csv.gz',
'wtn.15.31.csv.gz',
'wtn.15.32.csv.gz',
'wtn.15.33.csv.gz',
'wtn.15.34.csv.gz',
'wtn.15.35.csv.gz',
'wtn.15.36.csv.gz',
'wtn.15.37.csv.gz',
'wtn.15.38.csv.gz',
'wtn.15.39.csv.gz',
'wtn.15.4.csv.gz',
'wtn.15.40.csv.gz',
'wtn.15.41.csv.gz',
'wtn.15.42.csv.gz',
'wtn.15.43.csv.gz',
'wtn.15.44.csv.gz',
'wtn.15.45.csv.gz',
'wtn.15.46.csv.gz',
'wtn.15.47.csv.gz',
'wtn.15.48.csv.gz',
'wtn.15.49.csv.gz',
'wtn.15.5.csv.gz',
'wtn.15.50.csv.gz',
'wtn.15.51.csv.gz',
'wtn.15.52.csv.gz',
'wtn.15.53.csv.gz',
'wtn.15.54.csv.gz',
'wtn.15.55.csv.gz',
'wtn.15.6.csv.gz',

# 'wtn.15.7.csv.gz',
# 'wtn.15.8.csv.gz',
# 'wtn.15.9.csv.gz',
# 'wtn.16.1.csv.gz',
# 'wtn.16.10.csv.gz',
# 'wtn.16.11.csv.gz',
# 'wtn.16.12.csv.gz',
# 'wtn.16.13.csv.gz',
# 'wtn.16.14.csv.gz',
# 'wtn.16.15.csv.gz',
# 'wtn.16.16.csv.gz',
# 'wtn.16.17.csv.gz',
# 'wtn.16.18.csv.gz',
# 'wtn.16.19.csv.gz',
# 'wtn.16.2.csv.gz',
# 'wtn.16.20.csv.gz',
# 'wtn.16.21.csv.gz',
# 'wtn.16.22.csv.gz',
# 'wtn.16.23.csv.gz',
# 'wtn.16.24.csv.gz',
# 'wtn.16.25.csv.gz',
# 'wtn.16.26.csv.gz',
# 'wtn.16.27.csv.gz',
# 'wtn.16.28.csv.gz',
# 'wtn.16.29.csv.gz',
# 'wtn.16.3.csv.gz',
# 'wtn.16.30.csv.gz', # TODO - data not read in correctly, but probably OK 
# 'wtn.16.31.csv.gz',
# 'wtn.16.32.csv.gz',
# 'wtn.16.33.csv.gz',
# 'wtn.16.34.csv.gz',
# 'wtn.16.35.csv.gz',
# 'wtn.16.36.csv.gz',
# 'wtn.16.37.csv.gz',
# 'wtn.16.38.csv.gz',
# 'wtn.16.39.csv.gz',
# 'wtn.16.4.csv.gz',
# 'wtn.16.40.csv.gz',
# 'wtn.16.41.csv.gz',
# 'wtn.16.42.csv.gz',
# 'wtn.16.43.csv.gz',
# 'wtn.16.44.csv.gz',
# 'wtn.16.45.csv.gz',
# 'wtn.16.46.csv.gz',
# 'wtn.16.47.csv.gz',
# 'wtn.16.48.csv.gz',
# 'wtn.16.49.csv.gz',
# 'wtn.16.5.csv.gz',
# 'wtn.16.50.csv.gz',
# 'wtn.16.51.csv.gz',
# 'wtn.16.52.csv.gz',
# 'wtn.16.53.csv.gz',
# 'wtn.16.54.csv.gz',
# 'wtn.16.55.csv.gz',
# 'wtn.16.6.csv.gz',
# 'wtn.16.7.csv.gz',
# 'wtn.16.8.csv.gz',
# 'wtn.16.9.csv.gz',
# 'wtn.17.1.csv.gz',
# 'wtn.17.10.csv.gz',
# 'wtn.17.11.csv.gz',
# 'wtn.17.12.csv.gz',
# 'wtn.17.13.csv.gz',
# 'wtn.17.14.csv.gz',
# 'wtn.17.15.csv.gz',
# 'wtn.17.16.csv.gz',
# 'wtn.17.17.csv.gz',
# 'wtn.17.18.csv.gz',
# 'wtn.17.19.csv.gz',
# 'wtn.17.2.csv.gz',
# 'wtn.17.20.csv.gz',
# 'wtn.17.21.csv.gz',
# 'wtn.17.22.csv.gz',
# 'wtn.17.23.csv.gz',
# 'wtn.17.24.csv.gz',
# 'wtn.17.25.csv.gz',
# 'wtn.17.26.csv.gz',
# 'wtn.17.27.csv.gz',
# 'wtn.17.28.csv.gz',
# 'wtn.17.29.csv.gz',
# 'wtn.17.3.csv.gz',
# 'wtn.17.30.csv.gz',
# 'wtn.17.31.csv.gz',
# 'wtn.17.32.csv.gz',
# 'wtn.17.33.csv.gz',
# 'wtn.17.34.csv.gz',
# 'wtn.17.35.csv.gz',
# 'wtn.17.36.csv.gz',
# 'wtn.17.37.csv.gz',
# 'wtn.17.38.csv.gz',
# 'wtn.17.39.csv.gz',
# 'wtn.17.4.csv.gz',
# 'wtn.17.40.csv.gz',
# 'wtn.17.41.csv.gz',
# 'wtn.17.42.csv.gz',
# 'wtn.17.43.csv.gz',
# 'wtn.17.44.csv.gz',
# 'wtn.17.45.csv.gz',
# 'wtn.17.46.csv.gz',
# 'wtn.17.47.csv.gz',
# 'wtn.17.48.csv.gz',
# 'wtn.17.49.csv.gz',
# 'wtn.17.5.csv.gz',
# 'wtn.17.50.csv.gz',
# 'wtn.17.51.csv.gz',
# 'wtn.17.52.csv.gz',
# 'wtn.17.53.csv.gz',
# 'wtn.17.54.csv.gz',
# 'wtn.17.55.csv.gz',
# 'wtn.17.6.csv.gz',
# 'wtn.17.7.csv.gz',
# 'wtn.17.8.csv.gz',
# 'wtn.17.9.csv.gz',
# 'wtn.18.1.csv.gz',
# 'wtn.18.10.csv.gz',
# 'wtn.18.11.csv.gz',
# 'wtn.18.12.csv.gz',
# 'wtn.18.13.csv.gz',
# 'wtn.18.14.csv.gz',
# 'wtn.18.15.csv.gz',
# 'wtn.18.16.csv.gz',
# 'wtn.18.17.csv.gz',
# 'wtn.18.18.csv.gz',#Resetting 8 works fine, and the duplicates are OK; S12 is nearly flat, which means there are a lot of duplicates, but they are probably OK; made adjustments on S14 because the clock flag is set and the clock was adding a new record every time, including some repeated frames; FIXED
# 'wtn.18.19.csv.gz',
# 'wtn.18.2.csv.gz',
# 'wtn.18.20.csv.gz',
# 'wtn.18.21.csv.gz',
# 'wtn.18.22.csv.gz',
# 'wtn.18.23.csv.gz',
# 'wtn.18.24.csv.gz',
# 'wtn.18.25.csv.gz',
# 'wtn.18.26.csv.gz',
# 'wtn.18.27.csv.gz',
# 'wtn.18.28.csv.gz',
# 'wtn.18.29.csv.gz',
# 'wtn.18.3.csv.gz',
# 'wtn.18.30.csv.gz',
# 'wtn.18.31.csv.gz',
# 'wtn.18.32.csv.gz',
# 'wtn.18.33.csv.gz',
# 'wtn.18.34.csv.gz',
# 'wtn.18.35.csv.gz',
# 'wtn.18.36.csv.gz',
# 'wtn.18.37.csv.gz',
# 'wtn.18.38.csv.gz', #Data Check - 180 clock flag(s)
# 'wtn.18.39.csv.gz',
# 'wtn.18.4.csv.gz',
# 'wtn.18.40.csv.gz',
# 'wtn.18.41.csv.gz',
# 'wtn.18.42.csv.gz',
# 'wtn.18.43.csv.gz',
# 'wtn.18.44.csv.gz',
# 'wtn.18.45.csv.gz',
# 'wtn.18.46.csv.gz',
# 'wtn.18.47.csv.gz',
# 'wtn.18.48.csv.gz',
# 'wtn.18.49.csv.gz',
# 'wtn.18.5.csv.gz',
# 'wtn.18.50.csv.gz',
# 'wtn.18.51.csv.gz',
# 'wtn.18.52.csv.gz',
# 'wtn.18.53.csv.gz',
# 'wtn.18.54.csv.gz',
# 'wtn.18.55.csv.gz',
# 'wtn.18.6.csv.gz',
# 'wtn.18.7.csv.gz',
# 'wtn.18.8.csv.gz',
# 'wtn.18.9.csv.gz',
# 'wtn.19.1.csv.gz',
# 'wtn.19.10.csv.gz',
# 'wtn.19.11.csv.gz',
# 'wtn.19.12.csv.gz',
# 'wtn.19.13.csv.gz',
# 'wtn.19.14.csv.gz',
# 'wtn.19.15.csv.gz',
# 'wtn.19.16.csv.gz',
# 'wtn.19.17.csv.gz',
# 'wtn.19.18.csv.gz',
# 'wtn.19.19.csv.gz',
# 'wtn.19.2.csv.gz',
# 'wtn.19.20.csv.gz',
# 'wtn.19.21.csv.gz',
# 'wtn.19.22.csv.gz',
# 'wtn.19.23.csv.gz',
# 'wtn.19.24.csv.gz',
# 'wtn.19.25.csv.gz',
# 'wtn.19.26.csv.gz',
# 'wtn.19.27.csv.gz',
# 'wtn.19.28.csv.gz',
# 'wtn.19.29.csv.gz',
# 'wtn.19.3.csv.gz',
# 'wtn.19.30.csv.gz',
# 'wtn.19.31.csv.gz',
# 'wtn.19.32.csv.gz',
# 'wtn.19.33.csv.gz',
# 'wtn.19.34.csv.gz',
# 'wtn.19.35.csv.gz',
# 'wtn.19.36.csv.gz',
# 'wtn.19.37.csv.gz',
# 'wtn.19.38.csv.gz',
# 'wtn.19.39.csv.gz',
# 'wtn.19.4.csv.gz',
# 'wtn.19.40.csv.gz',
# 'wtn.19.41.csv.gz',
# 'wtn.19.42.csv.gz',
# 'wtn.19.43.csv.gz',
# 'wtn.19.44.csv.gz',
# 'wtn.19.45.csv.gz',
# 'wtn.19.46.csv.gz',
# 'wtn.19.47.csv.gz',
# 'wtn.19.48.csv.gz',
# 'wtn.19.49.csv.gz',
# 'wtn.19.5.csv.gz',
# 'wtn.19.50.csv.gz',
# 'wtn.19.51.csv.gz',
# 'wtn.19.52.csv.gz',
# 'wtn.19.53.csv.gz',
# 'wtn.19.54.csv.gz',
# 'wtn.19.55.csv.gz',
# 'wtn.19.6.csv.gz',
# 'wtn.19.7.csv.gz',
# 'wtn.19.8.csv.gz',
# 'wtn.19.9.csv.gz',
# 'wtn.2.1.csv.gz',
# 'wtn.2.10.csv.gz',
# 'wtn.2.11.csv.gz',
# 'wtn.2.12.csv.gz',
# 'wtn.2.13.csv.gz',
# 'wtn.2.14.csv.gz',
# 'wtn.2.15.csv.gz',
# 'wtn.2.16.csv.gz',
# 'wtn.2.17.csv.gz',
# 'wtn.2.18.csv.gz',
# 'wtn.2.19.csv.gz',
# 'wtn.2.2.csv.gz',
# 'wtn.2.20.csv.gz',
# 'wtn.2.21.csv.gz',
# 'wtn.2.22.csv.gz',
# 'wtn.2.23.csv.gz',
# 'wtn.2.24.csv.gz',
# 'wtn.2.25.csv.gz',
# 'wtn.2.26.csv.gz',
# 'wtn.2.27.csv.gz',
# 'wtn.2.28.csv.gz',
# 'wtn.2.29.csv.gz',
# 'wtn.2.3.csv.gz',
# 'wtn.2.30.csv.gz',
# 'wtn.2.31.csv.gz',
# 'wtn.2.32.csv.gz',
# 'wtn.2.33.csv.gz',
# 'wtn.2.34.csv.gz',
# 'wtn.2.35.csv.gz',
# 'wtn.2.36.csv.gz',
# 'wtn.2.37.csv.gz',
# 'wtn.2.38.csv.gz',
# 'wtn.2.39.csv.gz',
# 'wtn.2.4.csv.gz',
# 'wtn.2.40.csv.gz',
# 'wtn.2.41.csv.gz',
# 'wtn.2.42.csv.gz',
# 'wtn.2.43.csv.gz',
# 'wtn.2.44.csv.gz',
# 'wtn.2.45.csv.gz',
# 'wtn.2.46.csv.gz',
# 'wtn.2.47.csv.gz',
# 'wtn.2.48.csv.gz',
# 'wtn.2.49.csv.gz',
# 'wtn.2.5.csv.gz',
# 'wtn.2.50.csv.gz',
# 'wtn.2.51.csv.gz',
# 'wtn.2.52.csv.gz',
# 'wtn.2.53.csv.gz',
# 'wtn.2.54.csv.gz',
# 'wtn.2.55.csv.gz',
# 'wtn.2.6.csv.gz',
# 'wtn.2.7.csv.gz',
# 'wtn.2.8.csv.gz',
# 'wtn.2.9.csv.gz',
# 'wtn.20.1.csv.gz',
# 'wtn.20.10.csv.gz',
# 'wtn.20.11.csv.gz',
# 'wtn.20.12.csv.gz',
# 'wtn.20.13.csv.gz',
# 'wtn.20.14.csv.gz',
# 'wtn.20.15.csv.gz',
# 'wtn.20.16.csv.gz',
# 'wtn.20.17.csv.gz',
# 'wtn.20.18.csv.gz',
# 'wtn.20.19.csv.gz',
# 'wtn.20.2.csv.gz',
# 'wtn.20.20.csv.gz',
# 'wtn.20.21.csv.gz',
# 'wtn.20.22.csv.gz',
# 'wtn.20.23.csv.gz',
# 'wtn.20.24.csv.gz',
# 'wtn.20.25.csv.gz',
# 'wtn.20.26.csv.gz',
# 'wtn.20.27.csv.gz',
# 'wtn.20.28.csv.gz',
# 'wtn.20.29.csv.gz',
# 'wtn.20.3.csv.gz',
# 'wtn.20.30.csv.gz',
# 'wtn.20.31.csv.gz',
# 'wtn.20.32.csv.gz',
# 'wtn.20.33.csv.gz',
# 'wtn.20.34.csv.gz',
# 'wtn.20.35.csv.gz',
# 'wtn.20.36.csv.gz',
# 'wtn.20.37.csv.gz',
# 'wtn.20.38.csv.gz',
# 'wtn.20.39.csv.gz',
# 'wtn.20.4.csv.gz',
# 'wtn.20.40.csv.gz',
# 'wtn.20.41.csv.gz',
# 'wtn.20.42.csv.gz',
# 'wtn.20.43.csv.gz',
# 'wtn.20.44.csv.gz',
# 'wtn.20.45.csv.gz',
# 'wtn.20.46.csv.gz',
# 'wtn.20.47.csv.gz',
# 'wtn.20.48.csv.gz',
# 'wtn.20.49.csv.gz',
# 'wtn.20.5.csv.gz',
# 'wtn.20.50.csv.gz',
# 'wtn.20.51.csv.gz',
# 'wtn.20.52.csv.gz',
# 'wtn.20.53.csv.gz',
# 'wtn.20.54.csv.gz',
# 'wtn.20.6.csv.gz',
# 'wtn.20.7.csv.gz',
# 'wtn.20.8.csv.gz',
# 'wtn.20.9.csv.gz',
# 'wtn.21.1.csv.gz',
# 'wtn.21.10.csv.gz',
# 'wtn.21.11.csv.gz',
# 'wtn.21.12.csv.gz',
# 'wtn.21.13.csv.gz',
# 'wtn.21.14.csv.gz',
# 'wtn.21.15.csv.gz',
# 'wtn.21.16.csv.gz',
# 'wtn.21.17.csv.gz',
# 'wtn.21.18.csv.gz',
# 'wtn.21.19.csv.gz',
# 'wtn.21.2.csv.gz',
# 'wtn.21.20.csv.gz',
# 'wtn.21.21.csv.gz',
# 'wtn.21.22.csv.gz',
# 'wtn.21.23.csv.gz',
# 'wtn.21.24.csv.gz',
# 'wtn.21.25.csv.gz',
# 'wtn.21.26.csv.gz',
# 'wtn.21.27.csv.gz',
# 'wtn.21.28.csv.gz',
# 'wtn.21.29.csv.gz',
# 'wtn.21.3.csv.gz',
# 'wtn.21.30.csv.gz',
# 'wtn.21.31.csv.gz',
# 'wtn.21.32.csv.gz',
# 'wtn.21.33.csv.gz',
# 'wtn.21.34.csv.gz',
# 'wtn.21.35.csv.gz',
# 'wtn.21.36.csv.gz',
# 'wtn.21.37.csv.gz',
# 'wtn.21.38.csv.gz',
# 'wtn.21.39.csv.gz',
# 'wtn.21.4.csv.gz',
# 'wtn.21.40.csv.gz',
# 'wtn.21.41.csv.gz',
# 'wtn.21.42.csv.gz',
# 'wtn.21.43.csv.gz',
# 'wtn.21.44.csv.gz',
# 'wtn.21.45.csv.gz',
# 'wtn.21.46.csv.gz',
# 'wtn.21.47.csv.gz',
# 'wtn.21.48.csv.gz',
# 'wtn.21.5.csv.gz',
# 'wtn.21.6.csv.gz',
# 'wtn.21.7.csv.gz',
# 'wtn.21.8.csv.gz',
# 'wtn.21.9.csv.gz',
# 'wtn.3.1.csv.gz',
# 'wtn.3.10.csv.gz',

# 'wtn.3.11.csv.gz',
# 'wtn.3.12.csv.gz',
# 'wtn.3.13.csv.gz',
# 'wtn.3.14.csv.gz',
# 'wtn.3.15.csv.gz',
# 'wtn.3.16.csv.gz',
# 'wtn.3.17.csv.gz',
# 'wtn.3.18.csv.gz',
# 'wtn.3.19.csv.gz',
# 'wtn.3.2.csv.gz',
# 'wtn.3.20.csv.gz',
# 'wtn.3.21.csv.gz',
# 'wtn.3.22.csv.gz',
# 'wtn.3.23.csv.gz',
# 'wtn.3.24.csv.gz',
# 'wtn.3.25.csv.gz',
# 'wtn.3.26.csv.gz',
# 'wtn.3.27.csv.gz',
# 'wtn.3.28.csv.gz',
# 'wtn.3.29.csv.gz',
# 'wtn.3.3.csv.gz',
# 'wtn.3.30.csv.gz',
# 'wtn.3.31.csv.gz',
# 'wtn.3.32.csv.gz',
# 'wtn.3.33.csv.gz',
# 'wtn.3.34.csv.gz',
# 'wtn.3.35.csv.gz',
# 'wtn.3.36.csv.gz',
# 'wtn.3.37.csv.gz',
# 'wtn.3.38.csv.gz',
# 'wtn.3.39.csv.gz',
# 'wtn.3.4.csv.gz',
# 'wtn.3.40.csv.gz',
# 'wtn.3.41.csv.gz',
# 'wtn.3.42.csv.gz',
# 'wtn.3.43.csv.gz',
# 'wtn.3.44.csv.gz',
# 'wtn.3.45.csv.gz',
# 'wtn.3.46.csv.gz',
# 'wtn.3.47.csv.gz',
# 'wtn.3.48.csv.gz',
# 'wtn.3.49.csv.gz',
# 'wtn.3.5.csv.gz',
# 'wtn.3.50.csv.gz',
# 'wtn.3.51.csv.gz',
# 'wtn.3.52.csv.gz',
# 'wtn.3.53.csv.gz',
# 'wtn.3.54.csv.gz',
# 'wtn.3.55.csv.gz',
# 'wtn.3.6.csv.gz',
# 'wtn.3.7.csv.gz',
# 'wtn.3.8.csv.gz',
# 'wtn.3.9.csv.gz',
# 'wtn.4.1.csv.gz',
# 'wtn.4.10.csv.gz',
# 'wtn.4.11.csv.gz',
# 'wtn.4.12.csv.gz',
# 'wtn.4.13.csv.gz',
# 'wtn.4.14.csv.gz',
# 'wtn.4.15.csv.gz',
# 'wtn.4.16.csv.gz',
# 'wtn.4.17.csv.gz',
# 'wtn.4.18.csv.gz',
# 'wtn.4.19.csv.gz',
# 'wtn.4.2.csv.gz',
# 'wtn.4.20.csv.gz',
# 'wtn.4.21.csv.gz',
# 'wtn.4.22.csv.gz',
# 'wtn.4.23.csv.gz',
# 'wtn.4.24.csv.gz',
# 'wtn.4.25.csv.gz',
# 'wtn.4.26.csv.gz',
# 'wtn.4.27.csv.gz',
# 'wtn.4.28.csv.gz',
# 'wtn.4.29.csv.gz',
# 'wtn.4.3.csv.gz',
# 'wtn.4.30.csv.gz',
# 'wtn.4.31.csv.gz',
# 'wtn.4.32.csv.gz',
# 'wtn.4.33.csv.gz',
# 'wtn.4.34.csv.gz',
# 'wtn.4.35.csv.gz',
# 'wtn.4.36.csv.gz',
# 'wtn.4.37.csv.gz',
# 'wtn.4.38.csv.gz',
# 'wtn.4.39.csv.gz',
# 'wtn.4.4.csv.gz',
# 'wtn.4.40.csv.gz',
# 'wtn.4.41.csv.gz',
# 'wtn.4.42.csv.gz',
# 'wtn.4.43.csv.gz',
# 'wtn.4.44.csv.gz',
# 'wtn.4.45.csv.gz',
# 'wtn.4.46.csv.gz',
# 'wtn.4.47.csv.gz',
# 'wtn.4.48.csv.gz',
# 'wtn.4.49.csv.gz',
# 'wtn.4.5.csv.gz',
# 'wtn.4.50.csv.gz',
# 'wtn.4.51.csv.gz',
# 'wtn.4.52.csv.gz',
# 'wtn.4.53.csv.gz',
# 'wtn.4.54.csv.gz',
# 'wtn.4.55.csv.gz',
# 'wtn.4.6.csv.gz',
# 'wtn.4.7.csv.gz',
# 'wtn.4.8.csv.gz',
# 'wtn.4.9.csv.gz',
# 'wtn.5.1.csv.gz',
# 'wtn.5.10.csv.gz',
# 'wtn.5.11.csv.gz',
# 'wtn.5.12.csv.gz',
# 'wtn.5.13.csv.gz',
# 'wtn.5.14.csv.gz',
# 'wtn.5.15.csv.gz',
# 'wtn.5.16.csv.gz',
# 'wtn.5.17.csv.gz',
# 'wtn.5.18.csv.gz',
# 'wtn.5.19.csv.gz',
# 'wtn.5.2.csv.gz',
# 'wtn.5.20.csv.gz',
# 'wtn.5.21.csv.gz',
# 'wtn.5.22.csv.gz',
# 'wtn.5.23.csv.gz',
# 'wtn.5.24.csv.gz',
# 'wtn.5.25.csv.gz',
# 'wtn.5.26.csv.gz',
# 'wtn.5.27.csv.gz',
# 'wtn.5.28.csv.gz',
# 'wtn.5.29.csv.gz',
# 'wtn.5.3.csv.gz',
# 'wtn.5.30.csv.gz',
# 'wtn.5.31.csv.gz',
# 'wtn.5.32.csv.gz',
# 'wtn.5.33.csv.gz',
# 'wtn.5.34.csv.gz',
# 'wtn.5.35.csv.gz',
# 'wtn.5.36.csv.gz',
# 'wtn.5.37.csv.gz',
# 'wtn.5.38.csv.gz',
# 'wtn.5.39.csv.gz',
# 'wtn.5.4.csv.gz',
# 'wtn.5.40.csv.gz',
# 'wtn.5.41.csv.gz',
# 'wtn.5.42.csv.gz',
# 'wtn.5.43.csv.gz',
# 'wtn.5.44.csv.gz',
# 'wtn.5.45.csv.gz',
# 'wtn.5.46.csv.gz',
# 'wtn.5.47.csv.gz',
# 'wtn.5.48.csv.gz',
# 'wtn.5.49.csv.gz',
# 'wtn.5.5.csv.gz',
# 'wtn.5.50.csv.gz',
# 'wtn.5.51.csv.gz',
# 'wtn.5.52.csv.gz',
# 'wtn.5.53.csv.gz',
# 'wtn.5.54.csv.gz',
# 'wtn.5.55.csv.gz',
# 'wtn.5.6.csv.gz',
# 'wtn.5.7.csv.gz',
# 'wtn.5.8.csv.gz',
# 'wtn.5.9.csv.gz',
# 'wtn.6.1.csv.gz',
# 'wtn.6.10.csv.gz',
# 'wtn.6.11.csv.gz',
# 'wtn.6.12.csv.gz',
# 'wtn.6.13.csv.gz',
# 'wtn.6.14.csv.gz',
# 'wtn.6.15.csv.gz',
# 'wtn.6.16.csv.gz',
# 'wtn.6.17.csv.gz',
# 'wtn.6.18.csv.gz',
# 'wtn.6.19.csv.gz',
# 'wtn.6.2.csv.gz',
# 'wtn.6.20.csv.gz',
# 'wtn.6.21.csv.gz',
# 'wtn.6.22.csv.gz',
# 'wtn.6.23.csv.gz',
# 'wtn.6.24.csv.gz',
# 'wtn.6.25.csv.gz',
# 'wtn.6.26.csv.gz',
# 'wtn.6.27.csv.gz',
# 'wtn.6.28.csv.gz',
# 'wtn.6.29.csv.gz',
# 'wtn.6.3.csv.gz',
# 'wtn.6.30.csv.gz',
# 'wtn.6.31.csv.gz',
# 'wtn.6.32.csv.gz',
# 'wtn.6.33.csv.gz',
# 'wtn.6.34.csv.gz',
# 'wtn.6.35.csv.gz',
# 'wtn.6.36.csv.gz',
# 'wtn.6.37.csv.gz',
# 'wtn.6.38.csv.gz',
# 'wtn.6.39.csv.gz',
# 'wtn.6.4.csv.gz',
# 'wtn.6.40.csv.gz',
# 'wtn.6.41.csv.gz',
# 'wtn.6.42.csv.gz',
# 'wtn.6.43.csv.gz',
# 'wtn.6.44.csv.gz',
# 'wtn.6.45.csv.gz',
# 'wtn.6.46.csv.gz',
# 'wtn.6.47.csv.gz',
# 'wtn.6.48.csv.gz',
# 'wtn.6.49.csv.gz',
# 'wtn.6.5.csv.gz',
# 'wtn.6.50.csv.gz',
# 'wtn.6.51.csv.gz',
# 'wtn.6.52.csv.gz',
# 'wtn.6.53.csv.gz',
# 'wtn.6.54.csv.gz',
# 'wtn.6.55.csv.gz',
# 'wtn.6.6.csv.gz',
# 'wtn.6.7.csv.gz',
# 'wtn.6.8.csv.gz',
# 'wtn.6.9.csv.gz',
# 'wtn.7.1.csv.gz',
# 'wtn.7.10.csv.gz',
# 'wtn.7.11.csv.gz',
# 'wtn.7.12.csv.gz',
# 'wtn.7.13.csv.gz',
# 'wtn.7.14.csv.gz',
# 'wtn.7.15.csv.gz',
# 'wtn.7.16.csv.gz',
# 'wtn.7.17.csv.gz',
# 'wtn.7.18.csv.gz',
# 'wtn.7.19.csv.gz',
# 'wtn.7.2.csv.gz',
# 'wtn.7.20.csv.gz',
# 'wtn.7.21.csv.gz',
# 'wtn.7.22.csv.gz',
# 'wtn.7.23.csv.gz',
# 'wtn.7.24.csv.gz',
# 'wtn.7.25.csv.gz',
# 'wtn.7.26.csv.gz',
# 'wtn.7.27.csv.gz',
# 'wtn.7.28.csv.gz',
# 'wtn.7.29.csv.gz',
# 'wtn.7.3.csv.gz',
# 'wtn.7.30.csv.gz',
# 'wtn.7.31.csv.gz',
# 'wtn.7.32.csv.gz',
# 'wtn.7.33.csv.gz',
# 'wtn.7.34.csv.gz',
# 'wtn.7.35.csv.gz',
# 'wtn.7.36.csv.gz',
# 'wtn.7.37.csv.gz',
# 'wtn.7.38.csv.gz',
# 'wtn.7.39.csv.gz',
# 'wtn.7.4.csv.gz',
# 'wtn.7.40.csv.gz',
# 'wtn.7.41.csv.gz',
# 'wtn.7.42.csv.gz',
# 'wtn.7.43.csv.gz',
# 'wtn.7.44.csv.gz',
# 'wtn.7.45.csv.gz',
# 'wtn.7.46.csv.gz',
# 'wtn.7.47.csv.gz',
# 'wtn.7.48.csv.gz',
# 'wtn.7.49.csv.gz',
# 'wtn.7.5.csv.gz',
# 'wtn.7.50.csv.gz',
# 'wtn.7.51.csv.gz',
# 'wtn.7.52.csv.gz',
# 'wtn.7.53.csv.gz',
# 'wtn.7.54.csv.gz',
# 'wtn.7.55.csv.gz',
# 'wtn.7.56.csv.gz',
# 'wtn.7.6.csv.gz',
# 'wtn.7.7.csv.gz',
# 'wtn.7.8.csv.gz',
# 'wtn.7.9.csv.gz',
# 'wtn.8.1.csv.gz',
# 'wtn.8.10.csv.gz',
# 'wtn.8.11.csv.gz',
# 'wtn.8.12.csv.gz',
# 'wtn.8.13.csv.gz',
# 'wtn.8.14.csv.gz',
# 'wtn.8.15.csv.gz',
# 'wtn.8.16.csv.gz',
# 'wtn.8.17.csv.gz',
# 'wtn.8.18.csv.gz',
# 'wtn.8.19.csv.gz',
# 'wtn.8.2.csv.gz',
# 'wtn.8.20.csv.gz',
# 'wtn.8.21.csv.gz',
# 'wtn.8.22.csv.gz',
# 'wtn.8.23.csv.gz',
# 'wtn.8.24.csv.gz',
# 'wtn.8.25.csv.gz',
# 'wtn.8.26.csv.gz',
# 'wtn.8.27.csv.gz',
# 'wtn.8.28.csv.gz',
# 'wtn.8.29.csv.gz',
# 'wtn.8.3.csv.gz',
# 'wtn.8.30.csv.gz',
# 'wtn.8.31.csv.gz',
# 'wtn.8.32.csv.gz',
# 'wtn.8.33.csv.gz',
# 'wtn.8.34.csv.gz',
# 'wtn.8.35.csv.gz',
# 'wtn.8.36.csv.gz',
# 'wtn.8.37.csv.gz',
# 'wtn.8.38.csv.gz',
# 'wtn.8.39.csv.gz',
# 'wtn.8.4.csv.gz',
# 'wtn.8.40.csv.gz',
# 'wtn.8.41.csv.gz',
# 
# 'wtn.8.42.csv.gz',
# 'wtn.8.43.csv.gz',
# 'wtn.8.44.csv.gz',
# 'wtn.8.45.csv.gz',
# 'wtn.8.46.csv.gz',
# 'wtn.8.47.csv.gz',
# 'wtn.8.48.csv.gz',
# 'wtn.8.49.csv.gz',
# 'wtn.8.5.csv.gz',
# 'wtn.8.50.csv.gz',
# 'wtn.8.51.csv.gz',
# 'wtn.8.52.csv.gz',
# 'wtn.8.53.csv.gz',
# 'wtn.8.54.csv.gz',#180 clock flag records
# 'wtn.8.55.csv.gz',
# 'wtn.8.6.csv.gz',
# 'wtn.8.7.csv.gz',
# 'wtn.8.8.csv.gz',
# 'wtn.8.9.csv.gz',
# 'wtn.9.1.csv.gz',
# 'wtn.9.10.csv.gz',
# 'wtn.9.11.csv.gz',
# 'wtn.9.12.csv.gz',
# 'wtn.9.13.csv.gz',
# 'wtn.9.14.csv.gz',
# 'wtn.9.15.csv.gz',
# 'wtn.9.16.csv.gz',
# 'wtn.9.17.csv.gz',
# 'wtn.9.18.csv.gz',
# 'wtn.9.19.csv.gz',
# 'wtn.9.2.csv.gz',
# 'wtn.9.20.csv.gz',
# 'wtn.9.21.csv.gz',
# 'wtn.9.22.csv.gz',
# 'wtn.9.23.csv.gz',
# 'wtn.9.24.csv.gz',
# 'wtn.9.25.csv.gz',
# 'wtn.9.26.csv.gz',
# 'wtn.9.27.csv.gz',
# 'wtn.9.28.csv.gz',
# 'wtn.9.29.csv.gz',
# 'wtn.9.3.csv.gz',
# 'wtn.9.30.csv.gz',
# 'wtn.9.31.csv.gz',
# 'wtn.9.32.csv.gz',
# 'wtn.9.33.csv.gz',
# 'wtn.9.34.csv.gz',
# 'wtn.9.35.csv.gz',
# 'wtn.9.36.csv.gz',
# 'wtn.9.37.csv.gz',
# 'wtn.9.38.csv.gz',
# 'wtn.9.39.csv.gz',
# 'wtn.9.4.csv.gz',
# 'wtn.9.40.csv.gz',
# 'wtn.9.41.csv.gz',
# 'wtn.9.42.csv.gz',
# 'wtn.9.43.csv.gz',
# 'wtn.9.44.csv.gz',
# 'wtn.9.45.csv.gz',
# 'wtn.9.46.csv.gz',
# 'wtn.9.47.csv.gz',
# 'wtn.9.48.csv.gz',
# 'wtn.9.49.csv.gz',
# 'wtn.9.5.csv.gz',
# 'wtn.9.50.csv.gz',
# 'wtn.9.51.csv.gz',
# 'wtn.9.52.csv.gz',
# 'wtn.9.53.csv.gz',
# 'wtn.9.54.csv.gz',
# 'wtn.9.55.csv.gz',
# 'wtn.9.6.csv.gz',
# 'wtn.9.7.csv.gz',
# 'wtn.9.8.csv.gz',
# 'wtn.9.9.csv.gz'
]


    # this has -7777 in pairs 
    # filenames=['wtn.10.1.csv.gz',]
    # # call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,
    # #     filenames=filenames,single_station='S15',single_ground_station=9)



    # error with Fatal final check 
    # filenames=['wtn.20.1.csv.gz',]

# 'wtn.10.49.csv.gz',
# this is an excellent example where the clock flag is turning on and off
# and it doesn't end up deleting the -7777 


    # an example with even more nasty clock flags
    # filenames=['wtn.10.31.csv.gz',]

    # # an example with lots of nasty clock flags 
    # filenames=['wtn.7.49.csv.gz',]

# Gap: 4884.29, Corr_gap -8888
# df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.10.29.csv.gz', 
#               start_timestamp ='1976-12-05 14:30:01.227000+00:00', 
#               end_timestamp='',
#               start_station='S16',
#               end_station='',
#               corr_ground_station1=9)

    # filenames=['wtn.10.29.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,
    #   filenames=filenames,logging_level=logging.DEBUG,initial=False)
# single_station='S16',single_ground_station=9)
# 242,"1976-12-05T22:31:59.263000Z",40,"S14",0,1,497
# single_station='S16',single_ground_station=9)

# error 2 1237
# WARNING: Ground station/Station - 9 S15 Amending record backwards from valid_idx=18111 begin=0 error=1237, end=1238 len=19338 before_begin []
# # 9 S16 
# ,initial=False,single_station='S15',single_ground_station=9)
# cumsum_test=50)

# WARNING: Ground station/Station - 4 S15 Amending record forwards from valid_idx=2847 begin=0 error=1 end=4 len=1226
# cumsum_test=10,
      # 

# S15          0                    4  


# WARNING: Ground station/Station - 1 S12 Amending record backwards from valid_idx=55880 begin=142 end=37712 len=56011 before_begin 395
# DEBUG: station_fix_timeskips: checking backwards
# WARNING: Ground station/Station - 1 S12 Amending record backwards from valid_idx=55880 begin=0 end=37713 len=56011 before_begin 142
# 





    # run all (main version)
    call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=checked_dir_full,
      filenames=filenames_all,logging_level=logging.INFO)
    # 
    # exit()

# WARNING: Ground station/Station - 1 S12 Amending record forwards from valid_idx=2709 begin=2959 end=3481 len=2744

# 4 Station: S15
    # filenames=['wtn.1.1.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,
    #     filenames=filenames,logging_level=logging.DEBUG)





    # # logger = logging.getLogger('simple_example')
    # logger = logging.getLogger()
    # handlers = logger.handlers[:]
    # for handler in handlers:
    #     handler.close()
    #     logger.removeHandler(handler)
    # 
    # logging.basicConfig(filename='c.log', filemode='w', level=logging.INFO)
    # logger = logging.getLogger()
    # print(type(logger))
    # logging.info('info message c')
    # 
    # # close the handler
    # handlers = logger.handlers[:]
    # for handler in handlers:
    #     handler.close()
    #     logger.removeHandler(handler)
    # 
    # # note that the formatting message can have no extras in the log file 
    # logging.basicConfig(filename='d.log', filemode='w', level=logging.INFO, format='%(message)s')
    # logger = logging.getLogger()
    # print(type(logger))
    # logging.info('info message d')
    # 
    # # close the handler
    # handlers = logger.handlers[:]
    # for handler in handlers:
    #     handler.close()
    #     logger.removeHandler(handler)
    # 
    # exit()

    # # logger.setLevel(logging.INFO)
    # # create file handler which logs even debug messages
    # # fh = logging.FileHandler('a.log')
    # # fh.setLevel(logging.INFO)
    # 
    # a = logging.getLogger()
    # pr
    # 
    # # print(a)
    # 
    # # logger.addHandler(fh)
    # 
    # logging.info('info message a')

    # a.removeHandler()

    # ch = logging.FileHandler('b.log')
    # ch.setLevel(logging.INFO)
    # logger.addHandler(ch)
    # logger.info('info message b')


    # logger = logging.getLogger('simple_example')

    # logging.basicConfig(filename='a.log', filemode='w', level=logging.INFO)
    # logger = logging.getLogger()
    # logging.info('info message a')  
    # logger.close()
    # logging.shutdown()
    # # 
    # 
    # logging.basicConfig(filename='b.log', filemode='w', level=logging.INFO)
    # logging.info('info message b')


    # exit()
    # 
    # logger = logging.getLogger('simple_example')
    # 
    # # logging.basicConfig(filename=log_filename, filemode='w', level=logging_level)
    # logger.setLevel(logging.INFO)
    # # create file handler which logs even debug messages
    # fh = logging.FileHandler('a.log')
    # fh.setLevel(logging.INFO)
    # logger.addHandler(fh)
    # 
    # logger.info('info message a')
    # 
    # logger.removeHandler(fh)
    # 
    # ch = logging.FileHandler('b.log')
    # ch.setLevel(logging.INFO)
    # logger.addHandler(ch)
    # 
    # logger.info('info message b')
    # 
    # 
    # exit()


    # log_filename='logs/import_first_entry_S12_XXXXX.log'
    # 
    # # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    # #   filenames=['wtn.21.6.csv.gz'],log_filename=log_filename)
    # 
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=['wtn.20.6.csv.gz'],log_filename=log_filename)

    # filenames=['wtn.1.1.error.csv.gz']

    # I think this is working: 
    # filenames=['wtn.1.1.all.csv.gz']

    # I think this is working: 
    # filenames=['wtn.20.5.csv.gz']
    
    # filenames=['wtn.20.6.all.csv.gz']
    # filenames=['wtn.1.5.all.csv.gz']
    # filenames=['wtn.20.6.small.csv.gz']

    
    # filenames=['wtn.1.1.all.csv.gz', 'wtn.1.5.all.csv.gz',]
    # 
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames,single_station='S12',single_ground_station=1)

    # check 3 or 4 
    # filenames=['wtn.9.27.gz',]

    # filenames=['wtn.14.44.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames)

    # filenames=['wtn.6.30.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames)

    # filenames=['wtn.6.30.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames)

    # filenames=['wtn.11.22.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames)


    # filenames=['wtn.3.44.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames,single_ground_station=9,single_station='S14')

    # filenames=['wtn.7.10.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames,single_ground_station=9,single_station='S15')

    # Fix example with clock flag and exact gaps 
    # filenames=['wtn.4.31.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames,single_ground_station=4,single_station='S16')





    # filenames=['wtn.18.40.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames,single_ground_station=7,single_station='S12')



    # filenames=['wtn.11.7.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames)


    # filenames=['wtn.1.1.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames)
    
    # exit()

    # filenames=['wtn.7.39.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames)


    # filenames=['wtn.11.22.csv.gz',]
    # # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    # #   filenames=filenames)
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames,single_station='S17',single_ground_station=3)

    # An example of a whole section where the clock flag is set 
    # filenames=['wtn.8.54.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames,single_ground_station=7,single_station='S12')

    # An example of clock flags set (but this example is just wrong - and outside the range)
    # filenames=['wtn.9.40.csv.gz',]
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=filenames)



    # # run all of them
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=None)



    # 
    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=['wtna.20.6.csv.gz'],log_filename=log_filename)

    # call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,
    #   filenames=['wtn.20.6.csv.gz'],log_filename=log_filename)





# def run_test():
#     checked_dir='/Users/nunn/lunar_data/PDART_CHECKED_CSV/S12'
#     processed_dir='/Users/nunn/lunar_data/PDART_PROCESSED'
#     filelist='/Users/nunn/lunar_data/PDART_CHECKED_CSV/S12/filelist.pse.a12.1.29_05'
#     call_csv_import(checked_dir=checked_dir,processed_dir=processed_dir,filelist=filelist)
#
# def view_exampletest():
#     stream = read('/Users/nunn/lunar_data/PDART_PROCESSED/1970/XA/S12/pse.a12.1.99_11*.gz')
#
#     plot_from_stream(
#       stream=stream,
#       stations=['S12'],
#       channels=['AFR','ATT','MH1'],
#       merge_locations=False,
#       plot_type='normal',
#       time_interval=3,
#       outfile=None,
#       show=True)


if __name__ == "__main__":

    run_csv_check_work_tapes()

