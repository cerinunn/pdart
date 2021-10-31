#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to run initial checks on the csv files.

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
import pdart.config as config

import logging
# logging.handlers


import pandas as pd
import numpy as np

def run_csv_check_work_tapes():
    checked_dir='/Users/cnunn/lunar_data/PDART_CSV/S14'
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'


    # # these are all reset : 
    # filenames_resets = [
# pse.a12.1.35.csv.gz.log
# pse.a12.1.37.csv.gz.log
# pse.a12.1.39.csv.gz.log
# pse.a12.1.39.csv.gz.log
# pse.a12.1.46.csv.gz.log
# pse.a12.1.5.csv.gz.log
# pse.a12.1.5.csv.gz.log
# pse.a12.1.68.csv.gz.log
# pse.a12.1.68.csv.gz.log
# pse.a12.1.77.csv.gz.log
# pse.a12.1.83.csv.gz.log
# pse.a12.1.84.csv.gz.log
# pse.a12.1.9.csv.gz.log
# pse.a12.1.9.csv.gz.log
# pse.a12.2.126.csv.gz.log
# pse.a12.1.15.csv.gz.log
# pse.a12.1.21.csv.gz.log
    # ]

# still something wrong with these - after resetting
# failures:
# pse.a12.1.27.csv.gz.log
# pse.a12.1.41.csv.gz.log
# pse.a12.1.41.csv.gz.log

#     filenames_resets = [
# 
# 
#     'pse.a12.1.35.csv.gz',
#     'pse.a12.1.37.csv.gz',
#     'pse.a12.1.39.csv.gz',
#     'pse.a12.1.39.csv.gz',
#     'pse.a12.1.46.csv.gz',
#     'pse.a12.1.5.csv.gz',
#     'pse.a12.1.5.csv.gz',
#     'pse.a12.1.68.csv.gz',
#     'pse.a12.1.68.csv.gz',
#     'pse.a12.1.77.csv.gz',
#     'pse.a12.1.83.csv.gz',
#     'pse.a12.1.84.csv.gz',
#     'pse.a12.1.9.csv.gz',
#     'pse.a12.1.9.csv.gz',
#     'pse.a12.2.126.csv.gz',
#     'pse.a12.1.15.csv.gz',
#     'pse.a12.1.21.csv.gz',
# 
#     'pse.a12.1.27.csv.gz',
#     'pse.a12.1.41.csv.gz',
#     'pse.a12.1.41.csv.gz',
# 
# 
#     'pse.a12.1.27.csv.gz',
#     'pse.a12.1.41.csv.gz',
# 'pse.a12.1.103.csv.gz',
# 'pse.a12.1.14.csv.gz',
# 'pse.a12.1.146.csv.gz',
# 'pse.a12.1.15.csv.gz',
# 'pse.a12.1.21.csv.gz',
# 'pse.a12.1.27.csv.gz',
# 'pse.a12.1.31.csv.gz',
# 'pse.a12.1.33.csv.gz',
# 'pse.a12.1.35.csv.gz',
# 'pse.a12.1.37.csv.gz',
# 'pse.a12.1.39.csv.gz',
# 'pse.a12.1.41.csv.gz',
# 'pse.a12.1.43.csv.gz',
# 'pse.a12.1.43.csv.gz',
# 'pse.a12.1.44.csv.gz',
# 'pse.a12.1.45.csv.gz',
# 'pse.a12.1.46.csv.gz',
# 'pse.a12.1.5.csv.gz',
# 'pse.a12.1.68.csv.gz',
# 'pse.a12.1.77.csv.gz',
# 'pse.a12.1.8..csv.gz',
# 'pse.a12.1.83.csv.gz',
# 'pse.a12.1.84.csv.gz',
# 'pse.a12.1.9.csv.gz',
# 'pse.a12.2.120.csv.gz',
# 'pse.a12.2.126.csv.gz',
# 'pse.a12.2.132.csv.gz',
# 'pse.a12.2.47.csv.gz',
# 'pse.a12.2.60.csv.gz',
# 'pse.a12.2.91.csv.gz',
# 'pse.a12.2.97.csv.gz',
# 'pse.a12.3.125.csv.gz',
# 'pse.a12.3.133.csv.gz',
# 'pse.a12.3.45.csv.gz',
# 'pse.a12.3.58.csv.gz',
# 'pse.a12.5.119.csv.gz',
# 'pse.a12.8.43.csv.gz',
# ]

# exception_filenames = [
# 'pse.a12.1.162.csv.gz',
# 'pse.a12.1.7.csv.gz',
# 'pse.a12.2.171.csv.gz',
# 'pse.a12.2.49.csv.gz',
# 'pse.a12.3.139.csv.gz',
# 'pse.a12.3.139.csv.gz',
# 'pse.a12.4.18.csv.gz',
# 'pse.a12.4.67.csv.gz',
# 'pse.a12.5.142.csv.gz',
# 'pse.a12.6.222.csv.gz',
# 'pse.a12.6.255.csv.gz',
# 'pse.a12.7.59.csv.gz',
# 'pse.a12.7.76.csv.gz',
# 'pse.a12.8.188.csv.gz',
# 'pse.a12.9.256.csv.gz',
# 'pse.a12.9.256.csv.gz',
# 'pse.a12.9.76.csv.gz',
# ]

    # # example 1 with clock flag
    config.initial = False

    # from obspy.core.utcdatetime import UTCDateTime
    # overall_start_times = []
    # overall_start_times.append(UTCDateTime('1970-01-01'))
    # overall_start_times.append(UTCDateTime('1982-01-02'))
    # overall_start_times.append(UTCDateTime('1970-01-03'))
    # 
    # earlier = False
    # for i, st in enumerate(overall_start_times):
    #     if i > 0:
    #         print(st,overall_start_times[i-1])
    #         if st < overall_start_times[i-1]:
    #             earlier = True
    #         else:
    #             print('not')
    # if earlier == True:
    #     print('true')
    #     logging.info('SEVERE: File (overall) possibly contains records with the wrong timing')
    # 
    # exit()

    filenames=[


    # 'pse.a12.7.25.csv.gz',
    # 'pse.a12.7.26.csv.gz',

    # 'pse.a14.5.161.csv.gz',
    # 'pse.a14.5.162.csv.gz',
    # 
    # 'pse.a15.4.168.csv.gz',
    # 'pse.a15.4.169.csv.gz',
    # 
    # 'pse.a16.3.87.csv.gz',
    # 'pse.a16.3.88.csv.gz',

    # 'pse.a14.1.2.csv.gz',
    # 'pse.a14.1.3.csv.gz'
    'pse.a14.5.161.csv.gz'
    # 'pse.a14.5.162.csv.gz'

    # 'pse.a12.1.12.csv.gz',
    # 'pse.a12.1.13.csv.gz',
    # 'pse.a12.1.60.csv.gz',
    # 'pse.a12.1.95.csv.gz'
    # 'pse.a12.2.84.csv.gz'

# pse.a12.1.60.csv.gz.log:SEVERE: Ground station/Station - 10 S12 Dropping the section (45699 records)

# pse.a12.1.95.csv.gz.log:SEVERE: Ground station/Station - 2 S12 Dropping the section (34549 records)

# pse.a12.9.178.S12.13.1975.179.00_00_00_239821.csv.gz
    # 'pse.a16.1.83XXX.csv.gz'
    ]
    call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,log_dir=processed_dir,
      filenames=filenames,logging_level=logging.DEBUG
        # ,single_station='S12'
        # ,single_ground_station=10
        )

    # remember to rerun this at the end 

# 10 12



if __name__ == "__main__":

    run_csv_check_work_tapes()

