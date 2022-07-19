#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Final file to run initial checks on csv files

Make checked csv files (will drop some damaged data)
out of the raw csv files extracted from binary. 

This file will run with config.cumsum_final_test = 10.
This will recover more data, but in smaller 
file chucks - and they will need checking manaully.
Used for important events where the timing was 
sporadic.

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

    # processed_dir='/Users/cnunn/lunar_data/tmp_PDART_PROCESSED'
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'



# 47 1972 213 1972-07-31T18:08:15.720000Z S15 - pse.a15.3.8.S15 - fixed 
# 48 1972 213 1972-07-31T18:08:15.720000Z S16 - pse.a16.1.101.S16   - not fixed, investigate 

# 51 1972 242 1972-08-29T22:58:33.570000Z S15 - pse.a15.3.38.S15 - not fixed, ignore
# 52 1972 242 1972-08-29T22:58:33.570000Z S16 - pse.a16.1.131.S16, not fixed, ignore 
# 88 1975 124 1975-05-04T09:59:28.990000Z S16 - pse.a16.7.49.S16, not fixed, investigate briefly
# 99 1976 319 1976-11-14T23:13:06.670000Z S15 - wtn.9.47.S15.9, not an error 
# 103 1977 107 1977-04-17T23:32:06.030000Z S15, not fixed, ignore 
# 111 1972 261 1972-09-17T14:35:03.000000Z S15 - pse.a15.3.55 , not an event

    config.cumsum_final_test = 10

    checked_dir='/Users/cnunn/lunar_data/PDART_CSV/S12'
    filenames=[
        # large impact on this day 
        # 'pse.a12.3.96.csv.gz'
    ]
    call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,log_dir=processed_dir,
      filenames=filenames,logging_level=logging.INFO
        # ,single_station='S12'
        # ,single_ground_station=11
        )

    checked_dir='/Users/cnunn/lunar_data/PDART_CSV/S15'
    filenames=[
'pse.a15.3.8.csv.gz',
'pse.a15.3.38.csv.gz',
'pse.a15.3.55.csv.gz'
    ]
    call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,log_dir=processed_dir,
      filenames=filenames,logging_level=logging.INFO
        # ,single_station='S12'
        # ,single_ground_station=11
        )

    checked_dir='/Users/cnunn/lunar_data/PDART_CSV/S16'
    filenames=[
'pse.a16.1.101.csv.gz',
'pse.a16.1.131.csv.gz',
'pse.a16.7.49.csv.gz',
    ]
    call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,log_dir=processed_dir,
      filenames=filenames,logging_level=logging.INFO
        # ,single_station='S12'
        # ,single_ground_station=11
        )








if __name__ == "__main__":

    run_csv_check_work_tapes()

