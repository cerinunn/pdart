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
# from pdart.csv_check_work_tapes_jul2021 import call_csv_check_work_tapes
from pdart.csv_check_work_tapes import call_csv_check_work_tapes
# from pdart.view import plot_from_stream
from obspy.core.stream import read
import pdart.config as config



import logging
# logging.handlers


import pandas as pd
import numpy as np

def run_csv_check_work_tapes():
    checked_dir='/Users/cnunn/lunar_data/PDART_CSV/work_tapes'
    checked_dir_full='/Users/cnunn/lunar_data/PDART_CSV/work_tapes'
    processed_dir='/Users/cnunn/lunar_data/tmp_PDART_PROCESSED'

    print('using the new one')

    # fix the negative time gaps 
    config.initial = False



    filenames=[
    # 'wtn.10.10.csv.gz',
    # 'wtn.10.43.csv.gz',
    # 'wtn.10.31.csv.gz',
    # 'wtn.20.12.csv.gz',
    'wtn.20.14.csv.gz',
# wtn.10.10.csv.gz.log:WARNING: Dropping 508 repeated lines
# wtn.10.43.csv.gz.log:WARNING: Dropping 684 repeated lines
# wtn.10.31.csv.gz.log:WARNING: Dropping 2 repeated lines
# wtn.20.12.csv.gz.log:WARNING: Dropping 3791 repeated lines
# wtn.20.14.csv.gz.log:WARNING: Dropping 3707 repeated lines
    ]


    logging_level = logging.DEBUG
    # logging_level = logging.INFO
# WARNING: Ground station/Station - 7 S16 Unable to fix the timing
    call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=processed_dir,
      filenames=filenames,logging_level=logging_level
        
        # ,single_station='S15'
        # ,single_ground_station=4

# WARNING: Ground station/Station - 10 S15

        )



if __name__ == "__main__":

    run_csv_check_work_tapes()

