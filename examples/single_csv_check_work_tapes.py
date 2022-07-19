#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to check the work tapes (also works for main tapes)

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
    'wtn.20.14.csv.gz',
    ]

    # Only check data for S15 and ground station 4 - useful because the file is 
    # big - use this to check everything is working ok. 
    logging_level = logging.DEBUG
    call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=processed_dir,
      filenames=filenames,logging_level=logging_level
        ,single_station='S15'
        ,single_ground_station=4
        )

if __name__ == "__main__":

    run_csv_check_work_tapes()

