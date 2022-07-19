#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run csv_join_work_tapes for 1973. 
Includes some days which need rerunning because the automatic 
code wasn't  working for them.

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from pdart.csv_join_work_tapes import call_csv_join_work_tapes
# from pdart.view import plot_from_stream
from obspy.core.stream import read
import pdart.config as config

import logging
# logging.handlers

import pandas as pd
import numpy as np


def run_csv_join_work_tapes():
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'
    join_dir='/Users/cnunn/lunar_data/PDART_V2'
    config.combine_ground_stations=True
    config.clean_spikes=True
    print('MAKE SURE THAT THESE ARE RERUN PROPERLY ')

    config.view_corrected_traces = False
    config.fix_clock_error = True
    config.fix_jump_error = True
    call_csv_join_work_tapes(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    year_start=1973,
    year_end=1973,
    day_start=1,
    day_end=366,
    stations=['S12','S14','S15','S16'],
    manual_clock_correction='../manual_fix_files/manual_clock_fix.csv',
    manual_jump_correction='../manual_fix_files/manual_jump_fix.csv',
    manual_exclude='../manual_fix_files/manual_exclude.csv',
    manual_grab_before='../manual_fix_files/manual_grab_before.csv',
    manual_grab_after='../manual_fix_files/manual_grab_after.csv',
    logging_level=logging.INFO)

    # auto correction not working for this day - note where log file is 
    config.fix_jump_error = False
    call_csv_join_work_tapes(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir='/Users/cnunn/lunar_data/tmp_PDART_LOG',
    year_start=1973,
    year_end=1973,
    day_start=101,
    day_end=101,
    stations=['S16'],
    manual_clock_correction='../manual_fix_files/manual_clock_fix.csv',
    manual_jump_correction='../manual_fix_files/manual_jump_fix.csv',
    manual_exclude='../manual_fix_files/manual_exclude.csv',
    manual_grab_before='../manual_fix_files/manual_grab_before.csv',
    manual_grab_after='../manual_fix_files/manual_grab_after.csv',
    logging_level=logging.INFO)

    # auto correction not working for this day - note where log file is 
    config.fix_jump_error = False
    call_csv_join_work_tapes(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir='/Users/cnunn/lunar_data/tmp_PDART_LOG',
    year_start=1973,
    year_end=1973,
    day_start=87,
    day_end=87,
    stations=['S16'],
    manual_clock_correction='../manual_fix_files/manual_clock_fix.csv',
    manual_jump_correction='../manual_fix_files/manual_jump_fix.csv',
    manual_exclude='../manual_fix_files/manual_exclude.csv',
    manual_grab_before='../manual_fix_files/manual_grab_before.csv',
    manual_grab_after='../manual_fix_files/manual_grab_after.csv',
    logging_level=logging.INFO)

    # auto correction not working for this day - note where log file is 
    config.fix_jump_error = False
    call_csv_join_work_tapes(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir='/Users/cnunn/lunar_data/tmp_PDART_LOG',
    year_start=1973,
    year_end=1973,
    day_start=243,
    day_end=243,
    stations=['S16'],
    manual_clock_correction='../manual_fix_files/manual_clock_fix.csv',
    manual_jump_correction='../manual_fix_files/manual_jump_fix.csv',
    manual_exclude='../manual_fix_files/manual_exclude.csv',
    manual_grab_before='../manual_fix_files/manual_grab_before.csv',
    manual_grab_after='../manual_fix_files/manual_grab_after.csv',
    logging_level=logging.INFO)



if __name__ == "__main__":

    run_csv_join_work_tapes()
