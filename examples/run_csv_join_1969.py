#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run csv_join_work_tapes for 1969.

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
    wildcard_style='pse',
    year_start=1969,
    year_end=1969,
    day_start=1,
    day_end=322,
    stations=['S11'],
    manual_clock_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_clock_fix.csv',
    manual_jump_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_jump_fix.csv',
    manual_exclude='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_exclude.csv',
    manual_grab_before='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_before.csv',
    manual_grab_after='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_after.csv',
    logging_level=logging.INFO)
    
    config.view_corrected_traces = False
    call_csv_join_work_tapes(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    wildcard_style='pse',
    year_start=1969,
    year_end=1969,
    day_start=323,
    day_end=366,
    stations=['S12'],
    manual_clock_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_clock_fix.csv',
    manual_jump_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_jump_fix.csv',
    manual_exclude='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_exclude.csv',
    manual_grab_before='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_before.csv',
    manual_grab_after='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_after.csv',
    logging_level=logging.INFO)

if __name__ == "__main__":

    run_csv_join_work_tapes()
