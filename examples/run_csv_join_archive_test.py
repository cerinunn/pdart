#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to create PDART_ARCHIVE_TEST - used to 
create a subset of the MSEED files with the same 
structure as the whole archive. 

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
from obspy.core.stream import read
import pdart.config as config

import logging
# logging.handlers

import pandas as pd
import numpy as np

def run_csv_join_work_tapes():
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'
    join_dir='/Users/cnunn/lunar_data/PDART_ARCHIVE_TEST'
    config.combine_ground_stations=True
    config.clean_spikes=True

    config.view_corrected_traces = False
    call_csv_join_work_tapes(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    year_start=1976,
    year_end=1976,
    day_start=60,
    day_end=62,
    stations=['S12','S14','S15','S16'],
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
    day_start=202,
    day_end=202,
    stations=['S11'],
    manual_clock_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_clock_fix.csv',
    manual_jump_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_jump_fix.csv',
    manual_exclude='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_exclude.csv',
    manual_grab_before='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_before.csv',
    manual_grab_after='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_after.csv',
    logging_level=logging.INFO)

if __name__ == "__main__":

    run_csv_join_work_tapes()
