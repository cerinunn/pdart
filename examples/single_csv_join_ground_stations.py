#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to run initial checks on the csv files and extract 
to separate ground stations (normally they are combined)

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
from pdart.extra_plots.plot_timing_divergence import plot_timing_from_dir
from obspy.core.utcdatetime import UTCDateTime 


import logging
# logging.handlers

import pandas as pd
import numpy as np


def run_csv_join_work_tapes():
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'
    join_dir='/Users/cnunn/lunar_data/PDART_SEPARATE_GROUND_STATIONS'
    log_dir='/Users/cnunn/lunar_data/tmp_PDART_LOG'
    config.combine_ground_stations=False
    config.clean_spikes=True
    print('MAKE SURE THAT THESE ARE RERUN PROPERLY ')

    stations=['S12']
    year=1973
    # day=1

    print('log dir: ', log_dir)

    run_single=True
    plot_timing=True

    if run_single:
        config.view_corrected_traces = False
        config.fix_clock_error = True
        config.fix_jump_error = True
        call_csv_join_work_tapes(
        processed_dir=processed_dir,
        join_dir=join_dir,
        log_dir=log_dir,
        year_start=year,
        year_end=year,
        day_start=1,
        day_end=8,
        stations=stations,
        manual_clock_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_clock_fix.csv',
        manual_jump_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_jump_fix.csv',
        manual_exclude='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_exclude.csv',
        manual_grab_before='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_before.csv',
        manual_grab_after='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_after.csv',
        logging_level=logging.DEBUG)

    if plot_timing:
        plot_timing_from_dir(top_level_dir=join_dir, start_time=UTCDateTime(year=year,julday=day), stations=stations, include_line=True, out_dir='../extra_plots_output', save_fig=False, plot_fig=True)






if __name__ == "__main__":
    run_csv_join_work_tapes()
