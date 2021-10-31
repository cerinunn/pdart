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
    join_dir='/Users/cnunn/lunar_data/PDART'
    config.combine_ground_stations=True
    config.clean_spikes=True

    # processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED_MAIN_TAPES'
    # join_dir='/Users/cnunn/lunar_data/PDART_CONTINUOUS_MAIN_TAPES_ALL_LOCATIONS'
    # config.combine_ground_stations=False
    # config.clean_spikes=False

    config.view_corrected_traces = False
    call_csv_join_work_tapes(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=processed_dir,
    wildcard_style='pse',
    year_start=1970,
    year_end=1970,
    day_start=1,
    day_end=366,
# 238
    # stations=['S12','S14','S15','S16'],
    stations=['S12'],
    logging_level=logging.INFO)



if __name__ == "__main__":

    run_csv_join_work_tapes()
