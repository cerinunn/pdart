#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to join the work tapes (works on work tapes and main tapes).


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
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES'
    join_dir='/Users/cnunn/lunar_data/PDART_SEPARATE_GROUND_STATIONS'
    config.combine_ground_stations=False
    config.clean_spikes=False


    config.view_corrected_traces = True
    call_csv_join_work_tapes(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=processed_dir,
    year_start=1976,
    year_end=1976,
    day_start=64,
    day_end=64,
    stations=['S12'],
    logging_level=logging.DEBUG)





if __name__ == "__main__":

    run_csv_join_work_tapes()
