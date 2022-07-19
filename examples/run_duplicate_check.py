#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Call check_duplicates -  check there are no repeated traces.

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from pdart.check_duplicates import call_check_duplicates
# from pdart.view import plot_from_stream
from obspy.core.stream import read
import pdart.config as config

import logging
# logging.handlers

import pandas as pd
import numpy as np


def run_call_check_duplicates():

    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'
    join_dir='/Users/cnunn/lunar_data/PDART_CORRECTED'

    call_check_duplicates(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    year_start=1971,
    year_end=1971,
    day_start=1,
    day_end=366,
    stations=['S12','S14','S15','S16'],
    logging_level=logging.INFO)

    call_check_duplicates(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    year_start=1972,
    year_end=1972,
    day_start=1,
    day_end=366,
    stations=['S12','S14','S15','S16'],
    logging_level=logging.INFO)

    call_check_duplicates(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    year_start=1973,
    year_end=1973,
    day_start=1,
    day_end=366,
    stations=['S12','S14','S15','S16'],
    logging_level=logging.INFO)

    call_check_duplicates(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    year_start=1974,
    year_end=1974,
    day_start=1,
    day_end=365,
    stations=['S12','S14','S15','S16'],
    logging_level=logging.INFO)

    call_check_duplicates(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    year_start=1975,
    year_end=1975,
    day_start=1,
    day_end=366,
    stations=['S12','S14','S15','S16'],
    logging_level=logging.INFO)

    call_check_duplicates(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    year_start=1976,
    year_end=1976,
    day_start=1,
    day_end=366,
    stations=['S12','S14','S15','S16'],
    logging_level=logging.INFO)

    call_check_duplicates(
    processed_dir=processed_dir,
    join_dir=join_dir,
    log_dir=join_dir,
    year_start=1977,
    year_end=1977,
    day_start=1,
    day_end=273,
    stations=['S12','S14','S15','S16'],
    logging_level=logging.INFO)















if __name__ == "__main__":

    run_call_check_duplicates()
