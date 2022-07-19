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
from pdart.csv_check_work_tapes import call_csv_check_work_tapes
# from pdart.view import plot_from_stream
from obspy.core.stream import read
import pdart.config as config

import logging
# logging.handlers


import pandas as pd
import numpy as np

def run_csv_check_work_tapes():
    checked_dir='/Users/cnunn/lunar_data/tmp_PDART_CSV/S12'
    processed_dir='/Users/cnunn/lunar_data/tmp_PDART_PROCESSED'

    config.initial = False

    # can change the parameter in the config file
    config.cumsum_final_test = 180

    # In this case, just check a single station and a single ground station
    # (useful because the file are very large)
    filenames=[
        'pse.a11.1.1.csv.gz',
    ]
    call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,log_dir=processed_dir,
      filenames=filenames,logging_level=logging.DEBUG
        ,single_station='S12'
        ,single_ground_station=11
        )

if __name__ == "__main__":

    run_csv_check_work_tapes()

