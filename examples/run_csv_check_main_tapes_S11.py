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
from obspy.core.stream import read
import pdart.config as config

import logging
# logging.handlers


import pandas as pd
import numpy as np

def run_csv_check_work_tapes():
    checked_dir='/Users/cnunn/lunar_data/PDART_CSV/S11'
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'

    # example 1 with clock flag
    config.initial = False
    filenames=[
'pse.a11.1.1.csv.gz',
'pse.a11.1.10.csv.gz',
'pse.a11.1.11.csv.gz',
'pse.a11.1.12.csv.gz',
'pse.a11.1.13.csv.gz',
'pse.a11.1.14.csv.gz',
'pse.a11.1.15.csv.gz',
'pse.a11.1.16.csv.gz',
'pse.a11.1.17.csv.gz',
'pse.a11.1.18.csv.gz',
'pse.a11.1.19.csv.gz',
'pse.a11.1.2.csv.gz',
'pse.a11.1.20.csv.gz',
'pse.a11.1.21.csv.gz',
'pse.a11.1.22.csv.gz',
'pse.a11.1.3.csv.gz',
'pse.a11.1.4.csv.gz',
'pse.a11.1.5.csv.gz',
'pse.a11.1.6.csv.gz',
'pse.a11.1.7.csv.gz',
'pse.a11.1.8.csv.gz',
'pse.a11.1.9.csv.gz',

    ]
    call_csv_check_work_tapes(checked_dir=checked_dir,processed_dir=processed_dir,log_dir=checked_dir,
      filenames=filenames,logging_level=logging.INFO
        # ,single_station='12'
        # ,single_ground_station=10
        )

    print('Also remember to run final_csv_check_main_tapes_all.py!!!!!')
    # final_csv_check_main_tapes_all.py contains a few files with specific settings 
    # when the timing wasn't that great





if __name__ == "__main__":

    run_csv_check_work_tapes()

