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
from pdart.view import plot_from_stream
from obspy.core.stream import read
import pdart.config as config

import logging
# logging.handlers


import pandas as pd
import numpy as np

def run_csv_check_work_tapes():
    checked_dir='/Users/cnunn/lunar_data/PDART_CSV_WORK_TAPES'
    checked_dir_full='/Users/cnunn/lunar_data/PDART_CSV_WORK_TAPES_FULL'
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED_WORK_TAPES'

    config.initial = False

    # put back if this fails: XFAIL
# wtn.9.50.csv.gz
# wtn.1.10.csv.gz
# wtn.20.39.csv.gz
# wtn.19.55.csv.gz
# wtn.19.6.csv.gz

############ Ground station: 6 Station: S14 1976-03-05 05:50:51.276000+00:00
    call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=checked_dir_full,
      filenames=['wtn.1.10.csv.gz'],logging_level=logging.DEBUG
        # ,single_station='S14'
        # ,single_ground_station=6
    )

    exit()


    call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=checked_dir_full,
      filenames=['wtn.1.1.csv.gz'],logging_level=logging.DEBUG
        ,single_station='S16'
        ,single_ground_station=1
    )

    exit()

    # TEST 1 - works 
    # 
    # # 1) delete a crap section followed by a bit of a repeat
    call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=checked_dir_full,
      filenames=['wtn.9.50.csv.gz'],logging_level=logging.DEBUG
        ,single_station='S12'
        ,single_ground_station=7
    )
    # # df, rec_drop_damaged1 = all_drop_idx(df, gzip_filename, problem_gzip_filename='wtn.9.50.csv.gz', 
    # #           orig_idx_start=66, orig_idx_end=179)
    # 
    # exit()

    # this works
    # reset_all_ground_stations
    call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=checked_dir_full,
      filenames=['wtn.1.10.csv.gz'],logging_level=logging.DEBUG
        ,single_station='S12'
        ,single_ground_station=6
    )
    # # df = reset_all_ground_stations(df, gzip_filename, problem_gzip_filename='wtn.1.10.csv.gz',start_timestamp ='1976-03-05 01:09:47.019000+00:00',end_timestamp='1976-03-05 03:39:54.226000+00:00',start_station='S12',end_station='S17',start_orig_no=58358,end_orig_no=3970,orig_ground_station1=6)
    # # 58358,"1976-03-05T01:09:47.019000Z",41,"S12",0,6,525,483,4
    # 
    # exit()

    # 2) some -7777 errors
    call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=checked_dir_full,
      filenames=['wtn.20.39.csv.gz'],logging_level=logging.DEBUG
        ,single_station='S12'
        ,single_ground_station=5
    )

    # exit()  

    # 3) a real jump forward to a later time 
# 606,"1977-08-25T21:29:56.945000Z",60,"S17",0,2,0,0,0,514,508,0,0,0,0,0,0,0,"00011",59,"1110001001000011101101"
# 607,"1977-08-26T02:37:26.628000Z",1,"S12",0,2,514,513,513,514,513,513,514,513,513,514,513,513,"00011",20,"1110001001000011101101"
    # check it's not being manually removed 
    call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=checked_dir_full,
      filenames=['wtn.19.55.csv.gz'],logging_level=logging.DEBUG
        ,single_station='S12'
        ,single_ground_station=2
    )

    # there's severval clock jumps during a time period when clock error is set - could potentially be fixed
    # wtn.19.6.csv.gz.log:WARNING: wtn.19.6.csv.gz Removing 11189 known damaged record(s)

    # call_csv_check_work_tapes(checked_dir=checked_dir_full,processed_dir=processed_dir,log_dir=checked_dir_full,
    #   filenames=['wtn.19.6.csv.gz'],logging_level=logging.DEBUG
    #     ,single_station='S15'
    #     ,single_ground_station=5
    # )

    # do something with the repeating sections currently set up as new ground stations




if __name__ == "__main__":

    run_csv_check_work_tapes()

