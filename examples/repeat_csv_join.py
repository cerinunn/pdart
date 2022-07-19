#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Rerun csv_join_work_tapes

I reran csv_join_work_tapes for some files that I wasn't sure 
were working correctly. No need to run this file again.

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
    # processed_dir='/Users/cnunn/lunar_data/tmp_PDART_PROCESSED'
    processed_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'
    # join_dir='/Users/cnunn/lunar_data/tmp_PDART'
    join_dir='/Users/cnunn/lunar_data/PDART_CORRECTED'
    log_dir='/Users/cnunn/lunar_data/tmp_PDART_LOG'
    # log_dir='/Users/cnunn/lunar_data/PDART_CORRECTED'
    config.combine_ground_stations=True
    config.clean_spikes=True
    print('MAKE SURE THAT THESE ARE RERUN PROPERLY ')

    print('join dir: ', join_dir)
    print('log dir: ', log_dir)

    rerun_list = [
    ('S14',1971,51),
    ('S12',1971,71),
    ('S12',1971,136),
    ('S12',1971,312),

    ('S14',1971,223),

    ('S12',1973,55),
    ('S16',1973,62),
    ('S14',1973,69),
    ('S14',1973,80),

    ('S15',1973,109),



    ('S12',1973,116),


    ('S15',1973,120),

    ('S12',1973,143),
    ('S12',1973,159),

    ('S16',1973,208),

    ('S15',1973,221),

    ('S16',1973,224),

    ('S12',1973,274),

    ('S15',1973,285),

    ('S14',1973,294),
    ('S16',1973,294),
    ('S12',1973,294),

    ('S14',1973,331),

    ('S12',1973,338),

    ('S14',1974,46),

    ('S12',1974,129),

    ('S14',1974,317),

    ('S12',1975,34),

    ('S15',1976,37),
    ('S12',1976,92),
    ('S15',1976,92),

    ('S12',1976,93),
    ('S15',1976,93),

]

    run_single=True
    plot_timing=False

    for rerun in rerun_list:
        stations=[rerun[0]]
        year = rerun[1]
        day = rerun[2]

        print(stations,year, day)

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
            day_start=day,
            day_end=day,
            stations=stations,
            manual_clock_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_clock_fix.csv',
            manual_jump_correction='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_jump_fix.csv',
            manual_exclude='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_exclude.csv',
            manual_grab_before='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_before.csv',
            manual_grab_after='/Users/cnunn/lunar_data/PDART_MANUAL_FIX/manual_grab_after.csv',
            logging_level=logging.INFO)

        # if plot_timing:
        #     plot_timing_from_dir(top_level_dir=join_dir, start_time=UTCDateTime(year=year,julday=day), stations=stations, include_line=True, out_dir='../extra_plots_output', save_fig=False, plot_fig=True)






if __name__ == "__main__":
    # stream = read('/Users/cnunn/lunar_data/PDART_PDS_OUTPUT/data/xa/continuous_waveform/s12/1976/061/xa.s12..att.1976.061.0.mseed')
    # print(stream)
    # 
    # stream = read('/Users/cnunn/lunar_data/PDART_PDS_OUTPUT/data/xa/continuous_waveform/s12/1976/061/xa.s12.01.mhz.1976.061.0.mseed')
    # print(stream)
    # 
    # stream = read('/Users/cnunn/lunar_data/PDART/S12/1976/061/XA.S12..ATT.1976.061.0.MSEED')
    # print(stream)
    # 
    # stream = read('/Users/cnunn/lunar_data/PDART/S12/1976/061/XA.S12.01.MHZ.1976.061.0.MSEED')
    # print(stream)

    # print(UTCDateTime('1977-07-22').julday)
    # exit()

    
    # stream = read('/Users/cnunn/lunar_data/tmp_PDART/s16/1973/184/XA.S16..ATT.1973.184.0.MSEED')
    # print(stream)
    # stream.merge(method=1)
    # print(stream)    
    # import numpy.ma as ma
    # 
    # timestamp0 = 113528005.647
    # timestamps = []
    # for i in range(0,10):
    #     timestamps.append(timestamp0 + i*0.6038)
    # 
    # timestamps = np.array(timestamps)
    # mask_true = np.array([True, True, True, True, True, 
    #     True, True, True, True, True ])
    # mask_begin = np.array([True, True, True, True, True, 
    #     False, False, True, True, True])
    # print(mask_begin)
    # 
    # # 
    # # mask_timestamp = np.ma.masked_where(tr_CLK.data==1, tr_CLK.data)  
    # mask_timestamp = ma.masked_array(timestamps, mask=mask_true)
    # 
    # if mask_timestamp[-1] is ma.masked:
    #     C_false = np.where( mask_begin == False )
    #     print(C_false)
    #     if len(C_false[0]) > 0: 
    #         last_good = C_false[0][-1]
    #         logging.info(last_good)
    #         print(mask_begin)
    #         mask_begin[last_good:] = False
    #         print('reprint')
    #         print(mask_begin)
    # 
    #             # XXXX
    #         if mask_CLK[-1] is ma.masked:
    #             logging.info('true')
    #             logging.info('SEVERE: Last record contains clock flag - {}.{}'.format(date1.year,date1.julday)) 
    #             # retrieve the valid values 
    # 
    #             C_false = np.where( mask_CLK == False )
    #             logging.info(C_false)
    #             if len(C_false[0]) > 0: 
    #                 last_good = C_false[-1]
    #                 logging.info(last_good)
    #                 mask_CLK[last_good:] = False 
    # 
    #                 logging.info(mask_CLK.data[-10:])   

    # print(mask_timestamp)

    run_csv_join_work_tapes()

