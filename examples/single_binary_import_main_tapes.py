#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to import the csv files.

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from pdart.binary_import_main_tapes import import_csv_main_tapes
from pdart.binary_import_main_tapes_optim import import_csv_main_tapes as import_csv_main_tapes_optim
from datetime import datetime, timedelta

if __name__ == "__main__":




    t1 = datetime.now()

    # base_dir='/Users/cnunn/lunar_data/PDART_TAPES/S14'
    # out_base_dir_test='/Users/cnunn/lunar_data/del_PDART_test'
    # import_csv_main_tapes(base_dir=base_dir, out_base_dir=out_base_dir_test,filenames=['pse.a14.9.99.gz'])

    base_dir='/Users/cnunn/lunar_data/PDART_TAPES/S12'
    out_base_dir_test='/Users/cnunn/lunar_data/del_PDART_test'
    out_base_dir='/Users/cnunn/lunar_data/PDART_CSV/S12'

    import_csv_main_tapes(base_dir=base_dir, out_base_dir=out_base_dir_test,filenames=['pse.a12.1.43.gz'])
    # import_csv_main_tapes(base_dir=base_dir, out_base_dir=out_base_dir,filenames=['pse.a12.1.43.gz'])
    # 
    # 
    # t2 = datetime.now() 
    # print('t2 test : {}'.format((t2-t1).total_seconds())) 
    # # t2 test : 3.15166
    # 638K

    # import_csv_main_tapes(base_dir=base_dir, out_base_dir=out_base_dir,filenames=filenames_all)
    # new_format = 0 example

    # # test the optimized one 
    # t1 = datetime.now()
    # out_base_dir_optim='/Users/cnunn/lunar_data/del_PDART_test_optim'
    # import_csv_main_tapes_optim(base_dir=base_dir, out_base_dir=out_base_dir_optim,filenames=['pse.a12.1.140.gz'])
    # t2 = datetime.now() 
    # print('t2 optim : {}'.format((t2-t1).total_seconds())) 
    # t2 optim : 0.996489
    # 412K
    # doesn't contain shz

    # # test the optimized one 
    # t1 = datetime.now()
    # out_base_dir_optim='/Users/cnunn/lunar_data/del_PDART_test_optim2'
    # import_csv_main_tapes_optim(base_dir=base_dir, out_base_dir=out_base_dir_optim,filenames=['pse.a12.1.140.gz'])
    # t2 = datetime.now() 
    # print('t2 optim : {}'.format((t2-t1).total_seconds())) 


    # new_format = 1 example
    # import_csv_main_tapes(base_dir=base_dir, out_base_dir=out_base_dir,filenames=['pse.a12.1.144.gz'])


