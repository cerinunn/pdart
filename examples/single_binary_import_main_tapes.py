#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to import the csv files from the main tapes - run one or more files. 

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
# from pdart.binary_import_main_tapes_optim import import_csv_main_tapes as import_csv_main_tapes_optim
from datetime import datetime, timedelta

if __name__ == "__main__":

    t1 = datetime.now()

    base_dir='/Users/cnunn/lunar_data/PDART_TAPES/S12'
    out_base_dir_test='/Users/cnunn/lunar_data/del_PDART_test'
    out_base_dir='/Users/cnunn/lunar_data/PDART_CSV/S12'

    out_base_dir=out_base_dir_test

    filelist = [
        'pse.a12.10.91.gz',]

    import_csv_main_tapes(base_dir=base_dir, out_base_dir=out_base_dir,filenames=filelist)


