#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example file to import from the binary files. 

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from pdart.binary_import_work_tapes import import_csv_work_tapes

if __name__ == "__main__":

    base_dir='/Users/cnunn/lunar_data/PDART_TAPES/work_tapes/'
    out_base_dir='/Users/cnunn/lunar_data/PDART_CSV_WORK_TAPES'
    out_base_dir_full='/Users/cnunn/lunar_data/PDART_CSV_WORK_TAPES_FULL'

    import_csv_work_tapes(base_dir=base_dir, out_base_dir=out_base_dir_full,filenames=['wtn.19.45.gz'])


