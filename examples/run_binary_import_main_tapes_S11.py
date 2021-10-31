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
# from pdart.binary_import_work_tapes import import_csv_work_tapes
from pdart.binary_import_main_tapes import import_csv_main_tapes

if __name__ == "__main__":

    base_dir='/Users/cnunn/lunar_data/PDART_TAPES/S11'
    out_base_dir='/Users/cnunn/lunar_data/PDART_CSV/S11'

    filenames_all = [
'pse.a11.1.1.gz',
'pse.a11.1.10.gz',
'pse.a11.1.11.gz',
'pse.a11.1.12.gz',
'pse.a11.1.13.gz',
'pse.a11.1.14.gz',
'pse.a11.1.15.gz',
'pse.a11.1.16.gz',
'pse.a11.1.17.gz',
'pse.a11.1.18.gz',
'pse.a11.1.19.gz',
'pse.a11.1.2.gz',
'pse.a11.1.20.gz',
'pse.a11.1.21.gz',
'pse.a11.1.22.gz',
'pse.a11.1.3.gz',
'pse.a11.1.4.gz',
'pse.a11.1.5.gz',
'pse.a11.1.6.gz',
'pse.a11.1.7.gz',
'pse.a11.1.8.gz',
'pse.a11.1.9.gz',

]

    import_csv_main_tapes(base_dir=base_dir, out_base_dir=out_base_dir,filenames=filenames_all)