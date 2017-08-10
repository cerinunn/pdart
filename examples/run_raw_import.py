#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to run the download.

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from obspy.core.utcdatetime import UTCDateTime
from pdart.raw_import import call_raw_import

def run_call_raw_import():

    start_time=UTCDateTime('1977-09-30T18:00:00.000000Z')
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z')

    csv_dir='/Users/nunn/lunar_data/JAXA_CSV2'
    raw_dir='/Users/nunn/lunar_data/JAXA_RAW2'

    call_raw_import(
        stations = ['S11','S12','S14','S15','S16'],
        csv_dir=csv_dir,
        raw_dir=raw_dir,
        start_time=start_time,
        end_time=end_time)

    # or run all dates and stations
    # call_raw_import(
    #     csv_dir=csv_dir,
    #     raw_dir=raw_dir)

if __name__ == "__main__":
    run_call_raw_import()
