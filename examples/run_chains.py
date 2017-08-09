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
from pdart.chains import call_build_chains

def run_call_build_chains():

    start_time=UTCDateTime('1977-09-30T18:00:00.000000Z')
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z')

    raw_dir='/Users/nunn/lunar_data/JAXA_RAW2'
    chain_dir='/Users/nunn/lunar_data/JAXA_CHAINS2'

    call_build_chains(
        stations=['S12'],
        raw_dir=raw_dir,
        chain_dir=chain_dir,
        start_time=start_time,
        end_time=end_time)

    # or run all dates and stations
    # call_build_chains(
    #     raw_dir=raw_dir,
    #     chain_dir=chain_dir)

if __name__ == "__main__":
    run_call_build_chains()
