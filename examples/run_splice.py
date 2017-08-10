#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to run the splicing.

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
from pdart_apollo.splice import call_splice_chains

def run_call_build_chains():

    start_time=UTCDateTime('1977-09-30T18:00:00.000000Z')
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z')

    chain_dir='/Users/nunn/lunar_data/JAXA_CHAINS2'
    splice_dir='/Users/nunn/lunar_data/JAXA_SPLICE2'

    # It is very important to run this in date order - the timing
    # of each trace depends on the previous one.

    call_splice_chains(
      stations=['S12'],
      chain_dir=chain_dir,
      splice_dir=splice_dir,
      start_time=start_time,
      end_time=end_time)

    # or run all dates and stations
    # call_splice_chains(
    #   chain_dir=chain_dir,
    #   splice_dir=splice_dir)

if __name__ == "__main__":
    run_call_build_chains()
