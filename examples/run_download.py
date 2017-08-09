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
from pdart.download import download_from_jaxa

def run_download_from_jaxa():
    '''
    Example method to run the download.
    '''
    base_dir='/Users/nunn/lunar_data/JAXA_CSV'

    start_time=UTCDateTime('1977-09-30T18:00:00.000000Z')
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z')
    download_from_jaxa(base_dir=base_dir, start_time=start_time,
      end_time=end_time)

    # or run all dates
    # download_from_jaxa(base_dir=base_dir)

if __name__ == "__main__":
    run_download_from_jaxa()
