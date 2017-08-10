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
from pdart_apollo.view import plot_from_file, plot_from_stream


def run_plot_from_file():

    start_time=UTCDateTime('1977-09-30T18:00:00.000000Z')
    end_time=UTCDateTime('1977-09-30T21:00:00.000000Z')

    plot_from_file(
      stations=['S12'],
      channels=['AFR'],
      file_dir='/Users/nunn/lunar_data/JAXA_SPLICE2',
      start_time=start_time,
      end_time=end_time,
      merge_locations=True,
      plot_type='normal',
      time_interval=3,
      read_gzip=True)


if __name__ == "__main__":
    run_plot_from_file()
