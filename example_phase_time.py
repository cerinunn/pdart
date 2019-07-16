#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Helper file to update the XML manually 

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)



"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from pdart.diffusion.view_catalog_with_envelopes import (
  view_catalog_with_envelopes, view_section_with_envelopes)
from obspy.core.utcdatetime import UTCDateTime

# smi:jpl/filter/Apollo/standard_0.25_0.75
# smi:jpl/filter/Apollo/standard_0.75_1.25
# smi:jpl/filter/Apollo/standard_1.25_1.75

# phase_hint = t_onset
# phase_hint = t_max

  # <pick>
  #   <time>
  #     <value>1971-02-04T07:41:13.500000Z</value>
  #   </time>
  #   <waveformID networkCode="XA" stationCode="S12"></waveformID>
  #   <filterID>smi:jpl/filter/Apollo/standard_0.25_0.75</filterID>
  #   <phaseHint>t_max</phaseHint>
  # </pick>

def find_phase_time(origin,time_before,time_in_secs):
    phase_time = UTCDateTime(origin) + time_in_secs - time_before
    print(phase_time)
    return(phase_time)
    

if __name__ == "__main__":
    origin='1971-02-04T07:41:13.500000Z'
    time_before=600
    time_in_secs=200
    find_phase_time(origin,time_before,time_in_secs)
