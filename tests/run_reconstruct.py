#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests to check reconstructing streams. 

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
from pdart.splice import call_splice_chains
from pdart.view import plot_from_file, plot_from_stream
from pdart.chains import build_chains
import gzip
import shutil
import numpy as np
from pdart.reconstruct_frames import _reconstruct_streams, reconstruct_frames

from obspy.core import Stream, Trace, Stats, read

from pdart.raw_import import raw_import


def test_reconstruct_good_stream():
    DELTA = 0.15094
    INVALID = -99999
    stream = read()
    stream = stream[0:2]
    tr1 = stream[0]
    tr1.stats.delta=DELTA
    tr1.stats.channel='_TT'
    tr1.data = np.array([0,1*DELTA, 2*DELTA, 3*DELTA ])

    tr2 = stream[1]
    tr2.stats.delta=0.15094
    tr2.stats.channel='_FR'
    tr2.data = np.array([0,0.25, 0.5, 0.75 ])

    n = tr2.stats.npts
    pointer_array = np.full(n, INVALID, 'int32')
    pointer_array = np.array([0,1,2,3 ])
    re_stream, last_pointer = _reconstruct_streams(stream, pointer_array)

    print(re_stream)
    print(last_pointer)

def test_reconstruct_gap_stream():
    DELTA = 0.15094
    INVALID = -99999
    stream = read()
    stream = stream[0:2]
    tr1 = stream[0]
    tr1.stats.delta=DELTA
    tr1.stats.channel='_TT'
    tr1.data = np.array([0,1*DELTA, 4*DELTA, 5*DELTA,100 ])
    tr1.data += tr1.stats.starttime.timestamp

    tr2 = stream[1]
    tr2.stats.delta=0.15094
    tr2.stats.channel='_FR'
    tr2.data = np.array([0,0.25, 1, 1.25, 1.5])

    n = tr2.stats.npts
    pointer_array = np.full(n, INVALID, 'int32')
    pointer_array = np.array([0,1,4,5, INVALID ])
    re_stream, last_pointer = _reconstruct_streams(stream, pointer_array)

    print(re_stream)
    print(last_pointer)

def test_reconstruct_frames():

    DELTA = 0.15094
    INVALID = -99999
    stream = read()
    stream = stream[0:2]
    tr1 = stream[0]
    tr1.stats.delta=DELTA
    tr1.stats.channel='_TT'
    tr1.data = np.array([0,1*DELTA, 2*DELTA, 3*DELTA ])

    tr2 = stream[1]
    tr2.stats.delta=0.15094
    tr2.stats.channel='_FR'
    tr2.data = np.array([0,0.25, 0.5, 0.75 ])

    return_stream = reconstruct_frames(stream)
    print(return_stream)


def test_reconstruct_gap_frames():

    DELTA = 0.15094
    INVALID = -99999
    stream = read()
    stream = stream[0:2]
    tr1 = stream[0]
    tr1.stats.delta=DELTA
    tr1.stats.channel='_TT'
    tr1.data = np.array([0,1*DELTA, 4*DELTA, 5*DELTA,100 ])
    tr1.data += tr1.stats.starttime.timestamp

    tr2 = stream[1]
    tr2.stats.delta=DELTA
    tr2.stats.channel='_FR'
    tr2.data = np.array([0,0.25, 1, 1.25, 1.5])

    return_stream = reconstruct_frames(stream)
    print(return_stream)

if __name__ == "__main__":
    # test_reconstruct_streams()
    # test_reconstruct_gap_stream()
    # test_reconstruct_frames()
    test_reconstruct_gap_frames()
