#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Splices the records together.

:copyright:
    The pdart Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from datetime import timedelta
import os
import io
import gzip
import glob
import numpy as np
import numpy.ma as ma
import logging, logging.handlers
import csv
import fnmatch

from obspy.core import Stream, Trace, Stats, read

def stream_select(stream, network=None, station=None, location=None, channel=None,
           sampling_rate=None, npts=None, component=None, id=None,
           starttime=None, endtime=None):
    """
    Return new Stream object only with these traces that match the given
    stats criteria (e.g. all traces with ``channel="EHZ"``).

    This is based on the Stream.select() method from ObsPy, but it
    includes starttime and endtime (and is not a class method.)

    .. warning::
        A new Stream object is returned but the traces it contains are
        just aliases to the traces of the original stream. Does not copy
        the data but only passes a reference.

    All keyword arguments except for ``component`` are tested directly
    against the respective entry in the :class:`~obspy.core.trace.Stats`
    dictionary.

    If a string for ``component`` is given (should be a single letter) it
    is tested against the last letter of the ``Trace.stats.channel`` entry.

    Alternatively, ``channel`` may have the last one or two letters
    wildcarded (e.g. ``channel="EH*"``) to select all components with a
    common band/instrument code.

    All other selection criteria that accept strings (network, station,
    location) may also contain Unix style wildcards (``*``, ``?``, ...).
    """
    # make given component letter uppercase (if e.g. "z" is given)
    if component and channel:
        component = component.upper()
        channel = channel.upper()
        if channel[-1] != "*" and component != channel[-1]:
            msg = "Selection criteria for channel and component are " + \
                  "mutually exclusive!"
            raise ValueError(msg)
    new_stream = Stream()
    for trace in stream:
        # skip trace if any given criterion is not matched
        if id and not fnmatch.fnmatch(trace.id.upper(), id.upper()):
            continue
        if network is not None:
            if not fnmatch.fnmatch(trace.stats.network.upper(),
                                   network.upper()):
                continue
        if station is not None:
            if not fnmatch.fnmatch(trace.stats.station.upper(),
                                   station.upper()):
                continue
        if location is not None:
            if not fnmatch.fnmatch(trace.stats.location.upper(),
                                   location.upper()):
                continue
        if channel is not None:
            if not fnmatch.fnmatch(trace.stats.channel.upper(),
                                   channel.upper()):
                continue
        if sampling_rate is not None:
            if float(sampling_rate) != trace.stats.sampling_rate:
                continue
        if npts is not None and int(npts) != trace.stats.npts:
            continue
        if component is not None:
            if len(trace.stats.channel) < 3:
                continue
            if not fnmatch.fnmatch(trace.stats.channel[-1].upper(),
                                   component.upper()):
                continue
        if starttime is not None:
            if starttime != trace.stats.starttime:
                continue
        if endtime is not None:
            if endtime != trace.stats.endtime:
                continue

        new_stream.append(trace)
    return new_stream

def trace_eq(trace1, trace2):
    """
    Based on __eq__() method in obspy.trace for equality between traces.

    Implements rich comparison of Trace objects for "==" operator.

    Traces are the same, if both their data and stats are the same.
    Can also check for equality in masked traces.
    """
    # check if trace1 is a Trace
    if not isinstance(trace1, Trace):
        return False
    # check if trace2 is a Trace
    if not isinstance(trace2, Trace):
        return False
    # comparison of Stats objects is supported by underlying AttribDict
    if not trace1.stats == trace2.stats:
        return False
    # if the array is masked, check for equality of masked array
    if isinstance(trace1.data,np.ma.MaskedArray) and isinstance(trace2.data,np.ma.MaskedArray):
        if not ma.allequal(trace1.data, trace2.data, fill_value=True):
            return False
    # comparison of ndarrays is supported by NumPy
        if not np.array_equal(trace1.data, trace2.data):
            return False

    return True
