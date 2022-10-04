#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilities

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
import pandas as pd

from obspy.core import Stream, Trace, Stats, read
from obspy.core.utcdatetime import UTCDateTime

import matplotlib  
from matplotlib import pyplot as plt

DELTA = 0.1509433962

def timing_correction(stream,correction_time,interpolate=False):
    """Snippet to move the sample starttime so that it agrees with the timestamp 

    The data sampler had small variations in sampling rate, and also variations
    between the different stations. 
    The timestamp
    is recorded on the ATT trace. The data are provided as contiuous 
    traces with a nominal sampling rate. Over the course of the day 
    there will be a divergence between the recorded timestamp and the 
    sample time. This method shifts the traces, so that there is no 
    divergence between timestamp and sample time at the correction time. 
    If a moonquake event
    occurs later in the day, it may be useful to shift each of the data
    traces for each of the stations so that they are using a common reference
    time. The method will shift different stations different amounts. 
    Be aware that the timing trace is not always accurate.
    
    The trace displays sample time minus timestamp and the channel is 
    'ATT_DIVERGENCE'.

    :type stream: :class:`~obspy.core.Stream` 
    :param stream: A data stream
    :type correction_time: class:`~obspy.core.UTCDateTime`
    :param correction_time: The time where there will be no difference
    between the timestamp and the sample time.

    """

    print('need to make corrections for each channel in the trace') 
    # TODO copy some more data to the snippets_files directory for testing 
    # TODO think about whether the trace is in relative mode or not
    print('Warming - merge may not be necessary')
    stream.merge()
    station_set = set()
    for tr in stream:
        station_set.add(tr.stats.station)
    for station in station_set:
        stream_ATT = stream.select(station=station, channel='ATT')

        if len(stream_ATT) == 0:
            print('A timing trace (ATT) is required.')
            raise Exception
        tr_ATT = stream_ATT[0]
    
        sample_no_upper = None
        sample_no_lower = None

        found = np.where(tr_ATT.data > correction_time.timestamp)
        if len(found[0]) > 0:
            sample_no_upper = found[0][0]
        else: 
            sample_no_upper = None

        found_lower = np.where(tr_ATT.data < correction_time.timestamp)
        if len(found_lower[0]) > 0:
            sample_no_lower = found_lower[0][-1]
        else:
            sample_no_lower = None

        if sample_no_lower is None or sample_no_upper is None:
            print('Unable to correct timing at station ', station)
        else: 
            if (sample_no_upper - sample_no_lower) == 1:
                sample_no = sample_no_upper
                # find time at that sample number 
                sample_time = tr_ATT.stats.starttime + tr_ATT.times()[sample_no]
                timestamp_time = UTCDateTime(tr_ATT.data[sample_no])
            else:
                xp = [tr_ATT.data[sample_no_lower], tr_ATT.data[sample_no_upper]]
                fp = [sample_no_lower, sample_no_upper]
                sample_no_interp = np.interp([correction_time.timestamp], xp, fp )

                sample_no = int(sample_no_interp)
                
                print('Warning - Interpolating the sample number for the timing correction.')
                sample_time = tr_ATT.stats.starttime + sample_no * tr_ATT.stats.delta
                timestamp_time = correction_time

            # find difference between sample time and timestamp
            time_diff = timestamp_time - sample_time

            print('Correction Time for Station ', station, time_diff)
            
            # make the correction 
            for tr in stream.select(station=station):
                tr.stats.starttime += time_diff

def relative_timing_stream(stream,remove_original=True):
    """Snippet to change the stream to include relative timing traces  

    The data sampler had small variations in sampling rate, and also variations
    between the different stations. 
    The timestamp
    is recorded on the ATT trace. The data are provided as contiuous 
    traces with a nominal sampling rate. Over the course of the day 
    there will be a divergence between the recorded timestamp and the 
    sample time. The new trace contains sample time minus timestamp and the 
    channel is called 'ATT_DIVERGENCE'. It is also a lot easier to see 
    variations in the timing the trace in this format.

    :type stream: :class:`~obspy.core.Stream` 
    :param stream: A data stream
    :type remove_original: bool
    :param remove_original: If set to ``True``, the original ATT trace is removed.
        Otherwise, a new trace is added to the stream. Defaults to True.
    """
    # remove the negative ones where there are data gaps

    # stream.merge()
    
    for tr in stream:
        if tr.stats.channel == 'ATT' :
            remove_negative_ones_trace(tr)
            if remove_original: 
                tr_divergence = tr
            else:
                tr_divergence = tr.copy()
                stream.append(tr_divergence)
            tr_divergence.stats.channel = 'ATT_DIVERGENCE'
            # get the mask if there is one (when there are data gaps)
            if isinstance(tr_divergence.data, np.ma.MaskedArray):
                mask = np.ma.getmask(tr_divergence.data)
            else:
                mask = None

            # find the relative time
            timestamp0 = tr_divergence.stats.starttime.timestamp
            timestamp_arr = timestamp0 + np.arange(0,len(tr_divergence.data))* DELTA*4
            tr_divergence.data = timestamp_arr - tr.data

            # apply the mask back, if necessary
            if mask is not None:
                tr_divergence.data = ma.array(tr_divergence.data, mask=mask)
            

def relative_timing_trace(trace):
    """Snippet to calculate a relative timing (ATT) trace

    The data sampler had small variations in sampling rate, and also variations
    between the different stations. 
    The timestamp
    is recorded on the ATT trace. The data are provided as contiuous 
    traces with a nominal sampling rate. Over the course of the day 
    there will be a divergence between the recorded timestamp and the 
    sample time. The new trace contains sample time minus timestamp and the 
    channel is called 'ATT_DIVERGENCE'. It is also a lot easier to see 
    variations in the timing the trace in this format.

    :type trace: :class:`~obspy.core.Trace` 
    :param trace: A timing trace
    """
    if trace.stats.channel == 'ATT' :
        remove_negative_ones_trace(trace)
        trace.stats.channel = 'ATT_DIVERGENCE'
        # get the mask if there is one (when there are data gaps)
        if isinstance(trace.data, np.ma.MaskedArray):
            mask = np.ma.getmask(trace.data)
        else:
            mask = None

        # find the relative time
        timestamp0 = trace.stats.starttime.timestamp
        timestamp_arr = timestamp0 + np.arange(0,len(trace.data))* DELTA*4
        trace.data = timestamp_arr - trace.data

        # apply the mask back, if necessary
        if mask is not None:
            trace.data = ma.array(trace.data, mask=mask)    

def remove_negative_ones(stream,channels=['MH1','MH2','MHZ','SHZ','ATT']):
    """Snippet to remove the -1 values in the data traces. 

    The SHZ traces have missing data samples 3-4 times every 32 samples. 
    Providing the seed data with these missing data would mean using very 
    large files. Instead, we provide the data with -1 replacing the gaps. 
    To change the files to include the gaps, use this simple method to 
    replace the -1 values. 
    """

    for tr in stream:
        if tr.stats.channel in channels:
            if tr.stats.channel in ('MH1','MH2','MHZ','SHZ'):
                tr.data = np.ma.masked_where(tr.data==-1, tr.data)
            elif tr.stats.channel in ('ATT'):
                tr.data = np.ma.masked_values(tr.data,-1.0)

def remove_negative_ones_trace(trace):
    """Snippet to remove the -1 values in the data traces. 

    The SHZ traces have missing data samples 3-4 times every 32 samples. 
    Providing the seed data with these missing data would mean using very 
    large files. Instead, we provide the data with -1 replacing the gaps. 
    To change the files to include the gaps, use this simple method to 
    replace the -1 values. 
    """

    if trace.stats.channel in ('MH1','MH2','MHZ','SHZ'):
        trace.data = np.ma.masked_where(trace.data==-1, trace.data)
    elif trace.stats.channel in ('ATT'):
        trace.data = np.ma.masked_values(trace.data,-1.0)

def linear_interpolation(trace,interpolation_limit=1):
    """Snippet to interpolate missing data.  

    The SHZ traces have missing data samples 3-4 times every 32 samples. 
    Providing the seed data with these missing data would mean using very 
    large files. Instead, we provide the data with -1 replacing the gaps. 
    To change the files to interpolate across the gaps, use this simple method to 
    replace the -1 values. The trace is modified, and a mask is applied at 
    the end if necessary. 

    :type stream: :class:`~obspy.core.Trace` 
    :param trace: A data trace
    :type interpolation_limit: int 
    :param interpolation_limit: Limit for interpolation. Defaults to 1. For
      more information read the options for the `~pandas.Series.interpolate`
      method. 

    :return: original_mask :class:`~numpy.ndarray` or class:`~numpy.bool_`
       Returns the original mask, before any interpolation is made. 

    """

    trace.data = np.ma.masked_where(trace.data == -1, trace.data)
    original_mask = np.ma.getmask(trace.data)
    data_series = pd.Series(trace.data)
    # data_series.replace(-1.0, pd.NA, inplace=True)
    data_series.interpolate(method='linear', axis=0, limit=interpolation_limit, inplace=True, limit_direction=None, limit_area='inside', downcast=None)
    data_series.fillna(-1.0, inplace=True)
    trace.data=data_series.to_numpy(dtype=int)
    trace.data = np.ma.masked_where(trace.data == -1, trace.data)
    return original_mask

def maximize_plot(backend=None,fullscreen=False):
    """Maximize window independently of backend.
    Fullscreen sets fullscreen mode, that is same as maximized, but it doesn't have title bar (press key F to toggle full screen mode)."""
    if backend is None:
        backend=matplotlib.get_backend()
    mng = plt.get_current_fig_manager()

    if fullscreen:
        mng.full_screen_toggle()
    else:
        if backend == 'wxAgg':
            mng.frame.Maximize(True)
        elif backend == 'Qt4Agg' or backend == 'Qt5Agg':
            mng.window.showMaximized()
        elif backend == 'TkAgg':
            mng.window.state('zoomed') #works fine on Windows!
        else:
            print ("Unrecognized backend: ",backend) #not tested on different backends (only Qt)
    plt.show()

    plt.pause(0.1) #this is needed to make sure following processing gets applied (e.g. tight_layout)
