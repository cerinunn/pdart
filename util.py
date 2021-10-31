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

from obspy.core import Stream, Trace, Stats, read
from obspy.core.utcdatetime import UTCDateTime

import matplotlib  
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

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
