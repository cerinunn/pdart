#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
View the catalog

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)



"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from pdart.diffusion.view_catalog_with_envelopes import view_catalog_with_envelopes
# from builtins import input
# from obspy import read_events, read_inventory
# from obspy.core import read, Stream
# from obspy.core.util.misc import limit_numpy_fft_cache
# from obspy.signal.filter import envelope
# from future.utils import native_str
# import os.path
# import numpy as np
# from matplotlib.dates import date2num
# import math
# 
# 
# import matplotlib
# matplotlib.use('TkAgg')
# import matplotlib.pyplot as plt
# 
# # from matplotlib import rcParams
# # rcParams.update({'figure.autolayout': True})
# 


if __name__ == "__main__":

    top_level_processed_dir = '/Users/nunn/lunar_data/PDART_PROCESSED'
    inv_name = "/Users/nunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    pre_filt_main = [0.2,0.3,0.9,1.3]
    pre_filt_env = [[0.1,0.25,0.75,1],[0.5,0.75,1.25,1.5],[1,1.25,1.75,2.25]]

    view_catalog_with_envelopes(top_level_processed_dir,
      '../lunar_proposal/catalogs/S12_shallow.xml',inv_name,
     dir_type='processed_dir',pre_filt_main=pre_filt_main,
     pre_filt_env=pre_filt_env,output='ACC',smooth_periods=10,
     plot_type='with_seismogram')

