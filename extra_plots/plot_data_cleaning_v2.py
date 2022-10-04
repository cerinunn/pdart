#!/usr/bin/env python

'''
Plot data cleaning - used in the paper.

'''

from __future__ import print_function
import numpy as np
import os

import matplotlib
from matplotlib import pyplot as plt
# 
# from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
# from matplotlib import gridspec
# from pdart.view import stream_from_directory_new
# from obspy.imaging.util import (_set_xaxis_obspy_dates, _id_key, _timestring)
# from matplotlib.dates import date2num
# from pdart.util import relative_timing_trace, relative_timing_stream, remove_negative_ones
# from matplotlib.dates import AutoDateFormatter, AutoDateLocator, HOURLY
from obspy.core import Stream, read
# import pandas as pd
# import random
# # from pdart.csv_join_work_tapes import stream_import, initial_cleanup, merge_channel_stream, despike3, loose_frame_diff, read_file, process_list
# from pdart.csv_join_work_tapes import despike3, to_Int64, stream_import
# from pdart.csv_check_work_tapes import calculate_gaps, add_or_minus_frame
import pdart.config as config
import pdart.auth as auth
from obspy.clients.fdsn.client import Client
from pdart.util import linear_interpolation




def plot_data_cleaning_v2():
    onset = UTCDateTime('1971-02-07T00:45:25.79')
    # print(UTCDateTime('1971-02-07T00:45:25.79').julday)
    # exit()
    stream = read('/Users/cnunn/lunar_data/GEOSCOPE_lunar_data/1971/S12/MHZ/S12.XA..MHZ.1971.38')
    for tr in stream:
        tr.data = tr.data + 512
        tr.stats.starttime = tr.stats.starttime + 0.23

# 00:34 01:50 

    stream += read('/Users/cnunn/lunar_data/PDART_CORRECTED/s12/1971/038/XA.S12.00.MHZ.1971.038.0.MSEED')
    # stream += read('/Users/cnunn/lunar_data/GEOSCOPE_lunar_data/1976/S12/MH1/S12.XA..MH1.1976.64')
    stream.trim(starttime=onset-300,endtime=onset+3600)
    # stream.trim(starttime= UTCDateTime('1976-12-04T10:00:00'),endtime = UTCDateTime('1976-12-04T11:15:00'))
    for tr in stream:
        # interpolate across the gaps of one sample 
        linear_interpolation(tr,interpolation_limit=1)
    stream.merge()
    # stream.plot(equal_scale=False,size=(1000,600),method='full')

    fig = stream.plot(handle=True, show=False,size=(1200,600), method='full')
    for i, ax in enumerate(fig.get_axes()):
        # if i == 0:
            # ax.set_ylim(512-60,512+60)
        if i == 1:
            ax.set_ylim(512-70,512+70)

    plt.show()



if __name__ == "__main__":
    plot_data_cleaning_v2()