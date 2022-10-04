#!/usr/bin/env python

'''
Plot an overview
'''


from __future__ import print_function
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import plot_overview_from_file


'''
# for best results comment modify
# obspy.imaging.waveform.py
# this will plot the data as dots rather than a time series,
# which is usually easier to see what is going on

Comment this line:
ax.plot(x_values, trace.data, color=self.color,
        linewidth=self.linewidth, linestyle=self.linestyle)
and replace with this:
ax.plot(x_values, trace.data, color=self.color,
         linestyle='None', marker='o', markersize=1)

And comment this line:
ax.plot(x_values, y_values, color=self.color)
and replace with this:
ax.plot(x_values, y_values, color=self.color, marker='o',
    linestyle='None',markersize=1)

'''

def plot_overview_channels():
    '''
    Plot an overview of 1973 and MHZ

    For best results modify obspy.imaging.waveform.py - see note above.
    '''

    # start_time = UTCDateTime(year=1971,month=3,day=1,hour=3)
    # end_time = UTCDateTime(year=1971,month=3,day=7,hour=3)
    start_time = UTCDateTime(year=1973,month=1,day=1,hour=0)
    end_time = UTCDateTime(year=1974,month=1,day=1,hour=0)
    print(start_time.julday)
    stations = ['S12','S14','S15','S16']
    # stations = ['S14']
    channels = ['MHZ'] # the calculated time intervals
    # end_time = start_time + timedelta(days=5)
    # end_time = UTCDateTime(year=1971,month=3,day=2)

    plot_overview_from_file(
      stations=stations,
      channels=channels,
    #   ylim = (480 - 100, 480 + 100),
    #   ylim = (0.148,0.152),
      file_dir='/Users/nunn/lunar_data/PDART',
      pdf_dir='.',
      start_time=start_time,
      end_time=end_time,
      merge_locations=True,
      save_pdf=True)

def plot_march_S15():
    '''
    Plot an overview of 1973 and MHZ

    For best results modify obspy.imaging.waveform.py - see note above.
    '''

    # start_time = UTCDateTime(year=1971,month=3,day=1,hour=3)
    # end_time = UTCDateTime(year=1971,month=3,day=7,hour=3)
    start_time = UTCDateTime(year=1973,month=3,day=8,hour=0)
    end_time = UTCDateTime(year=1973,month=3,day=22,hour=0)
    print(start_time.julday)
    stations = ['S15',]
    # stations = ['S14']
    channels = ['MHZ'] # the calculated time intervals
    # end_time = start_time + timedelta(days=5)
    # end_time = UTCDateTime(year=1971,month=3,day=2)

    plot_overview_from_file(
      stations=stations,
      channels=channels,
    #   ylim = (480 - 100, 480 + 100),
    #   ylim = (0.148,0.152),
      file_dir='/Users/nunn/lunar_data/PDART',
      pdf_dir='.',
      start_time=start_time,
      end_time=end_time,
      merge_locations=True,
      save_pdf=False)

if __name__ == "__main__":
    # plot_overview_channels()
    plot_march_S15()
