#!/usr/bin/env python

'''
Plot an overview
'''


from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
# from pdart.view import calc_average_timing, plot_overview_from_file, plot_from_file
from pdart.view import plot_overview_from_file, calc_average_sampling_interval


'''
'''


def call_calculate_average_sampling_interval():
    '''
    Calculate the averages - they are saved to a separate directory

    '''

    start_time = UTCDateTime(year=1973,month=1,day=1,hour=0)
    end_time = UTCDateTime(year=1974,month=1,day=1,hour=0)
    stations = ['S12','S14','S15','S16']
    stations = ['S16']
    # stations = ['S14']
    channels = ['ASI'] # the calculated time intervals
    # end_time = start_time + timedelta(days=5)
    # end_time = UTCDateTime(year=1971,month=3,day=2)

    calc_average_sampling_interval(
      stations=stations,
      file_dir='/Users/nunn/lunar_data/PDART',
      av_dir='/Users/nunn/lunar_data/AV_TIMING_PDART',
      start_time=start_time,
      end_time=end_time,
      read_gzip=True,
      write_gzip=True,
    #   sampling_average=23851, # approx 1 hour
      sampling_average=5963 # approx 15 mins
    )


def plot_overview_sampling_interval():
    '''
    Plot an overview of 1973 - the timing interval
    Note - the timing interval is already calculated in the
    files.

    '''

    start_time = UTCDateTime(year=1973,month=1,day=1,hour=0)
    end_time = UTCDateTime(year=1974,month=1,day=1,hour=0)
    # stations = ['S12','S14','S15','S16']
    stations = ['S16']
    # stations = ['S14']
    channels = ['ASI'] # the calculated time intervals
    # end_time = start_time + timedelta(days=5)
    # end_time = UTCDateTime(year=1971,month=3,day=2)

    plot_overview_from_file(
      stations=stations,
      channels=channels,
    #   ylim = (480 - 100, 480 + 100),
    #   ylim = (0.148,0.152),
      file_dir='/Users/nunn/lunar_data/AV_TIMING_PDART',
      pdf_dir='.',
      start_time=start_time,
      end_time=end_time,
      merge_locations=True,
      save_pdf=True)


def plot_close_up_ground_statons():

    plot_overview_from_file(
        stations=['S12'],
        channels=['ASI'],
        file_dir='/Users/nunn/lunar_data/AV_TIMING_PDART',
        ylim=(0.150942,0.150946),
        pdf_dir='.',
        start_time=UTCDateTime('1973-01-01T00:00:00.000000Z'),
        end_time=UTCDateTime('1973-01-08T00:00:00.000000Z'),
        read_gzip=True,
        merge_locations=False,
        save_pdf=True,
        plot_solar_presence=True)

def plot_close_up():

    plot_overview_from_file(
        stations=['S1'],
        channels=['ASI'],
        file_dir='/Users/nunn/lunar_data/AV_TIMING_PDART',
        ylim=(0.150938,0.150948),
        pdf_dir='.',
        start_time=UTCDateTime('1973-01-01T00:00:00.000000Z'),
        end_time=UTCDateTime('1973-01-08T00:00:00.000000Z'),
        read_gzip=True,
        merge_locations=True,
        save_pdf=True,
        plot_solar_presence=False)

def plot_close_up_A14():

    plot_overview_from_file(
        stations=['S14'],
        channels=['ASI'],
        file_dir='/Users/nunn/lunar_data/AV_TIMING_PDART',
        ylim=(0.150932,0.150941),
        pdf_dir='.',
        start_time=UTCDateTime('1973-01-01T00:00:00.000000Z'),
        end_time=UTCDateTime('1973-01-08T00:00:00.000000Z'),
        read_gzip=True,
        merge_locations=False,
        save_pdf=False,
        plot_solar_presence=False)

def plot_temp():

    plot_overview_from_file(
        # this was plotting S14, now S16
        # stations=['S14'],
        stations=['S16'],
        channels=['ASI'],
        file_dir='/Users/nunn/lunar_data/AV_TIMING_PDART',
        # ylim=(0.150928,0.150945),
        # ylim=(0.150928,0.150945),
        pdf_dir='.',
        # start_time=UTCDateTime('1973-05-07T00:00:00.000000Z'),
        # end_time=UTCDateTime('1973-05-15T00:00:00.000000Z'),
        # start_time=UTCDateTime('1973-05-21T00:00:00.000000Z'),
        # end_time=UTCDateTime('1973-05-31T00:00:00.000000Z'),
        start_time=UTCDateTime('1973-05-21T00:00:00.000000Z'),
        end_time=UTCDateTime('1973-05-31T00:00:00.000000Z'),
        ylim=(0.150930,0.1509355),
        read_gzip=True,
        merge_locations=False,
        save_pdf=False,
        plot_solar_presence=False)

if __name__ == "__main__":
    call_calculate_average_sampling_interval()
    plot_overview_sampling_interval()
    # plot_close_up_ground_statons()
    # plot_close_up()
    # plot_close_up_A14()
    # plot_temp()
