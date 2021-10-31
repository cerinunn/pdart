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
from pdart.view import plot_availability, save_availability, save_availability_coarse
from obspy.core import read 


import matplotlib  
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

'''
'''

def call_ALL_save_availability():
    '''
    Plot the availability - takes at least 15 minutes to run
    '''

    # start_time = UTCDateTime(year=1970,julday=1)
    # start_time = UTCDateTime('1969-11-19T14:23:00Z')
    # end_time = UTCDateTime('1977-09-30T21:00:00.000000Z')

    # end_time = UTCDateTime('1969-10-30T21:00:00.000000Z')
    # end_time = start_time + timedelta(days=1)
    # start_time = UTCDateTime('1972-04-22T00:00:00.000000Z')
    # # end_time = UTCDateTime('1972-05-31T00:00:00.000000Z')
    # end_time = UTCDateTime('1972-06-30T00:00:00.000000Z')
    # end_time = UTCDateTime(year=1974,julday=1,hour=1)

    save_availability_coarse(
      # start_time=start_time,
      # end_time=end_time,
    #   file_dir='/Users/nunn/lunar_data/PDART_ELAPSED_TIME',
      # stations=['S11'],
      file_dir='/Users/cnunn/lunar_data/PDART',
      pdf_dir='.',
      read_gzip=False,
      save_pdf=True,
      outfile='../extra_plots_output/new_availability_coarse.csv')

# def call_ALL_save_availability_work_tapes():
#     '''
#     Plot the availability - takes at least 15 minutes to run
#     '''
# 
#     # start_time = UTCDateTime(year=1970,julday=1)
#     start_time = UTCDateTime('1976-01-01')
#     end_time = UTCDateTime('1977-10-01')
#     # end_time = start_time + timedelta(days=1)
#     # start_time = UTCDateTime('1972-04-22T00:00:00.000000Z')
#     # # end_time = UTCDateTime('1972-05-31T00:00:00.000000Z')
#     # end_time = UTCDateTime('1972-06-30T00:00:00.000000Z')
#     # end_time = UTCDateTime(year=1974,julday=1,hour=1)
# 
#     save_availability(
#       start_time=start_time,
#       end_time=end_time,
#     #   file_dir='/Users/nunn/lunar_data/PDART_ELAPSED_TIME',
#       stations=['S12','S14','S15','S16'],
#       # stations=['S16'],
#       file_dir='/Users/cnunn/lunar_data/PDART_CONTINUOUS_WORK_TAPES/',
#       pdf_dir='.',
#       read_gzip=False,
#       save_pdf=True,
#       outfile='availability_work_tapes.csv')


def call_ALL_plot_availability():
    '''
    Plot the availability - takes at least 15 minutes to run
    '''

    # start_time = UTCDateTime(year=1971,julday=1)
    # end_time = start_time + timedelta(days=1)
    # start_time = UTCDateTime('1969-11-01')
    # end_time = UTCDateTime('1977-09-30')

    plot_availability(
      filenames = ['../extra_plots_output/new_availability_coarse.csv'],
      # start_time=start_time,
      # end_time=end_time,
    #   file_dir='/Users/nunn/lunar_data/PDART_ELAPSED_TIME',
      # stations=['S12','S14','S15','S16'],
      # stations=['S16'],
      # channels=['MHZ'],
      pdf_dir='.',
      read_gzip=False,
      outfile='../extra_plots_output/data_availabilityXX.png',
      save_pdf=True)

# def call_ALL_plot_availability():
#     '''
#     Plot the availability - takes at least 15 minutes to run
#     '''
# 
#     # start_time = UTCDateTime(year=1971,julday=1)
#     # end_time = start_time + timedelta(days=1)
#     start_time = UTCDateTime('1969-11-01')
#     end_time = UTCDateTime('1977-09-30')
# 
#     plot_availability(
#       filenames = ['availability_main_tapes.csv','availability_work_tapes.csv'],
#       start_time=start_time,
#       end_time=end_time,
#     #   file_dir='/Users/nunn/lunar_data/PDART_ELAPSED_TIME',
#       # stations=['S12','S14','S15','S16'],
#       # stations=['S16'],
#       # channels=['MHZ'],
#       pdf_dir='.',
#       read_gzip=False,
#       outfile='../extra_plots_output/availX.png',
#       save_pdf=True)

def call_short_plot_availability():
    '''
    Plot the availability - takes at least 15 minutes to run
    '''

    start_time = UTCDateTime(year=1976,julday=293)
    end_time = UTCDateTime(year=1976,julday=295)

    plot_availability(
      filenames = ['tmp_availability_work_tapes.csv'],
      start_time=start_time,
      end_time=end_time,
    #   file_dir='/Users/nunn/lunar_data/PDART_ELAPSED_TIME',
      # stations=['S12','S14','S15','S16'],
      # stations=['S16'],
      # channels=['MHZ'],
      pdf_dir='.',
      read_gzip=False,
      outfile='../extra_plots_output/availX.png',
      save_pdf=True)

def call_short_save_availability_coarse_locations():
    '''
    Save the availability file - a short section which includes locations
    '''

    start_time = UTCDateTime(year=1974,julday=288)
    end_time = UTCDateTime(year=1974,julday=291)


# ./tmp_PDART_CONTINUOUS_WORK_TAPES/s15/1976/294/xa.s15..shz.1976.294.0.mseed
# ./tmp_PDART_CONTINUOUS_WORK_TAPES/s15/1976/293/xa.s15..shz.1976.293.0.mseed


    save_availability_coarse(
      start_time=start_time,
      end_time=end_time,
    #   file_dir='/Users/nunn/lunar_data/PDART_ELAPSED_TIME',
      # stations=['S12','S14','S15','S16'],
      stations=['S12'],
      file_dir='/Users/cnunn/lunar_data/tmp_PDART/',
      pdf_dir='.',
      read_gzip=False,
      save_pdf=True,
      outfile='../extra_plots_output/tmp_availability_coarse_short_locations.csv')

def call_short_save_availability_locations():
    '''
    Save the availability file - a short section which includes locations
    '''

    start_time = UTCDateTime(year=1974,julday=288)
    end_time = UTCDateTime(year=1974,julday=291)


# ./tmp_PDART_CONTINUOUS_WORK_TAPES/s15/1976/294/xa.s15..shz.1976.294.0.mseed
# ./tmp_PDART_CONTINUOUS_WORK_TAPES/s15/1976/293/xa.s15..shz.1976.293.0.mseed


    save_availability(
      start_time=start_time,
      end_time=end_time,
    #   file_dir='/Users/nunn/lunar_data/PDART_ELAPSED_TIME',
      # stations=['S12','S14','S15','S16'],
      stations=['S12'],
      file_dir='/Users/cnunn/lunar_data/tmp_PDART/',
      pdf_dir='.',
      read_gzip=False,
      save_pdf=True,
      outfile='../extra_plots_output/tmp_availability_short_locations.csv')

def call_short_plot_availability_locations():
    '''
    Plot the availability  - a short section which includes locations
    '''

    start_time = UTCDateTime(year=1974,julday=288)
    end_time = UTCDateTime(year=1974,julday=291)

    plot_availability(
      filenames = ['../extra_plots_output/tmp_availability_short_locations.csv'],
      start_time=start_time,
      end_time=end_time,
    #   file_dir='/Users/nunn/lunar_data/PDART_ELAPSED_TIME',
      # stations=['S12','S14','S15','S16'],
      # stations=['S16'],
      # channels=['MHZ'],
      pdf_dir='.',
      read_gzip=False,
      outfile='../extra_plots_output/availXX.png',
      save_pdf=False)

def call_test_all_plot_availability():
    '''
    Plot the availability  - test all time

    This is test data for all time to get the width of the plot right. 
    '''

    plot_availability(
      filenames = ['../extra_plots_output/tmp_all_time.csv'],
      # start_time=start_time,
      # end_time=end_time,
    #   file_dir='/Users/nunn/lunar_data/PDART_ELAPSED_TIME',
      # stations=['S12','S14','S15','S16'],
      # stations=['S16'],
      # channels=['MHZ'],
      pdf_dir='.',
      read_gzip=False,
      outfile='../extra_plots_output/avail_timeXX.png',
      save_pdf=False)

if __name__ == "__main__":

    # call_short_save_availability_work_tapes()
    # call_short_plot_availability()

    # test saving for three days 
    # call_short_save_availability_locations()
    # display the data 
    # call_short_plot_availability_locations()

    # get the width of the plot right with some simple test data 
    # call_test_all_plot_availability()

    # call_ALL_save_availability()
    # rerun the availability files before running this 
    call_ALL_plot_availability()

    # test saving for three days with the coarse method
    # call_short_save_availability_coarse_locations()
