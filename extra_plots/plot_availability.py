#!/usr/bin/env python

'''
Plot an overview of the availability. 

Start by building a text file which gets a coarse view of the 
availability of each track. Then run the plot. 
'''


from __future__ import print_function
import numpy as np
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import plot_availability, save_availability, save_availability_coarse
from obspy.core import read 


import matplotlib  
# matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

'''
'''

def call_ALL_save_availability():
    '''
    Plot the availability - takes at least 15 minutes to run
    '''

    save_availability_coarse(
      file_dir='/Users/cnunn/lunar_data/PDART_V2',
      pdf_dir='.',
      read_gzip=False,
      save_pdf=True,
      outfile='../extra_plots_output/new_availability_coarseXX.csv')


def call_ALL_plot_availability():
    '''
    Plot the availability - takes at least 15 minutes to run
    '''

    plot_availability(
      filenames = ['../extra_plots_output/new_availability_coarseXX.csv'],
      pdf_dir='.',
      read_gzip=False,
      outfile='../extra_plots_output/data_availabilityXX.png',
      save_pdf=True)

def call_short_plot_availability():
    '''
    Plot the availability - short test
    '''

    start_time = UTCDateTime(year=1976,julday=293)
    end_time = UTCDateTime(year=1976,julday=295)

    plot_availability(
      filenames = ['tmp_availability_work_tapes.csv'],
      start_time=start_time,
      end_time=end_time,
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


    save_availability_coarse(
      start_time=start_time,
      end_time=end_time,
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


    save_availability(
      start_time=start_time,
      end_time=end_time,
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
      pdf_dir='.',
      read_gzip=False,
      outfile='../extra_plots_output/avail_timeXX.png',
      save_pdf=False)

if __name__ == "__main__":

    # call_short_plot_availability()

    # test saving for three days 
    # call_short_save_availability_locations()
    # display the data 
    # call_short_plot_availability_locations()

    # get the width of the plot right with some simple test data 
    # call_test_all_plot_availability()

    call_ALL_save_availability()
    # rerun the availability files before running this 
    call_ALL_plot_availability()

    # test saving for three days with the coarse method
    # call_short_save_availability_coarse_locations()
