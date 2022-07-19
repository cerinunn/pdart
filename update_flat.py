#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Set location to flat when the seismometers were operating in 
flat mode. It retrieves both types of locations on the relevant days
and resets them. 

:copyright:
    The pdart Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)

Update the locations of streams with flat response files from 00 to 01.
There is no problem with rerunning the method - it will read in files 
with either location and write out files with the correct one. 

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from datetime import datetime, timedelta
import os
import io
import gzip
import glob
import math
import numpy as np
import numpy.ma as ma
import logging, logging.handlers
import csv
import fnmatch
import shutil
import pandas as pd
from random import choice, seed
from itertools import groupby
import sys, traceback
from collections import OrderedDict

from obspy.core.utcdatetime import UTCDateTime
from obspy.core import Stream, Trace, Stats, read

from pdart.view import find_filename_date_upper, find_dir_upper


DELTA = 0.1509433962

flat_times = [
['S12', UTCDateTime('1974-10-16T14:02:36.073000Z') , UTCDateTime('1975-04-09T15:31:03.702000Z')],
['S12', UTCDateTime('1975-06-28T13:48:23.124000Z') , UTCDateTime('1977-03-27T15:41:06.247000Z')],
['S14', UTCDateTime('1976-09-18T08:24:35.026000Z') , UTCDateTime('1976-11-17T15:34:34.524000Z')],
['S14', UTCDateTime('1976-05-27T08:24:35.026000Z') , UTCDateTime('1976-11-17T15:34:34.524000Z')],
['S15', UTCDateTime('1971-10-24T20:58:47.248000Z') , UTCDateTime('1971-11-08T00:34:39.747000Z')],
['S15', UTCDateTime('1975-06-28T14:36:33.034000Z') , UTCDateTime('1977-03-27T15:24:05.361000Z')],
['S16', UTCDateTime('1972-05-13T14:08:03.157000Z') , UTCDateTime('1972-05-14T14:47:08.185000Z')],
['S16', UTCDateTime('1975-06-29T02:46:45.610000Z') , UTCDateTime('1977-03-26T14:52:05.483000Z')],
]

PEAKED = '00'
FLAT = '01'

from datetime import date, timedelta

def daterange(start_date, end_date):
    start_date = start_date.date
    end_date = end_date.date
    for n in range(int((end_date - start_date).days)):
        # print(start_date + timedelta(n))
        out_date = UTCDateTime(start_date + timedelta(n)) 
        # print(start_date + timedelta(n), out_date)
        yield out_date



def update_flat(
    top_level_dir='.',
    log_dir='.',
    logging_level=logging.INFO,
    ):

    log_filename = 'flat_times.log'
    log_filename = os.path.join(log_dir,log_filename)
    logging.basicConfig(filename=log_filename, filemode='w', 
      level=logging_level,format='%(message)s')
    print('log file ', log_filename)


    for flat_time in flat_times:
        station = flat_time[0]
        flat_starttime = flat_time[1]
        flat_endtime = flat_time[2]

        logging.info('Flat range {} {}.{} - {}.{}'.format(station, flat_starttime.year, flat_starttime.julday, flat_endtime.year, flat_endtime.julday))    

        logging.info('############## Updating Full Days of Flat Mode')


        for i, single_date in enumerate(daterange(flat_starttime, flat_endtime)):
            if i==0:
                continue


            logging.info('Searching for files: {} {} {}.{}'.format(station, single_date, single_date.year, single_date.julday))

            section=True
            if section: 
                for channel in ['MH1', 'MH2', 'MHZ']:
            
                    # get filenames for whole days
                    filenames = get_filenames(
                        top_level_dir=top_level_dir,
                        start_time=UTCDateTime(year=single_date.year, julday=single_date.julday),
                        stations = [station],
                        channels = [channel],
                        merge_locations=False,
                        end_time=UTCDateTime(year=single_date.year, julday=single_date.julday)+ 24*3600
                    )
                    stream = Stream()
            
                    for file1 in filenames:
                        stream += read(file1)

                    if len(stream) > 0:
                        stream.merge()
                        for tr in stream:
                            tr.stats.location = FLAT

                        new_file1 = get_new_filename(file1,location=FLAT)
                        logging.info('New flat file {} '.format(new_file1))

                        stream = stream.split()

                        # write the new file
                        stream.write(new_file1, 'MSEED')
                        logging.info('Writing file: {}'.format(new_file1))
                
                        # delete any old files
                        for file1 in filenames:
                            if file1 != new_file1:
                                logging.info('Deleting file: {}'.format(file1))
                                os.remove(file1)

        section = True
        if section: 
            logging.info('############## Updating First Day of Flat Mode')
            logging.info('{} {}.{}'.format(station, flat_starttime.year, flat_starttime.julday))

            # ##################
            # update first day
            filenames = get_filenames(
                top_level_dir=top_level_dir,
                start_time=UTCDateTime(year=flat_starttime.year, julday=flat_starttime.julday),
                stations = [station],
                channels = ['ATT'],
                merge_locations=False,
                end_time=UTCDateTime(year=flat_starttime.year, julday=flat_starttime.julday) + 24*3600
            )
            
            stream = Stream()
            for file1 in filenames:
                stream += read(file1)
            stream.merge()
            start_time = None
            
            if len(stream) > 0:

                tr_ATT = stream.select(channel='ATT')[0]
                flat_starttime_timestamp = flat_starttime.timestamp
                start_idx = np.argmax(tr_ATT.data >= flat_starttime_timestamp )
                # irritatingly, start_idx will return 0 if it is not found at all, 
                # so check! 
                if tr_ATT.data[start_idx]  >= flat_starttime_timestamp:
                    start_time = get_index_time(tr_ATT.stats.starttime, start_idx)
                
                if start_time is not None:
                    for channel in ['MH1', 'MH2', 'MHZ']:
                        filenames = get_filenames(
                            top_level_dir=top_level_dir,
                            start_time=UTCDateTime(year=flat_starttime.year, julday=flat_starttime.julday),
                            stations = [station],
                            channels = [channel],
                            merge_locations=False,
                            end_time=UTCDateTime(year=flat_starttime.year, julday=flat_starttime.julday) + 24*3600
                        )
                
                        stream_channel = Stream()
                        for file1 in filenames:
                            stream_channel += read(file1)
                            for tr in stream_channel:
                                tr.stats.location = ''
                            stream_channel.merge()

                        if len(stream_channel) > 0:
                
                            stream_flat = stream_channel.copy()
                            stream_peaked = stream_channel.copy()
                    
                            new_filenames = []

                            stream_flat = stream_flat.trim(starttime=start_time)
                            flat_file1 = None
                            if len(stream_flat) > 0:
                                for tr in stream_flat:
                                    tr.stats.location = FLAT
                                stream_flat = stream_flat.split()
                                flat_file1 = get_new_filename(file1,FLAT)
                                new_filenames.append(flat_file1)
                                stream_flat.write(flat_file1, 'MSEED')
                                logging.info('Writing file: {}'.format(flat_file1))
                                
                    
                            stream_peaked = stream_peaked.trim(endtime=start_time)

                            peaked_file1 = None
                            if len(stream_peaked) > 0:
                                for tr in stream_peaked:
                                    tr.stats.location = PEAKED
                                stream_peaked = stream_peaked.split()
                                peaked_file1 = get_new_filename(file1,PEAKED)
                                new_filenames.append(peaked_file1)
                                stream_peaked.write(peaked_file1, 'MSEED')
                                logging.info('Writing file: {}'.format(peaked_file1))
                    
                            # delete any extra files (if it changes to flat during the 
                            # middle of the day, there won't be any)
                            for file1 in filenames:
                                if file1 not in new_filenames:
                                    logging.info('Deleting file: {}'.format(file1))
                                    os.remove(file1)
        


        section = True
        if section: 
            logging.info('############## Updating Last Day of Flat Mode')
            logging.info('{} {}.{}'.format(station, flat_endtime.year, flat_endtime.julday))

            # ##################
            # update last day
            filenames = get_filenames(
                top_level_dir=top_level_dir,
                start_time=UTCDateTime(year=flat_endtime.year, julday=flat_endtime.julday),
                stations = [station],
                channels = ['ATT'],
                merge_locations=False,
                end_time=UTCDateTime(year=flat_endtime.year, julday=flat_endtime.julday) + 24*3600
            )

            if len(stream) > 0:

                stream = Stream()
                for file1 in filenames:
                    stream += read(file1)
                stream.merge()
                start_time = None
                
                tr_ATT = stream.select(channel='ATT')[0]
                flat_endtime_timestamp = flat_endtime.timestamp
                start_idx = np.argmax(tr_ATT.data >= flat_endtime_timestamp )
                # irritatingly, start_idx will return 0 if it is not found at all, 
                # so check! 
                if tr_ATT.data[start_idx]  >= flat_endtime_timestamp:
                    start_time = get_index_time(tr_ATT.stats.starttime, start_idx)
                
                if start_time is not None:
                    for channel in ['MH1', 'MH2', 'MHZ']:
                        filenames = get_filenames(
                            top_level_dir=top_level_dir,
                            start_time=UTCDateTime(year=flat_endtime.year, julday=flat_endtime.julday),
                            stations = [station],
                            channels = [channel],
                            merge_locations=False,
                            end_time=UTCDateTime(year=flat_endtime.year, julday=flat_endtime.julday) + 24*3600
                        )
                
                        stream_channel = Stream()
                        for file1 in filenames:
                            stream_channel += read(file1)
                            for tr in stream_channel:
                                tr.stats.location = ''
                            stream_channel.merge()
                
                        stream_flat = stream_channel.copy()
                        stream_peaked = stream_channel.copy()
                
                        new_filenames = []


                
                        stream_peaked = stream_peaked.trim(starttime=start_time)
                        peaked_file1 = None
                        if len(stream_peaked) > 0:
                            for tr in stream_peaked:
                                tr.stats.location = PEAKED
                            stream_peaked = stream_peaked.split()
                            peaked_file1 = get_new_filename(file1,PEAKED)
                            new_filenames.append(peaked_file1)
                            stream_peaked.write(peaked_file1, 'MSEED')
                            logging.info('Writing file: {}'.format(peaked_file1))

                        stream_flat = stream_flat.trim(endtime=start_time)
                        flat_file1 = None
                        if len(stream_flat) > 0:
                            for tr in stream_flat:
                                tr.stats.location = FLAT
                            stream_flat = stream_flat.split()
                            flat_file1 = get_new_filename(file1,FLAT)
                            new_filenames.append(flat_file1)
                            stream_flat.write(flat_file1, 'MSEED')
                            logging.info('Writing file: {}'.format(flat_file1))
                            
                
                        # delete any extra files (if it changes to peaked during the 
                        # middle of the day, there won't be any)
                        for file1 in filenames:
                            if file1 not in new_filenames:
                                logging.info('Deleting file: {}'.format(file1))
                                os.remove(file1)



        

def get_index_time(starttime, idx):
    return starttime + idx * DELTA*4 






def get_new_filename(file1,location):
    # xa.s12.*.mh1.1970.289.*.mseed
    filename = os.path.basename(file1)
    split_filename = filename.split('.')
    if split_filename[0].lower() == 'xa' and split_filename[2] in ('*', PEAKED, FLAT):
        split_filename[2] = location
    new_filename = '.'.join(split_filename)
    new_file1 = os.path.join(os.path.dirname(file1),new_filename)
    return new_file1

def get_filenames(
    top_level_dir,
    start_time,
    stations,
    channels,
    merge_locations=False,
    end_time=None):

    # helper method to read files from the directory
    # and return a stream

    if end_time is None:
        end_time = start_time + timedelta(hours=3)

    time_interval = timedelta(hours=24)

    filenames = []

    start = start_time
    while start < end_time:
        for station in stations:
            directory = find_dir_upper(top_level_dir,station,start)
            for channel in channels:
                filename = find_filename_date_upper(station, '*', channel, start)
                # filename = "{}.gz".format(filename)


                # print(directory)
                filename1 = os.path.join(directory,filename)     
                # this makes sure the file is found, and gets the name without wildcards 
                for f1 in glob.glob(filename1):
                    filenames.append(f1)


        # increment the time interval
        start += time_interval

    return filenames


if __name__ == "__main__":
    top_level_dir = '/Users/cnunn/lunar_data/PDART_V2'
    log_dir='/Users/cnunn/lunar_data/PDART_PROCESSED'
    update_flat(top_level_dir=top_level_dir,log_dir=log_dir)