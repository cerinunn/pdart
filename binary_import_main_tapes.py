#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Import the data from the binary tape files.
Be careful to make sure that every part of the record is read - otherwise
the file positions will be out for the next record.


:copyright:
    The pdart Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA

from bitstring import BitArray, BitStream, Bits, ConstBitStream, ReadError
from obspy.core.utcdatetime import UTCDateTime
import csv
# from datetime import timedelta
import os
# import io
import gzip
import glob
import pandas as pd
import logging, logging.handlers

station_dict =	{
  1: 'S12',
  2: 'S15',
  3: 'S16',
  4: 'S14',
  5: 'S17'}

# import glob
# import numpy as np
# import numpy.ma as ma

# import csv
# import fnmatch
# import shutil

# from obspy.core.utcdatetime import UTCDateTime
# from obspy.core import Stream, Trace, Stats, read
#
# global DELTA
# DELTA = 0.15094
# INVALID = -99999


# /Users/nunn/lunar_data/JAXA_CSV :grep ",15," 1971-07-31T18:00:00.csv | head
# 8,15,9,1,1971-07-31 18:52:00.269000,1023,1023,1023
# 8,15,9,2,1971-07-31 18:52:00.420000,1023,1023,1023
# 8,15,9,3,1971-07-31 18:52:00.571000,1023,1023,1023
# 8,15,9,4,1971-07-31 18:52:00.722000,1023,1023,1023
# 9,15,9,1,1971-07-31 18:52:00.873000,1023,1023,1023
# 9,15,9,2,1971-07-31 18:52:01.024000,1023,1023,1023
# 9,15,9,3,1971-07-31 18:52:01.175000,1023,1023,1023
# 9,15,9,4,1971-07-31 18:52:01.326000,1023,1023,1023
# 10,15,9,1,1971-07-31 18:52:01.476000,1023,1023,1023
# 10,15,9,2,1971-07-31 18:52:01.627000,1023,1023,1023


# 1971-08-01T18:52:00.269000Z 8 9
# 1023 1023 1023
# 1023 1023 1023
# 1023 1023 1023
# 1023 1023 1023

# def call_raw_import(
#     csv_dir='.',
#     raw_dir='.',
#     start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
#     end_time=UTCDateTime('1977-09-30T:21:00.000000Z'),
#     write_gzip=True):
#     '''
#     Makes 'raw' stream files from gzipped csv files.
#     It works for all stations and all channels.
#     Calls raw_import()
#     '''
#
#     time_interval = timedelta(hours=3)
#     start = start_time
#
#     while start < end_time:
#
#         # find csv gzip importfile  with the start time
#         csv_file_gzip = '%s.csv.gz' % start.strftime("%Y-%m-%dT%H:%M:%S")
#         csv_file_gzip = os.path.join(csv_dir, csv_file_gzip)
#
#         # import the stream from the gzipped file
#         stream = raw_import(csv_file_gzip)
#
#         if len(stream) > 0:
#
#             stations = []
#             for tr in stream:
#                 stations.append(tr.stats.station)
#             stations = list(set(stations))
#
#             for station in stations:
#                 station_stream = stream.select(station=station)
#
#                 channels = []
#                 for tr in station_stream:
#                     channels.append(tr.stats.channel)
#                 channels = list(set(channels))
#
#                 for channel in channels:
#                     channel_stream = station_stream.select(channel=channel)
'''
Long Word 	Bits 	Data
1	0 	software time flag (set at JSC for clock substitution)
1	1-31 	31msb of 35-bit time of the year in ms
2	0-3 	4lsb of 35-bit time of the year in ms
2	4-7 	ALSEP tracking station ID
2	8-13 	bit error rate (from original tape)
2	14 	data rate (1=1060bps, 0=530bps)
2	15-21 	not used
For Old Format
2	22-31	ALSEP word 5
3-18	0-9	ALSEP words 1,4,9,12,16,22,26,29,33,36,40,43,46,52,58,61
3-18	11-20	ALSEP words 2,6,10,13,18,24,27,30,34,37,41,44,48,54,59,62
3-18	22-31	ALSEP words 3,8,11,14,20,25,28,32,35,38,42,45,50,57,60,64
3-18	10,21	not used
For New Format
2	22-31	not used
3-8	0-9	ALSEP words 1,9,25,33,41,46
3-8	11-20	ALSEP words 2,11,27,35,43,57
3-8	22-31	ALSEP words 3,13,29,37,45,59
3-8	10,21	not used
9	0-9	ALSEP word 61
9	10-31	not used
'''

def year_corrections(filename,year):
    '''
    The following files have an incorrect year specified
    '''
    file_list = [
        'pse.a12.10.91.gz',
        'pse.a12.10.92.gz',
        'pse.a12.10.93.gz',
        'pse.a12.10.94.gz',
        'pse.a12.10.95.gz',
        'pse.a12.10.97.gz',
        'pse.a12.10.98.gz',
        'pse.a12.10.102.gz',
        'pse.a12.10.103.gz',
        'pse.a12.10.104.gz',
        'pse.a12.10.106.gz',
        'pse.a12.10.107.gz',
        'pse.a12.10.108.gz',
        'pse.a12.10.109.gz',
        'pse.a12.10.111.gz']

    if os.path.basename(filename) in file_list:
        year = 1976
        print('Year corrected for ', filename, year)

    file_list = [
        'pse.a14.8.112.gz']

    if os.path.basename(filename) in file_list:
        year = 1974
        print('Year corrected for ', filename, year)

    return year

def binary_import_main_tapes(filename, out_filename_gz,
    log_filename='logs/binary.log',
    logging_level=logging.INFO):
    '''
    Long Word 	Bits 	Data
    1	0 	software time flag (set at JSC for clock substitution)
    1	1-31 	31msb of 35-bit time of the year in ms
    2	0-3 	4lsb of 35-bit time of the year in ms
    2	4-7 	ALSEP tracking station ID
    2	8-13 	bit error rate (from original tape)
    2	14 	data rate (1=1060bps, 0=530bps)
    2	15-21 	not used
    For Old Format
    2	22-31	ALSEP word 5
    3-18	0-9	ALSEP words 1,4,9,12,16,22,26,29,33,36,40,43,46,52,58,61
    3-18	11-20	ALSEP words 2,6,10,13,18,24,27,30,34,37,41,44,48,54,59,62
    3-18	22-31	ALSEP words 3,8,11,14,20,25,28,32,35,38,42,45,50,57,60,64
    3-18	10,21	not used
    For New Format
    2	22-31	not used
    3-8	0-9	ALSEP words 1,9,25,33,41,46
    3-8	11-20	ALSEP words 2,11,27,35,43,57
    3-8	22-31	ALSEP words 3,13,29,37,45,59
    3-8	10,21	not used
    9	0-9	ALSEP word 61
    9	10-31	not used
    '''

    logging.basicConfig(filename=log_filename, filemode='w', level=logging_level)

    # read the data as gzipped data
    in_file = gzip.GzipFile(filename, 'rb')
    gzipped_data = in_file.read()
    in_file.close()

    # Note difference - try the 0 in the ground station
    ground_stations_all = set(range(0,15))
    ground_stations = set()
    year = None

    # # open the output file
    # with gzip.open(out_filename_gz, 'wt') as f:
    #     writer = csv.writer(f)
    #     writer.writerow(HEADER)

    # open the BitStream from the file, into fixed-length chunks
    # of 72*8*270+(16*8) bits

    dict_list = []
    count_ground_station0 = 0

    chunks = ConstBitStream(gzipped_data).cut(155648)
    # if you need the files unzipped
    # chunks = ConstBitStream(filename=filename).cut(155648)

    # # for testing, view just a few chunks
    # import itertools
    # chunks = itertools.islice(chunks, 10)

    for s in chunks:

        # read the header
        tape = s.read(16).int
        station = s.read(16).int
        # tape_seq =  s.read(16).int
        # 5-6 original tape sequence number for PSE tapes; 2-digit station code plus 3-digit
        s.read(16).int
        record_no =  s.read(16).int
        year1 =  s.read(16).int
        # check this at the beginning of the file
        if year is None:
            # some years are not specified correctly
            year = year_corrections(filename,year1)
        new_format =  s.read(16).int
        # if new_format == 1:
        #     print('Using new format', filename)
            # return []
        phys_rec =  s.read(16).int
        read_flag = s.read(1).int
        extra_read_flags =  s.read(15).int

        # print(tape, station, record_no, year, new_format, phys_rec, read_flags)
# 1 15 10844 1 1971 0 3 0
        # s.pos = 128 bits

        # print(station, year)
        # 15 1971
        # print(tape, tape_seq, record_no, new_format, phys_rec, read_flags)
        # 1 10844 1 0 3 0

        try:
            if new_format == 1:
                while True:
                    # read in long words 1 and 2
                    clock_flag = s.read(1).uint
                    time = s.read(35).uint

                    # the time is miliseconds after the beginning of the year
                    # but with the first record at 86400 seconds
                    utcDateTime = UTCDateTime(year=year,julday=1)+time/1000-86400.

                    ground_station = s.read(4).uint
                    bit_error_rate= s.read(6).uint
                    data_rate = s.read(1).uint

                    if time != 0:
                        ground_stations.add(ground_station)

                    if ground_station not in ground_stations_all:
                        ground_station = pd.NA

                    s.read(7) # unused
                    s.read(10) # unused
                    # s.pos = 192 bits

                    word1 = s.read(10) # word1
                    s.read(1)  # not used
                    word2 = s.read(10) # word2
                    s.read(1)  # not used
                    word3 = s.read(10) # word3

                    sync = word1 + word2 + word3[0:2]
                    sync = sync.bin
                    frame = word3[2:9].uint
                    mode_bit = word3[9:10].uint

                    mh1_1 = s.read(10).uint # word9
                    s.read(1)  # not used
                    mh2_1 = s.read(10).uint # word11
                    s.read(1)  # not used
                    mhz_1 = s.read(10).uint # word13

                    mh1_2 = s.read(10).uint # word25
                    s.read(1)  # not used
                    mh2_2 = s.read(10).uint # word27
                    s.read(1)  # not used
                    mhz_2 = s.read(10).uint # word29

                    s.read(10) # word33
                    s.read(1)  # not used
                    s.read(10).uint # word35
                    s.read(1)  # not used
                    s.read(10) # word37

                    mh1_3 = s.read(10).uint # word41
                    s.read(1)  # not used
                    mh2_3 = s.read(10).uint # word43
                    s.read(1)  # not used
                    mhz_3 = s.read(10).uint # word45

                    word_46 = s.read(10)
                    # command1 = s.read(7).uint # word46, command (bits 3 - 9)
                    # map1 = s.read(1).uint # word46, map (bit 10)
                    s.read(1)  # not used
                    mh1_4 = s.read(10).uint # word57
                    s.read(1)  # not used
                    mh2_4 = s.read(10).uint # word59

                    mhz_4 = s.read(10).uint # word61
                    s.read(22)  # not used

                    # row = [record_no,frame,str(utcDateTime),clock_flag,frame,str(utcDateTime),station,ground_station,mh1_1,mh2_1,mhz_1,mh1_2,mh2_2,mhz_2, mh1_3,mh2_3,mhz_3, mh1_4,mh2_4,mhz_4]
                    # # write the record to the file if it is valid
                    # # if time != 0:
                    # #     writer.writerow(row)
                    # print('new')
                    # print(row)

                    # word46 used for SHZ at S14, and for command and 
                    # map at the other stations
                    if station == 14:
                        command1 = word5[2:9].uint
                        map1 = word5[9:10].uint
                    else:
                        command1 = word_46[2:9].uint
                        map1 = word_46[9:10].uint
    

                    dict = {
                      'orig_no' : record_no,
                      'orig_timestamp' : str(utcDateTime),
                      'orig_frame' : 0, # no longer used
                      'orig_station' : station,
                      'clock_flag' : clock_flag,
                      'orig_ground_station' : ground_station,
                      'orig_mh1_1' : mh1_1,
                      'orig_mh2_1' : mh2_1,
                      'orig_mhz_1' : mhz_1,
                      'orig_mh1_2': mh1_2,
                      'orig_mh2_2' : mh2_2,
                      'orig_mhz_2' : mhz_2, 
                      'orig_mh1_3' : mh1_3,
                      'orig_mh2_3' : mh2_3,
                      'orig_mhz_3' : mhz_3, 
                      'orig_mh1_4' : mh1_4,
                      'orig_mh2_4' : mh2_4,
                      'orig_mhz_4' : mhz_4,

                      # 'bit_synchronizer' : bit_synchronizer,
                      'frame' : frame,
                      'sync' : sync,
                      'bit_synchronizer' : read_flag, 

                      'command' : command1,
                      'map' : map1,
                      'data_rate' : data_rate
                    }

                    dict_list.append(dict)

                    # break when the physical record has been read
                    if s.pos == 155648:
                        break
            else: # old format
                while True:
                    # read in long words 1 and 2
                    clock_flag = s.read(1).uint
                    time = s.read(35).uint

                    # the time is miliseconds after the beginning of the year
                    # but with the first record at 86400 seconds
                    utcDateTime = UTCDateTime(year=year,julday=1)+time/1000-86400

                    ground_station = s.read(4).uint
                    bit_error_rate= s.read(6).uint
                    data_rate = s.read(1).uint

                    if time != 0:
                        ground_stations.add(ground_station)

                    # print(time, clock_flag, ground_station)

                    s.read(7) # unused
                    word5 = s.read(10) # word5
                    # s.pos = 192 bits

                    word1 = s.read(10) # word1
                    s.read(1)  # not used
                    word2 = s.read(10) # word2
                    s.read(1)  # not used
                    word3 = s.read(10) # word3

                    sync = word1 + word2 + word3[0:2]
                    sync = sync.bin
                    frame = word3[2:9].uint
                    mode_bit = word3[9:10].uint

                    # SHZ all even except 2, 46, 56*

                    shz_4 = s.read(10).uint # word4
                    # word 5 is already read above
                    s.read(1)  # not used
                    shz_6 = s.read(10).uint # word6
                    s.read(1)  # not used
                    shz_8 = s.read(10).uint # word8

                    mh1_1 = s.read(10).uint # word9
                    s.read(1)  # not used
                    shz_10 = s.read(10).uint # word10
                    s.read(1)  # not used
                    mh2_1 = s.read(10).uint # word11

                    shz_12 = s.read(10).uint # word12
                    s.read(1)  # not used
                    mhz_1 = s.read(10).uint # word13
                    s.read(1)  # not used
                    shz_14 = s.read(10).uint # word14

                    shz_16 = s.read(10).uint # word16
                    s.read(1)  # not used
                    shz_18 = s.read(10).uint # word18
                    s.read(1)  # not used
                    shz_20 = s.read(10).uint # word20

                    shz_22 = s.read(10).uint # word22
                    s.read(1)  # not used
                    shz_24 = s.read(10).uint # word24
                    s.read(1)  # not used
                    mh1_2 = s.read(10).uint # word25

                    shz_26 = s.read(10).uint # word26
                    s.read(1)  # not used
                    mh2_2 = s.read(10).uint # word27
                    s.read(1)  # not used
                    shz_28 = s.read(10).uint # word28

                    mhz_2 = s.read(10).uint # word29
                    s.read(1)  # not used
                    shz_30 = s.read(10).uint # word30
                    s.read(1)  # not used
                    shz_32 = s.read(10).uint # word32

                    s.read(10) # word33
                    s.read(1)  # not used
                    shz_34 = s.read(10).uint # word34
                    s.read(1)  # not used
                    s.read(10) # word35

                    shz_36 = s.read(10).uint # word36
                    s.read(1)  # not used
                    s.read(10) # word37
                    s.read(1)  # not used
                    shz_38 = s.read(10).uint # word38

                    shz_40 = s.read(10).uint # word40
                    s.read(1)  # not used
                    mh1_3 = s.read(10).uint # word41
                    s.read(1)  # not used
                    shz_42 = s.read(10).uint # word42

                    mh2_3 = s.read(10).uint # word43
                    s.read(1)  # not used
                    shz_44 = s.read(10).uint # word44
                    s.read(1)  # not used
                    mhz_3 = s.read(10).uint # word45

                    # s.read(2) # word46 first 2 bit
                    # command1 = s.read(7).uint # word46, command
                    # map1 = s.read(1).uint # word46, map
                    # word46 used for SHZ at S14, and for command and 
                    # map at the other stations
                    word_46 = s.read(10)
                    s.read(1)  # not used
                    shz_48 = s.read(10).uint # word48
                    s.read(1)  # not used
                    shz_50 = s.read(10).uint # word50



                    shz_52 = s.read(10).uint # word52
                    s.read(1)  # not used
                    shz_54 = s.read(10).uint # word54
                    s.read(1)  # not used
                    mh1_4 = s.read(10).uint # word57

                    shz_58 = s.read(10).uint # word58
                    s.read(1)  # not used
                    mh2_4 = s.read(10).uint # word59
                    s.read(1)  # not used
                    shz_60 = s.read(10).uint # word60

                    mhz_4 = s.read(10).uint # word61
                    s.read(1)  # not used
                    shz_62 = s.read(10).uint # word62
                    s.read(1)  # not used
                    shz_64 = s.read(10).uint # word64


                    # word46 used for SHZ at S14, and for command and 
                    # map at the other stations
                    if station == 14:
                        shz_46 = word_46.uint
                        command1 = word5[2:9].uint
                        map1 = word5[9:10].uint
                    else:
                        command1 = word_46[2:9].uint
                        map1 = word_46[9:10].uint
                        shz_46 = None

                    if station == 15:
                        # *also excepting ALSEP word 24 for station 15
                        shz_24 = None
    
                    # 16 * 8 = 128
                    # 18* 32 bits (long words) = 576
                    # = 704 bits

                    # row = [record_no,frame,str(utcDateTime),clock_flag,frame,str(utcDateTime),station,ground_station,mh1_1,mh2_1,mhz_1,mh1_2,mh2_2,mhz_2, mh1_3,mh2_3,mhz_3, mh1_4,mh2_4,mhz_4]
                    # write the record to the file if it is valid
                    # if time != 0:
                    #     writer.writerow(row)
                    # else:
                    #     print(s.pos, row)

                    # print('old')
                    # print(row)

                    dict = {
                      'orig_no' : record_no,
                      'orig_timestamp' : str(utcDateTime),
                      'orig_frame' : 0, # no longer used
                      'orig_station' : station,
                      'clock_flag' : clock_flag,
                      'orig_ground_station' : ground_station,
                      'orig_mh1_1' : mh1_1,
                      'orig_mh2_1' : mh2_1,
                      'orig_mhz_1' : mhz_1,
                      'orig_mh1_2': mh1_2,
                      'orig_mh2_2' : mh2_2,
                      'orig_mhz_2' : mhz_2, 
                      'orig_mh1_3' : mh1_3,
                      'orig_mh2_3' : mh2_3,
                      'orig_mhz_3' : mhz_3, 
                      'orig_mh1_4' : mh1_4,
                      'orig_mh2_4' : mh2_4,
                      'orig_mhz_4' : mhz_4, 

                    # SHZ - all even except 2, 46, 56
                    # (also excepting ALSEP word 24 for station 15)
                      'shz_4' : shz_4, 
                      'shz_6' : shz_6, 
                      'shz_8' : shz_8, 

                      'shz_10' : shz_10, 
                      'shz_12' : shz_12, 
                      'shz_14' : shz_14, 
                      'shz_16' : shz_16, 
                      'shz_18' : shz_18, 

                      'shz_20' : shz_20,
                      'shz_22' : shz_22, 
                      'shz_24' : shz_24, 
                      'shz_26' : shz_26,  
                      'shz_28' : shz_28, 

                      'shz_30' : shz_30,
                      'shz_32' : shz_32, 
                      'shz_34' : shz_34, 
                      'shz_36' : shz_36,  
                      'shz_38' : shz_38, 

                      'shz_40' : shz_40,
                      'shz_42' : shz_42, 
                      'shz_44' : shz_44, 
                      'shz_46' : shz_46, 
                      'shz_48' : shz_48, 

                      'shz_50' : shz_50,
                      'shz_52' : shz_52, 
                      'shz_54' : shz_54, 

                      'shz_58' : shz_58, 

                      'shz_60' : shz_60,
                      'shz_62' : shz_62, 
                      'shz_64' : shz_64, 
                      # 'bit_synchronizer' : bit_synchronizer,
                      'frame' : frame,
                      'sync' : sync,
                      'bit_synchronizer' : read_flag,

                      'command' : command1,
                      'map' : map1,
                      'data_rate' : data_rate

                    }

                    dict_list.append(dict)



                    # break when the physical record has been read
                    if s.pos == 155648:
                        break
        except ReadError:
            print('Position of the last line: ', s.pos)

    df = pd.DataFrame.from_dict(dict_list)
    logging.info(df.to_string())

    # initial cleanup
    df['orig_no'] = df['orig_no'].astype('Int64')
    df['orig_timestamp'] = df['orig_timestamp'].astype('datetime64[ns, UTC]')
    df['orig_frame'] = df['orig_frame'].astype('Int64')
    df['orig_station'] = df['orig_station'].astype('string')
    df['clock_flag'] = df['clock_flag'].astype('Int64')
    df['orig_ground_station'] = df['orig_ground_station'].astype('Int64')
    df['orig_mh1_1'] = df['orig_mh1_1'].astype('Int64')
    df['orig_mh2_1'] = df['orig_mh2_1'].astype('Int64')
    df['orig_mhz_1'] = df['orig_mhz_1'].astype('Int64')
    df['orig_mh1_2'] = df['orig_mh1_2'].astype('Int64')
    df['orig_mh2_2'] = df['orig_mh2_2'].astype('Int64')
    df['orig_mhz_2'] = df['orig_mhz_2'].astype('Int64')
    df['orig_mh1_3'] = df['orig_mh1_3'].astype('Int64')
    df['orig_mh2_3'] = df['orig_mh2_3'].astype('Int64')
    df['orig_mhz_3'] = df['orig_mhz_3'].astype('Int64')
    df['orig_mh1_4'] = df['orig_mh1_4'].astype('Int64')
    df['orig_mh2_4'] = df['orig_mh2_4'].astype('Int64')
    df['orig_mhz_4'] = df['orig_mhz_4'].astype('Int64')
    # df['bit_synchronizer'] = df['bit_synchronizer'].astype('string')
    df['sync'] = df['sync'].astype('string')
    df['frame'] = df['frame'].astype('Int64')
    df['bit_synchronizer'] = df['bit_synchronizer'].astype('Int64')

    df['command'] = df['command'].astype('Int64')
    df['map'] = df['map'].astype('Int64')
    df['data_rate'] = df['data_rate'].astype('Int64')

    
    if new_format == 0:
        df['shz_4'] = df['shz_4'].astype('Int64')
        df['shz_6'] = df['shz_6'].astype('Int64')
        df['shz_8'] = df['shz_8'].astype('Int64')


        df['shz_10'] = df['shz_10'].astype('Int64')
        df['shz_12'] = df['shz_12'].astype('Int64')
        df['shz_14'] = df['shz_14'].astype('Int64')
        df['shz_16'] = df['shz_16'].astype('Int64')
        df['shz_18'] = df['shz_18'].astype('Int64')


        df['shz_20'] = df['shz_20'].astype('Int64')
        df['shz_22'] = df['shz_22'].astype('Int64')
        df['shz_24'] = df['shz_24'].astype('Int64')
        df['shz_26'] = df['shz_26'].astype('Int64')
        df['shz_28'] = df['shz_28'].astype('Int64')


        df['shz_30'] = df['shz_30'].astype('Int64')
        df['shz_32'] = df['shz_32'].astype('Int64')
        df['shz_34'] = df['shz_34'].astype('Int64')
        df['shz_36'] = df['shz_36'].astype('Int64')
        df['shz_38'] = df['shz_38'].astype('Int64')


        df['shz_40'] = df['shz_40'].astype('Int64')
        df['shz_42'] = df['shz_42'].astype('Int64')
        df['shz_44'] = df['shz_44'].astype('Int64')
        df['shz_46'] = df['shz_46'].astype('Int64')
        df['shz_48'] = df['shz_48'].astype('Int64')


        df['shz_50'] = df['shz_50'].astype('Int64')
        df['shz_52'] = df['shz_52'].astype('Int64')
        df['shz_54'] = df['shz_54'].astype('Int64')

        df['shz_58'] = df['shz_58'].astype('Int64')


        df['shz_60'] = df['shz_60'].astype('Int64')
        df['shz_62'] = df['shz_62'].astype('Int64')
        df['shz_64'] = df['shz_64'].astype('Int64')

    print(filename)
    df = clean_year(filename,df)


    if count_ground_station0 > 0:
        logging.info('SEVERE: {} ground station =0'.format(count_ground_station0))
        print('SEVERE: {} ground station =0'.format(count_ground_station0))

    print(out_filename_gz)
    df.to_csv(out_filename_gz,index=False,date_format='%Y-%m-%dT%H:%M:%S.%fZ',quoting=csv.QUOTE_NONNUMERIC)
    
    return

def import_csv_main_tapes(base_dir, out_base_dir, filenames=None):
    # filenames (default = None) - if specified, then read only these filenames 

    # if no filenames have been passed, then look for them 
    if filenames is None:
        filenames = []
        for filename in glob.glob(os.path.join(base_dir,'wtn*.gz')):
            filenames.append(os.path.basename(filename))

    for filename in filenames: 
        print('Searching for :', os.path.join(base_dir,filename))
        for in_file in glob.glob(os.path.join(base_dir,filename)):    
            # remove the .gz from the filename
            out_filename_gz, _ = os.path.splitext(filename)
            out_filename_gz = os.path.join(out_base_dir, '{}.csv.gz'.format(out_filename_gz))
            print(in_file, out_filename_gz)
            binary_import_main_tapes(in_file, out_filename_gz)
