#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Import the data from the binary tape files - old format data.
Be careful to make sure that every part of the record is read - otherwise
the file positions will be out for the next record.


# These packages are currently being run on env36. 




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
from datetime import datetime, timedelta
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
        'pse.a12.10.91',
        'pse.a12.10.92',
        'pse.a12.10.93',
        'pse.a12.10.94',
        'pse.a12.10.95',
        'pse.a12.10.97',
        'pse.a12.10.98',
        'pse.a12.10.102',
        'pse.a12.10.103',
        'pse.a12.10.104',
        'pse.a12.10.106',
        'pse.a12.10.107',
        'pse.a12.10.108',
        'pse.a12.10.109',
        'pse.a12.10.111']


    # tape = os.path.basename(filename)
    # s1=tape.split('.')[1]
    # t1=int(tape.split('.')[2])
    # t2=int(tape.split('.')[3])
    #
    # # s12 1 43
    #
    # if s1 == 'a11':
    #     year_out = 1969
    # elif s1 == 'a12' and t1 ==1 and t2 in list(range(1,44)):
    #     year_out = 1969
    # elif s1 == 'a12' and t1 ==1:
    #     year_out = 1970
    # elif s1 == 'a12' and t1 ==2:
    #     year_out = 1970
    # elif s1 == 'a12' and t1 ==3 and t2 in list(range(1,62)):
    #     year_out = 1970
    # elif s1 == 'a12' and t1 ==3:
    #     year_out = 1971
    # elif s1 == 'a12' and t1 ==4:
    #     year_out = 1971
    # elif s1 == 'a12' and t1 ==5  and t2 in list(range(1,59)):
    #     year_out = 1971
    # elif s1 == 'a12' and t1 ==5:
    #     year_out = 1972
    # elif s1 == 'a12' and t1 ==6 and t2 in list(range(1,136)):
    #     year_out = 1972
    # elif s1 == 'a12' and t1 ==6:
    #     year_out = 1973
    # elif s1 == 'a12' and t1 ==7 and t2 in list(range(1,211)):
    #     year_out = 1973
    # elif s1 == 'a12' and t1 ==7:
    #     year_out = 1974
    # elif s1 == 'a12' and t1 ==8:
    #     year_out = 1974
    # elif s1 == 'a12' and t1 ==9:
    #     year_out = 1975
    # elif s1 == 'a12' and t1 ==10  and t2 in list(range(1,78)):
    #     year_out = 1975
    # elif s1 == 'a12' and t1 ==10:
    #     year_out = 1976
    # #######
    # elif s1 == 'a14' and t1 ==1:
    #     year_out = 1971
    # elif s1 == 'a14' and t1 ==2  and t2 in list(range(1,151)):
    #     year_out = 1971
    # elif s1 == 'a14' and t1 ==2:
    #     year_out = 1972
    # elif s1 == 'a14' and t1 ==3:
    #     year_out = 1972
    # elif s1 == 'a14' and t1 ==4 and t2 in list(range(1,166)):
    #     year_out = 1972
    # elif s1 == 'a14' and t1 ==4:
    #     year_out = 1973
    # elif s1 == 'a14' and t1 ==5:
    #     year_out = 1973
    # elif s1 == 'a14' and t1 ==6:
    #     year_out = 1973
    # elif s1 == 'a14' and t1 ==7:
    #     year_out = 1974
    # elif s1 == 'a14' and t1 ==8:
    #     year_out = 1974
    # elif s1 == 'a14' and t1 ==9 and t2 in list(range(1,7)):
    #     year_out = 1974
    # elif s1 == 'a14' and t1 ==9:
    #     year_out = 1975
    # elif s1 == 'a14' and t1 ==10 and t2 in list(range(1,20)):
    #     year_out = 1975
    # elif s1 == 'a14' and t1 ==10:
    #     year_out = 1976
    #
    # ########
    # elif s1 == 's15' and t1 ==1 and t2 in list(range(1,155)):
    #     year_out = 1971
    # elif s1 == 's15' and t1 ==1:
    #     year_out = 1972
    # elif s1 == 's15' and t1 ==2:
    #     year_out = 1972
    # elif s1 == 's15' and t1 ==3 and t2 in list(range(1,162)):
    #     year_out = 1972
    # elif s1 == 's15' and t1 ==3:
    #     year_out = 1973
    # elif s1 == 's15' and t1 ==4:
    #     year_out = 1973
    # elif s1 == 's15' and t1 ==5:
    #     year_out = 1973
    # elif s1 == 's15' and t1 ==6:
    #     year_out = 1974
    # elif s1 == 's15' and t1 ==7:
    #     year_out = 1974
    # elif s1 == 's15' and t1 ==8  and t2 in list(range(1,7)):
    #     year_out = 1974
    # elif s1 == 's15' and t1 ==8:
    #     year_out = 1975
    # elif s1 == 's15' and t1 ==9:
    #     year_out = 1975
    # elif s1 == 's15' and t1 ==10  and t2 in list(range(1,23)):
    #     year_out = 1975
    # elif s1 == 's15' and t1 ==10:
    #     year_out = 1976
    #
    # ##########
    #
    # elif s1 == 'a16' and t1 ==1:
    #     year_out = 1972
    # elif s1 == 'a16' and t1 ==2  and t2 in list(range(1,85)):
    #     year_out = 1972
    # elif s1 == 'a16' and t1 ==2:
    #     year_out = 1973
    # elif s1 == 'a16' and t1 ==3:
    #     year_out = 1973
    # elif s1 == 'a16' and t1 ==4 and t2 in list(range(1,94)):
    #     year_out = 1973
    # elif s1 == 'a16' and t1 ==4:
    #     year_out = 1974
    # elif s1 == 'a16' and t1 ==5:
    #     year_out = 1974
    # elif s1 == 'a16' and t1 ==6 and t2 in list(range(1,100)):
    #     year_out = 1974
    # elif s1 == 'a16' and t1 ==6:
    #     year_out = 1975
    # elif s1 == 'a16' and t1 ==7:
    #     year_out = 1975
    # elif s1 == 'a16' and t1 ==8 and t2 in list(range(1,116)):
    #     year_out = 1975
    # elif s1 == 'a16' and t1 ==8:
    #     year_out = 1976
    #
    # if year != year_out:
    #     print('Year does not match file ', os.path.basename(filename), year, year_out)
    #
    # print(year, year_out)

    if os.path.basename(filename) in file_list:
        year = 1976
        print('Year corrected for ', filename, year)

    return year

# XXXXXXXXXX
#   File "./run_binary_import_work_tapes.py", line 21, in <module>
#     import_csv_work_tapes(base_dir=base_dir, out_base_dir=out_base_dir)
#   File "/Users/cnunn/python_packages/pdart/binary_import_work_tapes.py", line 919, in import_csv_work_tapes
#     binary_import(filename, out_filename)
#   File "/Users/cnunn/python_packages/pdart/binary_import_work_tapes.py", line 428, in binary_import
#     tape_type = s.read(16).int
#   File "/opt/anaconda3/envs/env36/lib/python3.6/site-packages/bitstring.py", line 3883, in read
#     fmt, self.len - self._pos)
# bitstring.ReadError: Cannot read 16 bits, only 0 available.

def binary_import_work_tapes(filename, out_filename_gz,
    log_filename='logs/binary.log',
    logging_level=logging.INFO):
    """
    """

# Bytes Information (16-bit integers)
# 1-2 1 for PSE tapes; 2 for Event tapes
# 3-4 Apollo station number
# 5-6 original tape sequence number for PSE tapes; 2-digit station code plus 3-digit
# original event tape sequence number for Event tapes 7-8 record number
# 9-10 year
# 11-12 format (0 = old, 1 = new)
# 13-14 number of physical records from original tape
# 15-16 original tape read error flags (lsb set if read error occurred while reading the first
# physical record from the original tape, etc)

    logging.basicConfig(filename=log_filename, filemode='w', level=logging_level)

    t0 = datetime.now()

    # read the data as gzipped data
    in_file = gzip.GzipFile(filename, 'rb')
    gzipped_data = in_file.read()
    in_file.close()

    ground_stations_all = set(range(1,15))
    ground_stations = set()
    year = None

    # # open the output file
    # with gzip.open(out_filename_gz, 'wt') as f:
    #     writer = csv.writer(f)
    #     writer.writerow(HEADER)

    # open the BitStream from the file, into fixed-length chunks
    # Each file consists of two (occasionally one) 16-byte headers followed by 5760-byte data records.
    # 8*(5760+16+16)=46,336 bits
    # chunks = ConstBitStream(gzipped_data).cut(46336)
    # one frame = 8*96=768

    s = ConstBitStream(gzipped_data)

    # if you need the files unzipped
    # chunks = ConstBitStream(filename=filename).cut(46336)


    # read the header
    # 1-2 3 to identify Normal-Bit-Rate Work tape



    # if you need the files unzipped
    # chunks = ConstBitStream(filename=filename).cut(46336)
    # chunks = ConstBitStream(gzipped_data).cut(46336)

    # # for testing, view just a few chunks
    # import itertools
    # # chunks = itertools.islice(chunks, 3)
    # chunks = itertools.islice(chunks, 3)

    # search for hexadecimal '0x1'

    # for x, s in enumerate(chunks):
    #     print('Chunk ', x, s.pos)
        # if x == 1:
        #     print('more than one chunk, so exiting ')
        #     exit()

        # if x == 1:
        #     end_of_file = s.read(4)
        #     print('end_of_file' , end_of_file)
            # print(end_of_file.int)
            # 
            # hexy = s.read(32)
            # print('hexy2' , hexy)

            # end_of_file = s.read(4)
            # print('end_of_file' , end_of_file)
            # 
            # end_of_file = s.read(8)
            # print('end_of_file' , end_of_file)
            # exit()
            
        #     print('This is the result ', s.read(1).int)

        # # PSE Tapes 
        # Header format:
        # Bytes Information (16-bit integers)
        # 1-2 1 for PSE tapes; 2 for Event tapes
        # 3-4 Apollo station number
        # 5-6 original tape sequence number for PSE tapes; 2-digit station code plus 3-digit
        # original event tape sequence number for Event tapes 7-8 record number
        # 9-10 year
        # 11-12 format (0 = old, 1 = new)
        # 13-14 number of physical records from original tape
        # 15-16 original tape read error flags (lsb set if read error occurred while reading the first
        # physical record from the original tape, etc)

        ########

        # Work Tapes
        # Header format:
        # Bytes Information (16-bit integers except for last 6 bytes)
        # 1-2 3 to identify Normal-Bit-Rate Work tape
        # 3-4 active station code 1, 2, 3, 4, and 5 for Apollo 12, 15, 16, 14 and 17 stations,
        # respectively, placed in consecutive 3-bit positions starting at the second most
        # significant bit. The most significant bit is not used. 
        # 5-6 number of active stations
        # 7-8 original 9-track tape ID number (7001 through 8109)
        # 9-10 year
        # 11-14 and first 4 bits of byte 15 time of the year of the first data in msec, 36-bit integer Remaining bits of byte 15 and byte 16 not used

    # while True:

    physical_record = 0

    dict_list = []

    try:
        tape_type = s.read(16).int

        # 3-4 active station code 1, 2, 3, 4, and 5 for Apollo 12, 15, 16, 14 and 17 stations,
        # respectively, placed in consecutive 3-bit positions starting at the second most
        # significant bit. The most significant bit is not used. 
        s.read(1).int
        active_station_code = s.read(3).int
        # try:
        #     active_station = station_dict[active_station_code]
        # except KeyError:
        #     active_station = 'NaN'

        # XXX1
            
        s.read(12).int
        # 5-6 number of active stations
        active_startion_no = s.read(16).int
        if active_startion_no != 5:
            # the code is written for 5 active stations - will need 
            # modifying if there is a different number 
            print("No of active stations ", active_startion_no)
            # return
        # original 9-track tape ID number (7001 through 8109)
        s.read(16).int
        # 9-10 year
        year =  s.read(16).int
        # check this at the beginning of the file
        # if year is None:
        #     # some years are not specified correctly
        #     year = year_corrections(filename,year1)
        # 11-14 and first 4 bits of byte 15 time of the year of the first data in msec, 36-bit integer Remaining bits of byte 15 and byte 16 not used

        time =  s.read(36).int
        # the time is miliseconds after the beginning of the year
        # but with the first record at 86400 seconds
        orig_timestamp_file = UTCDateTime(year=year,julday=1)+time/1000-86400.

        s.read(12).int
        # 3 1 5 1977 20826354812 1977-08-29T01:05:54.812000Z 241
        # print(tape_type, active_station_code, active_startion_no, year, time, orig_timestamp_file, orig_timestamp_file.julday )       
        # 3 S12 5 1976 5270460000 1976-03-01T00:01:00.000000Z 61

        # from Jaxa webservice
        # id: 3
        # active station: 1,2,3,4,5, ( )
        # Number of active stations: 5
        # Original 9 track ID: 7001
        # Year: 1976
        # First msec: 5270460000 top prev next last
        # Active Station Offset: 0 1 2 3 4

        if 'wtn.6.30.gz' in filename: 
            # second header is missing from this file 
            pass
        else: 
            # read the second header
            s.read(128).int
        
        # s.pos is 256 bits
        # print(s.pos)
        

        orig_frame = 1
        orig_no = None
        dict = None

# 'orig_no','orig_timestamp','orig_frame','orig_station',
#   'clock_flag','orig_ground_station',
# 'orig_mh1_1','orig_mh2_1','orig_mhz_1',
# 'orig_mh1_2','orig_mh2_2','orig_mhz_2',
# 'orig_mh1_3','orig_mh2_3','orig_mhz_3',
# 'orig_mh1_4','orig_mh2_4','orig_mhz_4']


        t1 = datetime.now()

        if tape_type == 3:
# 46080
            while True:

                try: 
                    clock_flag = s.read(1)
                except ReadError:
                    print('Last row', dict)
                    print ('eof', s.pos)
                    break



                    # Example of last real record:
                    # 4369,60,1976-02-22T00:34:58.449000Z,0,60,1976-02-22T00:34:58.449000Z,S17,1,273,546,68,68,273,546,546,68,273,273,546,68

                # read in long words 1 and 2
                clock_flag = clock_flag.uint

                time = s.read(35)
                bin_time = time.bin
                time = time.uint



                if time != 0:
                    # the time is miliseconds after the beginning of the year
                    # but with the first record at 86400 seconds
                    orig_timestamp = UTCDateTime(year=year,julday=1)+time/1000-86400.
                    
                else:
                    # if the time is missing, set it to Not a Time 
                    # orig_timestamp = 'NaT'
                    orig_timestamp = pd.NaT


                orig_ground_station = s.read(4)
                orig_ground_station = orig_ground_station.uint
                if orig_ground_station not in ground_stations_all:
                    # orig_ground_station = 'NaN'
                    orig_ground_station = pd.NA

                package_id= s.read(3)
                # bin_package_id = package_id.bin
                package_id = package_id.uint
                
                bit_synchronizer = s.read(5).bin

                record_no = s.read(16)
                # bin_record_no = record_no.bin
                record_no = record_no.uint
                # record_no is only set on frame 1, but it is useful for all frames
                if orig_no is None:
                    orig_no = record_no

                
                # print('temporary')
                # if record_no > 1:
                #     orig_timestamp = 'NaT'
                #     print('timestamping')

            


                # if frames == 1: 
                # flag bit: 0
                # msec of year: 5270504696 (doy=61 time=00:01:44.696) ===> 1976-03-01T00:01:44.696
                # ALSEP tracking ID: 9
                # ALSEP package id: 1 (Apollo12)
                # Bit synchronizer status: search:0 verify:0 confirm:0 lock:1 input level:1
                # Original record: 1 sync pattern barker code: 1810 (11100010010)
                # sync pattern barker code compliment: 237 (00011101101)
                # frame count: 1
                # mode bit: 1

                # print(clock_flag, orig_timestamp, ground_station, package_id,record_no)
                # 0 1976-03-01T00:01:44.696000Z 9 1 3 1
                # 0 1976-03-01T00:01:44.638000Z 9 2 3 0
                # 0 1976-03-01T00:01:44.378000Z 9 3 3 0
                # 0 1976-03-01T00:01:44.373000Z 9 4 3 0
                # 0 1976-03-01T00:01:44.671000Z 9 5 3 0
    
                if time != 0:
                    ground_stations.add(orig_ground_station)
    
                # s.pos is 192 bits
    
                # tape - long word 3 
                word1 = s.read(10) # word1
                word1a = s.read(1)  # not used
                word2 = s.read(10) # word2
                word2a = s.read(1)  # not used
                word3 = s.read(10) # word3

                sync = word1 + word2 + word3[0:2]
                sync = sync.bin
                frame = word3[2:9].uint
                mode_bit = word3[9:10].uint

                # tape - long word 4
                word4 = s.read(10) # word4
                word4a = s.read(1)  # not used
                word5 = s.read(10) # word5
                word5a = s.read(1)  # not used
                word6 = s.read(10) # word6

                # tape - long word 5
                word7 = s.read(10) # word7
                word7a = s.read(1)  # not used
                word8 = s.read(10) # word8
                word8a = s.read(1)  # not used
                word9 = s.read(10) # word9

                # tape - long word 6
                word10 = s.read(10) # word10
                word10a = s.read(1)  # not used
                word11 = s.read(10) # word11
                word11a = s.read(1)  # not used
                word12 = s.read(10) # word12

                # tape - long word 7
                word13 = s.read(10) # word13
                word13a = s.read(1)  # not used
                word14 = s.read(10) # word14
                word14a = s.read(1)  # not used
                word15 = s.read(10) # word15

                # tape - long word 8
                word16 = s.read(10) # word16
                word16a = s.read(1)  # not used
                word17 = s.read(10) # word17
                word17a = s.read(1)  # not used
                word18 = s.read(10) # word18

                # tape - long word 9
                word19 = s.read(10) # word19
                word19a = s.read(1)  # not used
                word20 = s.read(10) # word20
                word20a = s.read(1)  # not used
                word21 = s.read(10) # word21

                # tape - long word 10
                word22 = s.read(10) # word22
                word22a = s.read(1)  # not used
                word23 = s.read(10) # word23
                word23a = s.read(1)  # not used
                word24 = s.read(10) # word24

                # tape - long word 11
                word25 = s.read(10) # word25
                word25a = s.read(1)  # not used
                word26 = s.read(10) # word26
                word26a = s.read(1)  # not used
                word27 = s.read(10) # word27

                # tape - long word 12
                word28 = s.read(10) # word28
                word28a = s.read(1)  # not used
                word29 = s.read(10) # word29
                word29a = s.read(1)  # not used
                word30 = s.read(10) # word30

                # tape - long word 13
                word31 = s.read(10) # word31
                word31a = s.read(1)  # not used
                word32 = s.read(10) # word32
                word32a = s.read(1)  # not used
                word33 = s.read(10) # word33

                # tape - long word 14
                word34 = s.read(10) # word34
                word34a = s.read(1)  # not used
                word35 = s.read(10) # word35
                word35a = s.read(1)  # not used
                word36 = s.read(10) # word36

                # tape - long word 15
                word37 = s.read(10) # word37
                word37a = s.read(1)  # not used
                word38 = s.read(10) # word38
                word38a = s.read(1)  # not used
                word39 = s.read(10) # word39

                # tape - long word 16
                word40 = s.read(10) # word40
                word40a = s.read(1)  # not used
                word41 = s.read(10) # word41
                word41a = s.read(1)  # not used
                word42 = s.read(10) # word42

                # tape - long word 17
                word43 = s.read(10) # word43
                word43a = s.read(1)  # not used
                word44 = s.read(10) # word44
                word44a = s.read(1)  # not used
                word45 = s.read(10) # word45

                word46 = s.read(10) # word46
                word46a = s.read(1)  # not used
                word47 = s.read(10) # word47
                word47a = s.read(1)  # not used
                word48 = s.read(10) # word48

                # tape - long word 19
                word49 = s.read(10) # word49
                word49a = s.read(1)  # not used
                word50 = s.read(10) # word50
                word50a = s.read(1)  # not used
                word51 = s.read(10) # word51

                # tape - long word 20
                word52 = s.read(10) # word52
                word52a = s.read(1)  # not used
                word53 = s.read(10) # word53
                word53a = s.read(1)  # not used
                word54 = s.read(10) # word54

                # tape - long word 21
                word55 = s.read(10) # word55
                word55a = s.read(1)  # not used
                word56 = s.read(10) # word56
                word56a = s.read(1)  # not used
                word57 = s.read(10) # word57

                # tape - long word 22
                word58 = s.read(10) # word58
                word58a = s.read(1)  # not used
                word59 = s.read(10) # word59
                word59a = s.read(1)  # not used
                word60 = s.read(10) # word60

                # tape - long word 23
                word61 = s.read(10) # word61
                word61a = s.read(1)  # not used
                word62 = s.read(10) # word62
                word62a = s.read(1)  # not used
                word63 = s.read(10) # word63

                # tape - long word 24
                word64 = s.read(10) # word64
                word64a = s.read(22)  # not used

                orig_mh1_1 = word9.uint
                orig_mh2_1 = word11.uint
                orig_mhz_1 = word13.uint

                orig_mh1_2 = word25.uint
                orig_mh2_2 = word27.uint
                orig_mhz_2 = word29.uint

                orig_mh1_3 = word41.uint
                orig_mh2_3 = word43.uint
                orig_mhz_3 = word45.uint

                orig_mh1_4 = word57.uint
                orig_mh2_4 = word59.uint
                orig_mhz_4 = word61.uint

                # SHZ - all even except 2, 46, 56
                # (also excepting ALSEP word 24 for station 15)

                shz_4 = word4.uint
                shz_6 = word6.uint
                shz_8 = word8.uint

                shz_10 = word10.uint
                shz_12 = word12.uint
                shz_14 = word14.uint
                shz_16 = word16.uint
                shz_18 = word18.uint

                shz_20 = word20.uint
                shz_22 = word22.uint
                shz_24 = word24.uint
                shz_26 = word26.uint
                shz_28 = word28.uint

                shz_30 = word30.uint
                shz_32 = word32.uint
                shz_34 = word34.uint
                shz_36 = word36.uint
                shz_38 = word38.uint

                shz_40 = word40.uint
                shz_42 = word42.uint
                shz_44 = word44.uint
                # shz_46 is below
                shz_48 = word48.uint

                shz_50 = word50.uint
                shz_52 = word52.uint
                shz_54 = word54.uint

                shz_58 = word58.uint

                shz_60 = word60.uint
                shz_62 = word62.uint
                shz_64 = word64.uint

                try:
                    orig_station = station_dict[package_id]
                except KeyError:
                    # orig_station = 'NaN'
                    orig_station = pd.NA

                    # word46 used for SHZ at S14, and for command and 
                    # map at the other stations
                try:
                    if orig_station == 'S14':
                        shz_46 = word46.uint
                        command1 = word5[2:9].uint
                        map1 = word5[9:10].uint
                    else:
                        command1 = word46[2:9].uint
                        map1 = word46[9:10].uint
                        shz_46 = None
                except:
                    command1 = word46[2:9].uint
                    map1 = word46[9:10].uint   
                    shz_46 = None     

                try:
                    if orig_station == 'S15':
                        # *also excepting ALSEP word 24 for station 15
                        shz_24 = None
                except:
                    pass

                # if orig_station == 'S16':
                #     orig_station = ''
                    # if the station code is missing, we can set it manually
                    # from the frame number 
                    # print('Setting station manually')
                    # if orig_frame in (1,6,11,16,21,26,31,36,41,46,51,56):
                    #     corr_station = 'S12'
                    # elif orig_frame in (2,7,12,17,22,27,32,37,42,47,52,57):
                    #     corr_station = 'S15'
                    # elif orig_frame in (3,8,13,18,23,28,33,38,43,48,53,58):
                    #     corr_station = 'S16'
                    # elif orig_frame in (4,9,14,19,24,29,34,39,44,49,54,59):
                    #     corr_station = 'S14'
                    # elif orig_frame in (5,10,15,20,25,30,35,40,45,50,55,60):
                    #     corr_station = 'S17'

                dict = {
                  'orig_no' : orig_no,
                  'orig_timestamp' : str(orig_timestamp),
                  'orig_frame' : orig_frame,
                  'orig_station' : orig_station,
                  'clock_flag' : clock_flag,
                  'orig_ground_station' : orig_ground_station,
                  'orig_mh1_1' : orig_mh1_1,
                  'orig_mh2_1' : orig_mh2_1,
                  'orig_mhz_1' : orig_mhz_1,
                  'orig_mh1_2': orig_mh1_2,
                  'orig_mh2_2' : orig_mh2_2,
                  'orig_mhz_2' : orig_mhz_2, 
                  'orig_mh1_3' : orig_mh1_3,
                  'orig_mh2_3' : orig_mh2_3,
                  'orig_mhz_3' : orig_mhz_3, 
                  'orig_mh1_4' : orig_mh1_4,
                  'orig_mh2_4' : orig_mh2_4,
                  'orig_mhz_4' : orig_mhz_4, 
                  'bit_synchronizer' : bit_synchronizer,
                  'frame' : frame,
                  'sync' : sync,

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

                  'command' : command1,
                  'map' : map1 
                    
                }
                dict_list.append(dict)

                # print(frames, s.pos)
                if orig_frame == 60:
                    # print('End of record')            
                    orig_frame = 1
                    orig_no = None
                else: 
                    orig_frame+=1

                # # to run just a short test 
                # if physical_record == 179:
                #     print('Temporary break')
                #     break

                # to run just a long test 
                # if physical_record == 3000:
                #     print('Temporary break')
                #     break


                    # if s.pos > (16*8 + physical_record*5760*8):
                    # print(s.pos)



                physical_record += 1



    except ReadError as ex:
        print(ex)
        print('Position of the last line: ', s.pos)

    print('len at end ', len(dict_list))
    df = pd.DataFrame.from_dict(dict_list)

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
    df['bit_synchronizer'] = df['bit_synchronizer'].astype('string')
    df['sync'] = df['sync'].astype('string')
    df['frame'] = df['frame'].astype('Int64')

    df['command'] = df['command'].astype('Int64')
    df['map'] = df['map'].astype('Int64')

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
    

    # print(df.dtypes)
    # print(out_filename_gz)
    df.to_csv(out_filename_gz,index=False,date_format='%Y-%m-%dT%H:%M:%S.%fZ',quoting=csv.QUOTE_NONNUMERIC)

    t2 = datetime.now()
    print('t2 normal : {}'.format((t2-t1).total_seconds()))

    print('t2 ALL : {}'.format((t2-t0).total_seconds()))

    
    return

# def split_and_delete_work_tapes(ground_stations,out_filename):
# 
#     # make a gzipped file for each of the ground stations and each of the stations
#     for station in ['S12','S15', 'S16', 'S14', 'S17']:
#         for ground_staton in ground_stations:
#             out_filename_gs = out_filename.replace('.csv', '{}_{:02d}.csv.gz'.format(station, ground_staton))
#             # temporary
#             out_filename1 = out_filename.replace('.csv', '{}_{:02d}.csv'.format(station, ground_staton))
# 
# 
#             with open(out_filename, 'r') as r1:
#                 with gzip.open(out_filename_gs, 'wt') as w1:
#                 # with open(out_filename1, 'wt') as w1:
#                     writer1 = csv.DictWriter(w1, fieldnames=HEADER)
#                     writer1.writeheader()
#                     reader = csv.DictReader(r1)
#                     for row in reader:
#                         # use the ground station and corrected station to put 
#                         # record into the right file 
#                         if row['orig_ground_station'] == str(ground_staton) and row['orig_station'] == station:
#                             writer1.writerow(row)
# 
# # XXXXXXX
# 
#     # delete the csv file
#     os.remove(out_filename)

def import_csv_work_tapes(base_dir, out_base_dir, filenames=None):
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
            binary_import_work_tapes(in_file, out_filename_gz)
                

    # if links_filename is not None:
    #     links = []
    # 
    #     print(os.path.join(base_dir,links_filename))
    # 
    #     with open(os.path.join(base_dir,links_filename), 'rt') as l1:
    #         for line in l1:
    #             links.append(line.strip('\n'))
    # 
    #     print(links)
    # 
    #     # print('temporary - just read one')
    #     # for link in links[0:1]:
    #     for link in links:
    #         filename=os.path.join(base_dir,link)
    #         out_filename_gz = os.path.basename(filename)
    #         # remove the .gz from the filename
    #         out_filename_gz, _ = os.path.splitext(out_filename_gz)
    #         out_filename_gz=os.path.join(out_base_dir,'{}.csv.gz'.format(out_filename_gz))
    #         print(filename, out_filename_gz)
    #         binary_import_work_tapes(filename, out_filename_gz)
    # 
    # else:
    #     # just read one
    #     # for filename in glob.glob(os.path.join(base_dir,'wtn*.gz'))[0:1]:




# /Users/cnunn/lunar_data/PDART_TAPES/work_tapes/wtn.1.1
# /Users/cnunn/lunar_data/PDART_CSV_WORK_TAPES/wtn.1.1.csv

        # test_end(filename, out_filename)
        # 
        # print(ground_stations   )
        # 
        # # temporary 
        # print('temporary')
        # split_and_delete_work_tapes(ground_stations,out_filename)

# not checked for this code 
# def test_end(filename, out_filename):
#     """
#     """
# 
#     # read the data as gzipped data
#     gzip_filename = '{}.gz'.format(filename)
#     in_file = gzip.GzipFile(gzip_filename, 'rb')
#     gzipped_data = in_file.read()
#     in_file.close()
# 
#     ground_stations = set()
#     year = None
# 
#     # open the output file
#     with open(out_filename, 'w') as f:
#         writer = csv.writer(f)
#         writer.writerow(HEADER)
# 
#         # open the BitStream from the file, into fixed-length chunks
#         # Each file consists of two (occasionally one) 16-byte headers followed by 5760-byte data records.
#         # 8*(5760+16+16)=46,336 bits
#         number_to_advance = 20000*1000
#         chunks = ConstBitStream(gzipped_data).cut(463368000)
# 
# 
#         # for testing, view just a few chunks
#         import itertools
#         chunks = itertools.islice(chunks, 3)
# 
#         for x, s in enumerate(chunks):
# 
#             if x > 0:
#                 exit()
# 
#             current_pos = s.pos
#             for advance in range(0,number_to_advance):
#                 s.pos = current_pos + advance
# 
#                 try: 
#                 # read the header
#                 # 1-2 3 to identify Normal-Bit-Rate Work tape
#                     tape_type = s.read(16).int
#                     # 3-4 active station code 1, 2, 3, 4, and 5 for Apollo 12, 15, 16, 14 and 17 stations,
#                     # respectively, placed in consecutive 3-bit positions starting at the second most
#                     # significant bit. The most significant bit is not used. 
#                     s.read(1).int
#                     active_station_code = s.read(3).int
#                     station_code = station_dict[active_station_code]
#                     s.read(12).int
#                     # active_station_code = s.read(16).int
#                     # 5-6 number of active stations
#                     active_startion_no = s.read(16).int
#                     # original 9-track tape ID number (7001 through 8109)
#                     s.read(16).int
#                     # 9-10 year
#                     year =  s.read(16).int
# 
# 
#                     # check this at the beginning of the file
#                     # if year is None:
#                     #     # some years are not specified correctly
#                     #     year = year_corrections(filename,year1)
#                     # 11-14 and first 4 bits of byte 15 time of the year of the first data in msec, 36-bit integer Remaining bits of byte 15 and byte 16 not used
# 
#                     time =  s.read(36).int
# 
#                     # TODO - THIS MAY NOT BE TRUE FOR THE WORK TAPE
#                     #             # the time is miliseconds after the beginning of the year
#                     #             # but with the first record at 86400 seconds
#                     #             orig_timestamp = UTCDateTime(year=year,julday=1)+time/1000-86400.
#                     orig_timestamp = UTCDateTime(year=year,julday=1)+time/1000-86400.
# 
#                     s.read(12).int
#                     # if year == 1976:
#                     #     print(advance, tape_type, station_code, active_startion_no, year, time, orig_timestamp, orig_timestamp.julday )
#                 except:
#                     pass


if __name__ == "__main__":


    pass

    # print(UTCDateTime(year=1971,julday=1))
    # print(UTCDateTime('1975-10-17T00:00:00.344000Z').julday)
    # exit()

    # filename='/Users/nunn/lunar_data/PDART_TAPES/S15/pse.a15.1.1'
    # out_filename='/Users/nunn/lunar_data/PDART_CSV/S15/a15.1.1.csv'
    # out_filename='/Users/nunn/lunar_data/PDART_CSV/S15/a15.1.1XX.csv'


    # base_dir='/Users/nunn/lunar_data/PDART_TAPES/S15/'
    # out_base_dir='/Users/nunn/lunar_data/PDART_CSV/S15/'
    # # binary_import(filename, out_filename)
    # import_csv(base_dir=base_dir, out_base_dir=out_base_dir)


    # base_dir='/Users/nunn/lunar_data/PDART_TAPES/S12/'
    # out_base_dir='/Users/nunn/lunar_data/PDART_CSV/S12/'
    # # binary_import(filename, out_filename)
    # import_csv(base_dir=base_dir, out_base_dir=out_base_dir)

    # base_dir='/Users/nunn/lunar_data/PDART_TAPES/S15/temp/'
    # out_base_dir='/Users/nunn/lunar_data/PDART_TAPES/S15/temp/'
    # links_filename = '/Users/nunn/lunar_data/PDART_TAPES/S15/links2.txt'
    # links_filename = None
    # # binary_import(filename, out_filename)
    # import_csv(base_dir=base_dir, out_base_dir=out_base_dir, link_file=link_file)


# 
# TODO this worktape is defintely incorrect:
# /Users/cnunn/lunar_data/PDART_TAPES/work_tapes/wtn.1.12
# 3 S12 5 1976 5709736631 1976-03-06T02:02:16.631000Z 66
# 256
# 0 1976-03-06T02:02:16.631000Z 6 1 3 2201
# 0 1976-03-06T02:02:16.775000Z 6 2 3 0
