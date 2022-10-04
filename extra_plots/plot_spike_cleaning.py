#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import os

import matplotlib
from matplotlib import pyplot as plt

from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import stream_from_directory_new
from obspy.imaging.util import (_set_xaxis_obspy_dates, _id_key, _timestring)
from matplotlib.dates import date2num
from pdart.util import relative_timing_trace, relative_timing_stream, remove_negative_ones
from matplotlib.dates import AutoDateFormatter, AutoDateLocator, HOURLY
from obspy.core import Stream
import pandas as pd
import random
# from pdart.csv_join_work_tapes import stream_import, initial_cleanup, merge_channel_stream, despike3, loose_frame_diff, read_file, process_list
from pdart.csv_join_work_tapes import despike3, to_Int64, stream_import
from pdart.csv_check_work_tapes import calculate_gaps, add_or_minus_frame
import pdart.config as config

# default date format 
DATEFORMAT='%Y-%m-%dT%H:%M:%S.%fZ'

ORIG_FRAMES = list(range(1,61))

# consecutive stations (the order is important)
STATIONS = ['S12', 'S15', 'S16', 'S14', 'S17']

def basic_timeseries_station(no_records=6):
    df = basic_timeseries(no_records=no_records)

    # make a timeseries with only 'S12' data in it 
    df = local_initial_cleanup(df)
    # detailed_report(df,'')
    grouped_st = df.groupby('orig_station')
    df_st = grouped_st.get_group('S12')
    df_st = df_st.copy(deep=True)

    # # make a new index
    df_st.reset_index(inplace=True,drop=True)

    return df_st

def basic_timeseries(no_records=3):
    random.seed(a=2)
    # orig_no,orig_timestamp,orig_frame,orig_station,clock_flag,orig_ground_station,orig_mh1_1,orig_mh2_1,orig_mhz_1,orig_mh1_2,orig_mh2_2,orig_mhz_2,orig_mh1_3,orig_mh2_3,orig_mhz_3,orig_mh1_4,orig_mh2_4,orig_mhz_4
    # 1,1976-03-01T00:01:44.696000Z,1,S12,0,9,531,491,481,531,491,481,531,491,481,531,491,481

    orig_timestamp = [UTCDateTime('1976-03-01T00:01:44.696000Z'), 
      UTCDateTime('1976-03-01T00:01:44.696000Z'),
      UTCDateTime('1976-03-01T00:01:44.378000Z'),
      UTCDateTime('1977-03-01T00:01:44.373000Z'),
      UTCDateTime('1976-03-01T00:01:44.671000Z')
    ]

    orig_mh1_1 = [531, 534, 534, 560, 530]
    orig_mh2_1 = [490, 491, 492, 493, 494]
    orig_mhz_1 = [490, 491, 492, 493, 494]

    orig_mh1_2 = [531, 534, 534, 560, 530]
    orig_mh2_2 = [490, 491, 492, 493, 494]
    orig_mhz_2 = [490, 491, 492, 493, 494]

    orig_mh1_3 = [531, 534, 534, 560, 530]
    orig_mh2_3 = [490, 491, 492, 493, 494]
    orig_mhz_3 = [490, 491, 492, 493, 494]

    orig_mh1_4 = [531, 534, 534, 560, 530]
    orig_mh2_4 = [490, 491, 492, 493, 494]
    orig_mhz_4 = [490, 491, 492, 493, 494]
    

# 1,1976-03-01T00:01:44.696000Z,1,S12,0,9,531,491,481,531,491,481,531,491,481,531,491,481
# 1,1976-03-01T00:01:44.638000Z,2,S15,0,9,504,505,491,504,505,491,504,505,491,504,505,491
# 1,1976-03-01T00:01:44.378000Z,3,,0,9,511,512,469,511,512,469,510,511,469,510,511,470
# 1,1976-03-01T00:01:44.373000Z,4,S14,0,9,511,490,509,511,490,509,511,490,508,511,490,508
# 1,1976-03-01T00:01:44.671000Z,5,S17,0,9,0,0,0,514,504,0,0,0,0,0,0,0

    start_no = 150

    for i in range(60*no_records-5):
        # add a random number onto the last random.choice (which has an index of 
        # -5 because there are 5 stations)
        orig_timestamp.append(orig_timestamp[-5] + random.choice([0.603, 0.604, 0.604, 0.604, 0.604]))
    orig_timestamp = [x.strftime(DATEFORMAT) for x in orig_timestamp]

    # logging.info(STATIONS)
    # logging.info(ORIG_FRAMES)

    # make the data frame from orig_timestamp
    df = pd.DataFrame(orig_timestamp, columns=['orig_timestamp'])
    df['orig_timestamp'] = pd.to_datetime(df['orig_timestamp'], format=DATEFORMAT)

    # tile the frames until the end of the record
    df['orig_frame'] = np.tile(ORIG_FRAMES, len(df)//len(ORIG_FRAMES) + 1)[:len(df)]
    # tile the stations until the end of the record
    df['orig_station'] = np.tile(STATIONS, len(df)//len(STATIONS) + 1)[:len(df)]    

    # repeat the origin number for a whole physical record 
    df['orig_no'] = list(np.repeat(np.arange(start_no,start_no+no_records), 60))
    # set the software clock flag to zero 
    df['clock_flag'] = 0 
    df['orig_ground_station'] = 9
    df['bit_synchronizer'] = '00011'
    df['sync'] = '1110001001000011101101'
    frame = [3,10,13,20,23]

    for i in range(60*no_records-5):
        # add a random number onto the last random.choice (which has an index of 
        # -5 because there are 5 stations)
        orig_mh1_1.append(orig_mh1_1[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mh2_1.append(orig_mh2_1[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mhz_1.append(orig_mhz_1[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))

        orig_mh1_2.append(orig_mh1_2[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mh2_2.append(orig_mh2_2[-5] + random.choice([0, 0, -1, 1, 2, -2, 0, 0]))
        orig_mhz_2.append(orig_mhz_2[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))

        orig_mh1_3.append(orig_mh1_3[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mh2_3.append(orig_mh2_3[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mhz_3.append(orig_mhz_3[-5] + random.choice([0, 0, 1, -1, -2, 2, 0, 0]))

        orig_mh1_4.append(orig_mh1_4[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mh2_4.append(orig_mh2_4[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))
        orig_mhz_4.append(orig_mhz_4[-5] + random.choice([0, 0, -1, 1, -2, 2, 0, 0]))

        frame.append(add_or_minus_frame(frame[-5],1))

    df['orig_mh1_1'] = orig_mh1_1
    df['orig_mh2_1'] = orig_mh2_1
    df['orig_mhz_1'] = orig_mhz_1

    df['orig_mh1_2'] = orig_mh1_2
    df['orig_mh2_2'] = orig_mh2_2
    df['orig_mhz_2'] = orig_mhz_2

    df['orig_mh1_3'] = orig_mh1_3
    df['orig_mh2_3'] = orig_mh2_3
    df['orig_mhz_3'] = orig_mhz_3

    df['orig_mh1_4'] = orig_mh1_4
    df['orig_mh2_4'] = orig_mh2_4
    df['orig_mhz_4'] = orig_mhz_4

    df['frame'] = frame

    df['corr_ground_station'] = df['orig_ground_station']
    df['corr_timestamp'] = df['orig_timestamp']

    return df

def local_initial_cleanup(df):

    # df['orig_timestamp'] = df['orig_timestamp'].astype('datetime64[ns, UTC]')

    df['orig_station'] = df['orig_station'].astype('string')
    df['bit_synchronizer'] = df['bit_synchronizer'].astype('string')
    df['sync'] = df['sync'].astype('string')

    df['orig_no'] = to_Int64(df['orig_no'])
    df['clock_flag'] = to_Int64(df['clock_flag'])
    df['orig_frame'] = to_Int64(df['orig_frame'])
    df['orig_ground_station'] = to_Int64(df['orig_ground_station'])
    df['orig_mh1_1'] = to_Int64(df['orig_mh1_1'])
    df['orig_mh2_1'] = to_Int64(df['orig_mh2_1'])
    df['orig_mhz_1'] = to_Int64(df['orig_mhz_1'])
    df['orig_mh1_2'] = to_Int64(df['orig_mh1_2'])
    df['orig_mh2_2'] = to_Int64(df['orig_mh2_2'])
    df['orig_mhz_2'] = to_Int64(df['orig_mhz_2'])
    df['orig_mh1_3'] = to_Int64(df['orig_mh1_3'])
    df['orig_mh2_3'] = to_Int64(df['orig_mh2_3'])
    df['orig_mhz_3'] = to_Int64(df['orig_mhz_3'])
    df['orig_mh1_4'] = to_Int64(df['orig_mh1_4'])
    df['orig_mh2_4'] = to_Int64(df['orig_mh2_4'])
    df['orig_mhz_4'] = to_Int64(df['orig_mhz_4'])
    df['frame'] = to_Int64(df['frame'])

    # replace empty ground station and orig station with values
    df['orig_ground_station'].fillna(-1,inplace=True)
    df['orig_station'].fillna('S0',inplace=True)
    df['orig_frame'].fillna(-99,inplace=True)
    df['frame'].fillna(-99,inplace=True)

    df['corr_timestamp'] = df['corr_timestamp'].astype('datetime64[ns, UTC]')
    df['corr_ground_station'] = to_Int64(df['corr_ground_station'])

    df['time_index'] = 0

    # add new columns for timestamps
    # df['corr_timestamp'] = df['orig_timestamp']   
    df['corr_gap'] = np.NaN
    df['frame_change'] = 1
    df['corr_gap_count'] = 1
    df['corr_ground_station'] = df['orig_ground_station']
    df['cumsum'] = 0
    
    df['begin'] = False
    df['end'] = False

    df['delta4'] = 0.6038

    return df

def plot_spikes_removal():
    ##############################
    # add a single downward spike to an array with a few gaps
    df_gst = basic_timeseries_station()  
    starttime0 = df_gst.corr_timestamp.iloc[0]
    # starttime0 = pd.Timestamp(starttime0) 
    starttime0 = UTCDateTime(df_gst.corr_timestamp.iloc[0])
    print(starttime0)
    
    # add a spike 
    spikes = 0
    for row_i in range(25,26):
        df_gst.at[row_i,'orig_mh1_1'] = 100
        spikes += 1

    # print(df_gst.to_string())
    # exit()
    # delete some records 
    # to_drop = [28, 31, 32]
    # df_gst.drop(to_drop,inplace=True)
    
    # df_gst.reset_index(inplace=True,drop=True)
    
    df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst)  
    
    # Skipped timestamps should be added 
    # df_gst, df_dropped, rec_adjusted_timestamps, rec_deleted_timestamps = station_fix_timeskips(df_gst,GZIP_FILENAME)
    
    # df_gst, gaps_long, gaps_8888 = calculate_gaps(df_gst) 
    df_gst = local_initial_cleanup(df_gst)
    
    df_gst['time_index'] = np.arange(len(df_gst))
    
    stream = stream_import(df_gst,sample_time0=starttime0,index0=0,attempt_merge=False)
    stream_before = stream.select(channel='MH1')
    # replace with a sine wave 
    for tr in stream_before:
        tr.data = (511.2+1.5*np.sin(tr.times())).astype(int)
        # mask the array
        tr.data = np.ma.masked_array(data=tr.data,
                     mask=False)
        tr.data = tr.data[0:201]

    # add a small spike 
    # spikes = 0
    tr.data[24] = 498
        # spikes += 1

    # add a small spike before a gap 
    tr.data[49] = 520
    tr.data[50] = np.ma.masked

    # add a spike after a gap 
    tr.data[75] = np.ma.masked
    tr.data[76] = 498

    # add an up and down spike 
    tr.data[150] = 498 
    tr.data[151] = 520   

    # add a double spike 
    tr.data[175] = 514 
    tr.data[176] = 520   


    # print(tr.data)
    stream_despike = stream_before.copy()
    for tr in stream_despike:
        tr.stats.location = 'despike'

    config.clean_spikes = True
    # if required, clean the spikes 
    if config.clean_spikes:
        for tr in stream_despike:
            # print('trace ', len(tr.data), tr.id)
            despike3(tr)

    total_stream = stream_before + stream_despike

    fig = total_stream.plot(method='full',size=(1200,600),type='relative',handle=True,show=False)
    axes = plt.gcf().get_axes()
    for i, ax in enumerate(axes):
        ax.tick_params(axis='x', labelsize=14 )
        ax.tick_params(axis='y', labelsize=14 )
        if i == 1:
            ax.set_xlabel(xlabel='[s]',size=14)    
        ax.set_ylabel(ylabel='Digital Units',size=14) 
    plt.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.10)
    fig.savefig('spike_cleaningXX.png')
    plt.show()




if __name__ == "__main__":
    plot_spikes_removal()