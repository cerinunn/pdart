#!/usr/bin/env python

'''
Plot the Statistics. The figure didn't make the paper, 
but the numbers were included in a table. 

grep 'Statistics' /Users/cnunn/lunar_data/PDART_CSV/S12/*.log > /Users/cnunn/lunar_data/PDART_CSV/S12/Statistics.txt
grep 'Statistics' /Users/cnunn/lunar_data/PDART_CSV/S14/*.log > /Users/cnunn/lunar_data/PDART_CSV/S14/Statistics.txt
grep 'Statistics' /Users/cnunn/lunar_data/PDART_CSV/S15/*.log > /Users/cnunn/lunar_data/PDART_CSV/S15/Statistics.txt
grep 'Statistics' /Users/cnunn/lunar_data/PDART_CSV/S16/*.log > /Users/cnunn/lunar_data/PDART_CSV/S16/Statistics.txt
grep 'Statistics' /Users/cnunn/lunar_data/PDART_CSV/work_tapes/*.log > /Users/cnunn/lunar_data/PDART_CSV/work_tapes/Statistics.txt
'''


from __future__ import print_function
import numpy as np
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
from matplotlib import gridspec
from pdart.view import plot_availabilty, save_availabilty
from obspy.core import read 
from obspy.core.utcdatetime import UTCDateTime



import matplotlib  
# matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt


def plot_statistics():



    # main tapes, S12 
    dict_12 = {'station': 'S12',
        'mission_start': UTCDateTime('1969-11-19 14:23:22.94300'),
        # 'mission_end': UTCDateTime('1977-09-30'),
        # TEMP MISSION END for testing 
        'mission_end': UTCDateTime('1976-03-01 00:01:44.696000'),
        'no_of_stations' : 1,

        'file' :'/Users/cnunn/lunar_data/PDART_CSV/S12/statistics.txt',
        'final_valid' : 0,
        'sync' : 0,
        'timing' : 0,
        'initial': 0,
        'duration' : 0
}

    # main tapes, S14
    dict_14 = {'station': 'S14',
        'mission_start': UTCDateTime('1971-02-05 17:44:33.349000'),
        # 'mission_end': UTCDateTime('1977-09-30'),
        # TEMP MISSION END for testing 
        'mission_end': UTCDateTime('1976-03-01 00:01:44.696000'),
        'no_of_stations' : 1,

        'file' :'/Users/cnunn/lunar_data/PDART_CSV/S14/statistics.txt',
        'final_valid' : 0,
        'sync' : 0,
        'timing' : 0,
        'initial': 0,
        'duration' : 0
}

    # main tapes, S15 
    dict_15 = {'station': 'S15',
        'mission_start': UTCDateTime('1971-07-31 18:52:00.269000'),
        # 'mission_end': UTCDateTime('1977-09-30'),
        # TEMP MISSION END for testing 
        'mission_end': UTCDateTime('1976-03-01 00:01:44.696000'),
        'no_of_stations' : 1,

        'file' :'/Users/cnunn/lunar_data/PDART_CSV/S15/statistics.txt',
        'final_valid' : 0,
        'sync' : 0,
        'timing' : 0,
        'initial': 0,
        'duration' : 0
}

    # main tapes, S16
    dict_16 = {'station': 'S16',
        'mission_start': UTCDateTime('1972-04-21 20:20:00.407000'),
        # 'mission_end': UTCDateTime('1977-09-30'),
        # TEMP MISSION END for testing 
        'mission_end': UTCDateTime('1976-03-01 00:01:44.696000'),
        'no_of_stations' : 1,

        'file' :'/Users/cnunn/lunar_data/PDART_CSV/S16/statistics.txt',
        'final_valid' : 0,
        'sync' : 0,
        'timing' : 0,
        'initial': 0,
        'duration' : 0
}

    # work tapes, S12, S14, S15, S16
    dict_work = {'station': 'Work',
        'mission_start': UTCDateTime('1976-03-01 00:01:44.696000'),
        # 'mission_end': UTCDateTime('1977-09-30'),
        # TEMP MISSION END for testing 
        'mission_end': UTCDateTime('1977-09-30T05:51:41.866000Z'),
        'no_of_stations' : 4,

        'file' :'/Users/cnunn/lunar_data/PDART_CSV_WORK_TAPES_FULL/Statistics.txt',
        'final_valid' : 0,
        'sync' : 0,
        'timing' : 0,
        'initial': 0,
        'duration' : 0
}



    # list_dict = [dict_12]
    # list_dict = [dict_14]
    # list_dict = [dict_15]
    # list_dict = [dict_16]
    # list_dict = [dict_work]

    
    list_dict = [
        dict_12, 
        dict_14,
        dict_15,
        dict_16,
        dict_work
    ]

    total_final_valid=0
    total_sync=0
    total_timing=0
    total_initial=0
    total_number_records=0
    total_overall_samples=0
    total_not_recovered=0

    for l in list_dict:
        duration = l['mission_end'] - l['mission_start']
        duration_days=duration/(3600*24)

        # l[duration] = l['no_of_stations']*duration_days

        final_valid=0
        sync=0
        timing=0
        initial=0
        number_records=0
        with open(l['file']) as f:
            for line in f:
                if 'Final Valid' in line:
                    recover = line.split()
                    final_valid1 = int(recover[17].split('=')[1])
                    sync1 = int(recover[5].split('=')[1])
                    timing1 = int(recover[12].split('=')[1])
                    initial1 = int(recover[2].split('=')[1])

                    final_valid+=final_valid1
                    sync+=sync1
                    timing+=timing1
                    initial+=initial1
                    number_records+=1

        # l['final_valid'] = final_valid
        # l['sync'] = sync
        # l['timing'] = timing
        # l['not_recovered'] = not_recovered
        overall_samples = duration * (1/0.6038) * l['no_of_stations']
        not_recovered = overall_samples - (final_valid + sync + timing)
        
        valid_percent = 100*final_valid/overall_samples
        sync_percent = 100*sync/overall_samples
        timing_percent = 100*timing/overall_samples
        not_recovered_percent=100*not_recovered/overall_samples

        # print(final_valid+sync+not_recovered+timing,overall_samples)
        # print(not_recovered,overall_samples-final_valid)
        

        print('{} Duration={:.0f} Overall samples={:.0f} Final Valid={} Sync Errors={} Timing Errors={} Not Recovered={:.0f} Percent Valid={:.1f}%, Percent Sync Error={:.1f}%,, Percent Not Recovered={:.1f}%'.format(
            l['station'],duration_days,overall_samples,final_valid,sync,timing,not_recovered,valid_percent,sync_percent,not_recovered_percent))

        total_final_valid+=final_valid
        total_sync+=sync
        total_timing+=timing
        total_initial+=initial
        total_number_records+=number_records
        total_overall_samples+=overall_samples
        total_not_recovered+=not_recovered


    total_valid_percent = 100*total_final_valid/total_overall_samples
    total_sync_percent = 100*total_sync/total_overall_samples
    total_timing_percent = 100*total_timing/total_overall_samples
    total_not_recovered_percent=100*total_not_recovered/total_overall_samples

    # print('100?',total_valid_percent+total_sync_percent+total_timing_percent+total_not_recovered_percent )
        

    total_valid_percent = 100*total_final_valid/total_overall_samples
    print('Total Overall samples={:.0f} Final Valid={} Sync Errors={} Timing Errors={} Not Recovered={:.0f} Percent Valid={:.2f}%, Sync Error Percent={:.2f}%, Timing Error Percent={:.2f}%,, Not Recovered Percent={:.2f}%'.format(
            total_overall_samples,total_final_valid,total_sync,total_timing,total_not_recovered,total_valid_percent,total_sync_percent,total_timing_percent, total_not_recovered_percent))


    # From the Tol color palette
    # https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40
    tol_dark_blue ='#341B88' 
    tol_green = '#357932'
    tol_blue_green = '#4CAA9A'
    tol_pale_blue = '#82C0E1'

    
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = 'Valid Records', 'Not Recovered (Damaged Transmission)', 'Not Recovered (Timing Issues)', 'Not Recorded'
    sizes = [total_final_valid,total_sync,total_timing, total_not_recovered ]
    print(sizes)
    explode = (0.1, 0.1, 0.1, 0.1)  

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, explode=None, labels=labels, 
            # autopct='%1.1f%%',
            shadow=False, startangle=90, wedgeprops=dict(width=0.4), 
            colors=[tol_dark_blue,tol_pale_blue,tol_blue_green,tol_green],
            labeldistance=None)
    # Note that labeldistance = None supresses the labels
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.legend(title='Data Recovery')

    labels = [f'{label} {np.round(size/sum(sizes)*100,1)}%' for label,size in zip(labels,sizes)]
    plt.legend(title='Data Recovery', labels=labels, fontsize=14, title_fontsize=14)
               # bbox_to_anchor=(1,1))
    # plt.suptitle('Data Recovery',fontsize=16)
    plt.tight_layout()
    plt.savefig('../extra_plots_output/DataRecovery_XX.png')
    plt.show()

# possibly think about something like this for the legend.
# plt.legend(labels=[f'{x} {np.round(y/sum(totalAmount_sample)*100,1)}%' for x,y in crimeTypes.items()], 
#            bbox_to_anchor=(1,1))


if __name__ == "__main__":
    plot_statistics()





