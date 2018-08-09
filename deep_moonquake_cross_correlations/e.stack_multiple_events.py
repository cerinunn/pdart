#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

:copyright:
    The PDART Development Team & Katja Heger & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""


import obspy
from obspy.signal.cross_correlation import xcorr_pick_correction
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
from obspy import read_events
from obspy import read, Catalog
import os.path
import select



base_dir = '/Users/nunn/lunar_data/PDART'
# base_dir = '/Users/nunn/Google Drive/for_Katja/PDART2'
catalog_file = '../LunarCatalog_Nakamura_1981_and_updates_v1/LunarCatalog_Nakamura_1981_and_updates_v1_A01.xml'



catalog = read_events(catalog_file)

stack=0
i=0

for ev in catalog:

    picks = ev.picks

    for pick in picks:
        pick_time = pick.time
        t0 = pick.time
        startday=UTCDateTime(year=pick_time.year,julday=pick_time.julday)
        year = pick.time.year
        julday = startday.julday
        station = pick.waveform_id.station_code
        channel = pick.waveform_id.channel_code
        # if station == 'S12' and year==1973 and channel == 'MHZ':
        if station == 'S12' and year==1973:
            try:
                filename='XA.%s..MHZ.%s.%03d.gz' % ( station, str(year),julday )
                filepath = os.path.join(base_dir,str(year),'XA',station,'MHZ',filename)
            except Exception as e: print(e)

            st = read(filepath)

            st=st.trim(starttime=pick_time-1200, endtime=pick_time+3600)
            tr = st.select(component="Z")[0]
            tr.detrend()
            tr.filter("bandpass", freqmin = 0.3, freqmax=0.9,corners=3)

            st1 = obspy.read(os.path.join(base_dir, '1973/XA/S12/MHZ/XA.S12..MHZ.1973.006.gz'))
            t1= UTCDateTime("1973-01-06T05:37")
            st1=st1.trim(starttime=t1-1200, endtime=t1+3600)
            tr1 = st1.select(component="Z")[0]
            tr1.detrend()
            tr1.filter("bandpass", freqmin = 0.3, freqmax=0.9,corners=3)

            print('---> Start')
            print('this is trimmed')
            print(tr1)
            print(tr)

            # try:
            #     dt, coeff = xcorr_pick_correction(t1, tr1, t0, tr, -600, 600, 600, plot=False)
            #     print('made it past xcorr')
            # except Exception as e : print('cross correlation error')
            dt, coeff = xcorr_pick_correction(t1, tr1, t0, tr, -600, 600, 600, plot=False)
            if coeff>0.6:
                print('  Filename:{0}'.format(filename))
                print("  Time correction for pick 2: %.6f" % dt)
                print("  Correlation coefficient: %.2f" % coeff)
                tnew=t0+dt
                trnew=tr.trim(starttime=tnew-10, endtime=tnew+1200)
                tr2=tr1.trim(starttime=t1-10, endtime=t1+1197.4)
                times1=tr2.times()

                print('tr-NEW')
                print(trnew)

                try:
                    stack+= trnew[:8000]
                except Exception as e : print(e)
                trnew=trnew[:8000]
                i=i+1
                print('i = {0}'.format(i))
                print()

                try:
                    fig, axs = plt.subplots(10,1, figsize=(15, 6), facecolor='w', edgecolor='k')
                    fig.subplots_adjust(hspace = .001, wspace=.001)

                    axs = axs.ravel()

                    fig.axes.get_xaxis().set_visible(False)
                    fig.axes.get_yaxis().set_visible(False)

                    plt.plot(times1, trnew.data, linestyle="-", marker=None, c='k',
                    linewidth=1.5, label='julday={0}'.format(pick_time.julday))
                    plt.ylim(-10,10)
                    plt.grid()
                    plt.legend()
                    axs[i].set_title(str(250+i))
                except Exception as e : print(e)

# plt.subplot(10,1,10)
# plt.plot(times1, stack.data, linestyle="-", marker=None, c='g',linewidth=1)
# plt.grid()
# plt.xlabel('time in seconds')

plt.show()
