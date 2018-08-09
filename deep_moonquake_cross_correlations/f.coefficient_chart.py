#!/usr/bin/env python

"""

Plot example

:copyright:
    The PDART Development Team & Katja Heger & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

import obspy
from obspy import read
from obspy.core import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.cross_correlation import xcorr_pick_correction, correlate, xcorr_max

plt.figure(figsize=(16,8))

# base_dir = '/Users/nunn/lunar_data/PDART'
base_dir = '/Users/nunn/Google Drive/for_Katja/PDART2'

st10 = obspy.read(os.path.join(base_dir,'1973/XA/S12/MHZ/XA.S12..MHZ.1973.006.gz'))
st20 = obspy.read(os.path.join(base_dir,'1973/XA/S12/MHZ/XA.S12..MHZ.1973.060.gz'))

st1=st10.copy()
st2=st20.copy()

sta=st10.copy()
stb=st20.copy()

t1= UTCDateTime("1973-01-06T05:39:12.209122Z")
t2= UTCDateTime("1973-03-01T07:18:03.829105Z")-1.020545 #-correction from x-correlation

start_time=250
window_length=3600-start_time

st1.merge()
st1=st1.trim(starttime=t1-start_time, endtime=t1+window_length)
st2.merge()
st2=st2.trim(starttime=t2-start_time, endtime=t2+window_length)#+0.1)

tr1 = st1.select(component="Z")[0]
tr2 = st2.select(component="Z")[0]

tr1.detrend()
tr2.detrend()

tr1.filter("bandpass", freqmin = 0.3, freqmax=0.9,corners=3)
tr2.filter("bandpass", freqmin = 0.3, freqmax=0.9,corners=3)

times1=tr1.times()

ax1=plt.subplot(211)
plt.plot(times1, tr1.data, linestyle="-", marker=None, c='r',linewidth=1, label='Day=006')
plt.plot(times1, tr2.data, linestyle="-", marker=None, c='g',linewidth=1, label='Day=060')
plt.ylim(-5,5)
plt.grid(axis='y')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=2, ncol=2, mode="expand", borderaxespad=0.)

k=10

for j in range (0,window_length+start_time,k):
    sta.merge()
    sta.trim(starttime=t1-2000, endtime=t1+10000)
    stb.merge()
    stb.trim(starttime=t2-2000, endtime=t2+10000)

    tra = sta.select(component="Z")[0]
    trb = stb.select(component="Z")[0]

    tra.detrend()
    trb.detrend()

    tra.filter("bandpass", freqmin = 0.1, freqmax=1.0,corners=3)
    trb.filter("bandpass", freqmin = 0.1, freqmax=1.0,corners=3)

#    print(sta)
#    print(stb)

    ta= UTCDateTime("1973-01-06T05:39:12.209122Z")+j-start_time
    tb= UTCDateTime("1973-03-01T07:18:03.829105Z")-1.020545+j-start_time

    trc=tra.copy()
    trd=trb.copy()

    trc.trim(starttime=ta, endtime=ta+k)
    trd.trim(starttime=tb, endtime=tb+k)

    l=j+k/2

    cc2 = correlate(trc, trd, 2)
    shift, value2 = xcorr_max(cc2)

    plt.subplot(212,sharex=ax1)
    plt.bar(l,value2,width=9)
    plt.grid()
    plt.axvline(x=start_time,color='k', linestyle="--", marker=None, linewidth=1.0)
    plt.grid(b=None, which='major', axis='both')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=2, ncol=2, mode="expand", borderaxespad=0.)



plt.xlabel('Time in seconds')
plt.subplots_adjust(hspace=0.01)
#plt.subplots(sharex=True)
plt.show()
