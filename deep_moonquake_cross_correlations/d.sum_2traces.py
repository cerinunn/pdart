#!/usr/bin/env python

"""

Plot example of cross correlated traces.

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
from obspy.signal.cross_correlation import xcorr_pick_correction
from os.path import join

dir1 = '/Users/nunn/lunar_data/PDART'
# dir1 = '/Users/nunn/Google Drive/for_Katja/PDART2'



st1 = obspy.read(join(dir1,'1973/XA/S12/MHZ/XA.S12..MHZ.1973.006.gz'))
st2 = obspy.read(join(dir1,'1973/XA/S12/MHZ/XA.S12..MHZ.1973.060.gz'))

#print(st1)
#print(st2)

t1= UTCDateTime("1973-01-06T05:39:12.209122Z")
t2= UTCDateTime("1973-03-01T07:18:03.829105Z")-1.020545 #-correction from x-correlation

st1.merge()
st1.trim(starttime=t1-100, endtime=t1+599.8)
st2.merge()
st2.trim(starttime=t2-100, endtime=t2+600)

tr1 = st1.select(component="Z")[0]
tr2 = st2.select(component="Z")[0]

tr1.detrend()
tr2.detrend()

tr1.filter("bandpass", freqmin = 0.3, freqmax=0.5,corners=3)
tr2.filter("bandpass", freqmin = 0.3, freqmax=0.5,corners=3)

times1=tr1.times()
#times2=tr2.times()


print(tr1)
print(tr2)


stack = np.sum([tr1.data, tr2.data] ,axis=0)

plt.figure(figsize=(10,6))

plt.subplot(411)
plt.plot(times1, tr1.data, linestyle="-", marker=None, c='r',linewidth=1.5, label='Day=006')
plt.ylim(-10,10)
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1, ncol=2, mode="expand", borderaxespad=0.)

plt.subplot(412)
plt.plot(times1, tr2.data, linestyle="-", marker=None, c='g',linewidth=1.5, label='Day=060')
plt.ylim(-10,10)
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1, ncol=2, mode="expand", borderaxespad=0.)

plt.subplot(413)
plt.plot(times1, tr1.data, linestyle="-", marker=None, c='r',linewidth=1, label='Day=006')
plt.plot(times1, tr2.data, linestyle="-", marker=None, c='g',linewidth=1, label='Day=060')
plt.ylim(-10,10)
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1, ncol=2, mode="expand", borderaxespad=0.)

plt.subplot(414)
plt.plot(times1, stack.data, linestyle="-", marker=None, c='k',linewidth=1.5, label='Stack')
plt.ylim(-10,10)
plt.grid()
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1, ncol=2, mode="expand", borderaxespad=0.)
plt.xlabel('Time in seconds')


plt.show()
