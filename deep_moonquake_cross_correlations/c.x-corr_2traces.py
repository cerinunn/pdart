#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Calculate the pick correction. 

:copyright:
    The PDART Development Team & Katja Heger & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

import obspy
from obspy.signal.cross_correlation import xcorr_pick_correction
from obspy import UTCDateTime
from obspy import read_events
from obspy import read
from os.path import join

dir1 = '/Users/nunn/lunar_data/PDART'
# dir1 = '/Users/nunn/Google Drive/for_Katja/PDART2'

st1 = obspy.read(join(dir1,'1973/XA/S12/MHZ/XA.S12..MHZ.1973.006.gz'))
st2 = obspy.read(join(dir1,'1973/XA/S12/MHZ/XA.S12..MHZ.1973.060.gz'))


print(st1)
print(st2)

t1= UTCDateTime("1973-01-06T05:39:12.209122Z")
t2= UTCDateTime("1973-03-01T07:18:03.829105Z")

# I'm using large windows in case I want to change the parameters
st1.merge()
st1.trim(starttime=t1-1200, endtime=t1+3600)
st2.merge()
st2.trim(starttime=t2-1200, endtime=t2+3600)

tr1 = st1.select(component="Z")[0]
tr2 = st2.select(component="Z")[0]

tr1.filter("bandpass", freqmin = 0.3, freqmax=0.5, corners=3)
tr2.filter("bandpass", freqmin = 0.3, freqmax=0.5, corners=3)

print(tr1)
print(tr2)

dt, coeff = xcorr_pick_correction(t1, tr1, t2, tr2, -10, 600, 10, plot=True)

print("  Time correction for pick 2: %.6f" % dt)
print("  Correlation coefficient: %.2f" % coeff)
