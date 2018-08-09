#!/usr/bin/env python

"""

Plot example seismograms, including bad examples.

:copyright:
    The PDART Development Team & Katja Heger & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

from obspy import read
import matplotlib.pyplot as plt
from os.path import join

dir1 = '/Users/nunn/lunar_data/PDART'

plt.figure(figsize=(8,8))

plt.subplot(5,2,1)
st1 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.03.MH1.1973.156.gz'))
print(st1)
st1.detrend()
plt.plot(st1[0])
plt.ylim(-10,10)


plt.subplot(5,2,2)
st2 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.01.MH1.1973.213.gz'))
print(st2)
st2.detrend()
plt.plot(st2[1])
plt.ylim(-10,10)

plt.subplot(5,2,3)
st3 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.06.MH1.1973.321.gz'))
print(st3)
st3.detrend()
plt.plot(st3[0])
plt.ylim(-10,10)

plt.subplot(5,2,4)
st4 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.05.MH1.1973.006.gz'))
print(st4)
st4.detrend()
plt.plot(st4[0])
plt.ylim(-10,10)

plt.subplot(5,2,5)
st5 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.05.MH1.1973.022.gz'))
print(st5)
st5.detrend()
plt.plot(st5[0])
plt.ylim(-10,10)

plt.subplot(5,2,6)
st6 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.14.MH1.1973.032.gz'))
print(st6)
st6.detrend()
plt.plot(st6[2])
plt.ylim(-10,10)

plt.subplot(5,2,7)
st7 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.02.MH1.1973.034.gz'))
print(st7)
st7.detrend()
plt.plot(st7[0])
plt.ylim(-10,10)

plt.subplot(5,2,8)
st8 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.02.MH1.1973.036.gz'))
print(st8)
st8.detrend()
plt.plot(st8[5])
plt.ylim(-10,10)

plt.subplot(5,2,9)
st9 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.14.MH1.1973.050.gz'))
print(st9)
st9.detrend()
plt.plot(st9[0])
plt.ylim(-10,10)

plt.subplot(5,2,10)
st10 = read(join(dir1,'1973/XA/S12/MH1/XA.S12.09.MH1.1973.060.gz'))
print(st10)
st10.detrend()
plt.plot(st10[0])
plt.ylim(-10,10)

plt.show()
