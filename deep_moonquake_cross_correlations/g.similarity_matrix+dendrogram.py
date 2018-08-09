#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculate and plot similarity matrix for A1 cluster (in 1973)

:copyright:
    The PDART Development Team & Katja Heger & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""


import obspy
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
from obspy import read_events
from obspy import read, Catalog, Trace, Stream
import os.path
import select
import io, urllib
from scipy.cluster import hierarchy
from scipy.spatial import distance
from obspy.signal.rotate import rotate2zne




def plot_similarity_matrix():

    # base_dir = '/Users/nunn/lunar_data/PDART'
    base_dir = '/Users/nunn/Google Drive/for_Katja/PDART2'
    catalog_file = '../LunarCatalog_Nakamura_1981_and_updates_v1/LunarCatalog_Nakamura_1981_and_updates_v1_A01.xml'
    catalog = read_events(catalog_file)
    stream1 = Stream(traces=[])
    stream2 = Stream(traces=[])
    i=0
    j=0
    k=0
    list=[]

    nx = 24
    ny = 24
    x = range(0,nx+1)
    y = range(0,ny+1)


    X, Y = np.meshgrid(x, y, indexing='xy')
    correlation = np.zeros(shape=(len(x),len(y)))


    for ev in catalog:
        picks = ev.picks
        for pick in picks:
            pick_time = pick.time
            t01 = pick.time
            startday=UTCDateTime(year=pick_time.year,julday=pick_time.julday)
            year = pick.time.year
            julday = startday.julday
            station = pick.waveform_id.station_code
            channel = pick.waveform_id.channel_code
            if station == 'S12' and year==1973 and channel == 'MHZ':
                try:
                    filename='XA.%s..MHZ.%s.%03d.gz' % ( station, str(year),julday )
                    filepath = os.path.join(base_dir,str(year),'XA',station,'MHZ',filename)
                except Exception as e: print(e)

                st1 = read(filepath)
                st1 = st1.trim(starttime=pick_time-1200, endtime=pick_time+3600)
            #    st1 = st1.merge()
                tra = st1.select(component="Z")[0]
                tr1 = tra.copy()
                tr1.detrend()
                tr1.filter("bandpass", freqmin = 0.3, freqmax=0.9,corners=3)
                i=i+1
                print('    i = {0}'.format(i))
                print(tr1)
                stream1.append(tr1)
                stream2.append(tr1)

    print(('    STREAM = {0}'.format(stream1)))

    for trj in stream1:
        for trk in stream2:
            cc = correlate(trj, trk, 31802)
            shift, value = xcorr_max(cc)
            list.append(value)
            j=j+1

    # generate random correlation matrix (only one half required because
    # it is symmetric)
    for x1 in x:
        for y1 in y[x1:]:
                k=k+1
                if x1 == y1:
                    correlation[x1,y1] = 1
                else:
                    correlation[x1,y1] = abs(list[k-1])
                # copy to the bottom half of the matrix
                correlation[y1,x1] = correlation[x1,y1]

    fig, (ax0, ax1) = plt.subplots(ncols=2)

    im = ax0.pcolormesh(X, Y, correlation)
    fig.colorbar(im, ax=ax0)
    ax0.set_title('Similarity')
    ax0.set_aspect('equal')

    dissimilarity = distance.squareform(1 - correlation)
    threshold = 0.3
    linkage = hierarchy.linkage(dissimilarity, method="single")
    clusters = hierarchy.fcluster(linkage, threshold, criterion="distance")

    ax1 = plt.subplot(122)
    ax1.set_title('Dendrogram')
    # hierarchy.dendrogram(linkage, color_threshold=0.3)
    hierarchy.dendrogram(linkage)
    plt.xlabel("Event number")
    plt.ylabel("Dissimilarity")
    plt.show()


if __name__ == "__main__":
    plot_similarity_matrix()
