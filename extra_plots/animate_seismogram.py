#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Code to make a video of an Apollo seismogram.



:copyright:
    Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from datetime import timedelta
from obspy.core import read
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory import read_inventory
import numpy as np
import pandas as pd
from pdart.util import linear_interpolation, remove_negative_ones
from pdart.view import read_solar_presence
import os
from datetime import datetime
import pdart.auth as auth


import matplotlib  
# matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import mlab
from matplotlib.colors import Normalize
import matplotlib.animation as animation

from obspy.imaging.cm import obspy_sequential
from obspy.imaging.spectrogram import _nearest_pow_2

DELTA = 0.1509433962

def get_stream():
    """Snippet to read in and remove the instrument response for Apollo.

    Example code, with example filtering parameters, to read in and 
    remove the instrument response from an Apollo seismogram. 
    
    """

    from obspy import UTCDateTime
    from obspy.clients.fdsn.client import Client
    
    # Apollo 15
    # onset
    starttime= UTCDateTime('1971-07-29T20:59:37.000000Z')
    # start 6 mins earlier (1 min will be trimmed later from beginning and end)
    starttime=starttime - 60*6
    # 2 hour record plus two mins for trimming
    endtime = starttime + 3600*2 + 2*60
    network='XA'
    station='S12'
    channel='MHZ'
    location='*'

    user = auth.user
    auth_password=auth.auth_password
    client = Client("IRIS",user=user,password=auth_password)

    if channel not in ['ATT']:
        # get the response file (wildcards allowed)
        inv = client.get_stations(starttime=starttime, endtime=endtime,
            network=network, sta=station, loc=location, channel=channel,
            level="response")

    stream = client.get_waveforms(network=network, station=station, channel=channel, location=location, starttime=starttime, endtime=endtime)
    for tr in stream:
        # interpolate across the gaps of one sample 
        linear_interpolation(tr,interpolation_limit=1)
    stream.merge()

    # plot removing the response, if required
    plot_response = False
    for tr in stream:
        # optionally interpolate across any gap 
        # for removing the instrument response from a seimogram, 
        # it is useful to get a mask, then interpolate across the gaps, 
        # then mask the trace again. 
        if tr.stats.channel in ['MH1', 'MH2', 'MHZ']:

            # add linear interpolation but keep the original mask
            original_mask = linear_interpolation(tr,interpolation_limit=None)
            # remove the instrument response
            pre_filt = [0.1,0.3,0.9,1.1]
            tr.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP",
                       water_level=None, plot=plot_response)
            if plot_response:
                plt.show()
            # apply the mask back to the trace 
            tr.data = np.ma.masked_array(tr, mask=original_mask)

        elif tr.stats.channel in ['SHZ']:

            # add linear interpolation but keep the original mask
            original_mask = linear_interpolation(tr,interpolation_limit=None)
            # remove the instrument response
            pre_filt = [1,2,11,13] 
            tr.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP",
                       water_level=None, plot=plot_response)
            if plot_response:
                plt.show()
            
            # apply the mask back to the trace 
            tr.data = np.ma.masked_array(tr, mask=original_mask)

    stream.trim(starttime=starttime+60,endtime=endtime-60)

    return stream

def animate_seismogram():

    # returns a 2 hour record 
    stream = get_stream()
    print(stream)
    # stream.plot(size=(600,800))
    # exit()


    # short version for testing 
    # stream[0].data = stream[0].data[0:310]

    xmin = 0
    xmax = stream[0].times()[-1]
    max = 1.05*(abs(stream[0].data)).max()
    ymin = -max
    ymax = max

    # creating a blank window
    # for the animation 
    fig = plt.figure(figsize=(8, 3)) 
    axis = plt.axes(xlim =(xmin, xmax),
                    ylim =(ymin, ymax)) 
    # turn the axis off (remove this line for testing)
    plt.axis('off')

    line, = axis.plot([], [], lw = 1) 
       
    # initiate the line
    def init(): 
        line.set_data([], []) 
        return line, 
       
    # animation function 
    def animate(i): 
        # add 100 new points every sample 
        line.set_data(stream[0].times()[0:i*100], stream[0].data[0:i*100])
        return line,

    print('End time ', stream[0].times()[-1])

    # normal speed for video is 24 fps
    # which is 41 ms interval, so we need something around this speed 
    # don't try something where the fps is much higher than 24 fps, 
    # because the file will be large and slow, and we won't even 
    # see the improvement

    # The real sampling interval is 0.1509433962 s
    # I want to display 2 hours in 18 s (400 times faster)
    # 0.1509433962/400*100 = 0.03774 s = 37.74 ms 

    speedup_factor = 400 
    interval = 100 * DELTA/speedup_factor # interval in miliseconds

    print('speedup_factor ', speedup_factor)
    print('interval ', interval)

    interval=37.74 # ms

    # calling the animation function     
    anim = animation.FuncAnimation(fig, animate, init_func = init, 
                                   interval=interval, blit = True, repeat=False, frames=(int(len(stream[0].data)/100)+1))

    # view the animation
    plt.draw()
    # plt.show()

    # saves the animation
    anim.save('../extra_plots_output/seismoApollo15XX.mp4', writer = 'ffmpeg')


if __name__ == "__main__":

    animate_seismogram()

