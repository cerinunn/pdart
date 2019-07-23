#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
View the catalog

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)



"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
# from builtins import input
from obspy import read_events, read_inventory
from obspy.core import read, Stream
from obspy.core.event.origin import Pick
from obspy.core.event.base import WaveformStreamID
from obspy.core.util.misc import limit_numpy_fft_cache
from obspy.signal.filter import envelope
from obspy.core.utcdatetime import UTCDateTime
from future.utils import native_str
import os.path
import numpy as np
from matplotlib.dates import date2num
import math

from obspy.geodetics.base import gps2dist_azimuth
MOON_RADIUS = 1737 * 1000.
MOON_FLATTENING = 1/825

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})

STATION_LIST = ['S12','S14','S15','S16']
LATITUDES = [-3.01084,-3.6445,26.13407,-8.97577]
LONGITUDES = [-23.42456,-17.47753,3.62981,15.49649]
MH1_AZIMUTHS = [180,0,0,334.5]
MH2_AZIMUTHS = [270,90,90,64.5]

colors=['r','b','g','yellow','purple','pink']
pick_colors=['purple','pink']
 
filter_ids = ['smi:jpl/filter/Apollo/standard_0.25_0.75',
  'smi:jpl/filter/Apollo/standard_0.75_1.25',
  'smi:jpl/filter/Apollo/standard_1.25_1.75',
]

xdata = None

ampl_dict =	{
  "purple": 3e-9,
  "pink": 6e-9,
  "darkblue": 12e-9,
  "lightblue": 24e-9,
  "green": 48e-9,
  "yellow": 96e-9,
  "orange": 192e-9,
  "red": 384e-9,
  "n": 'Not detected',
  "na": 'Not available',
}


# ampl_dict_disp =	{
#   "purple": 2e-10,
#   "pink": 4e-10,
#   "darkblue": 8e-10,
#   "lightblue": 16e-10,
#   "green": 32e-10,
#   "yellow": 64e-10,
#   "orange": 256e-10,
#   "red": 512e-10,
#   "n": 'Not detected',
#   "na": 'Not available',
# }

AUTHORITY_ID = 'nakamura81' 

"""
View the catalog files
"""



def find_dir(top_level_dir,year,station,channel):
    return os.path.join(top_level_dir, str(year), 'XA', station, channel)

def find_processed_dir(top_level_dir,year,station):
    return os.path.join(top_level_dir, str(year), 'XA', station)

def find_seismogram(top_level_dir,starttime,endtime,
  stations=['S12','S14','S15','S16'],channels=['MH1', 'MH2', 'MHZ'], 
  dir_type='pdart_dir'):
    # TODO add more complex code to get the right days and the 
    # including if it is early in the morning


    for station in stations: 
        stream = Stream()
        channel='*'
        if dir_type=='processed_dir':
            dir = find_processed_dir(top_level_dir,starttime.year,station)
            filename = '*%s.%s.%s.%s.%s.%03d*.gz' % ('XA',station, '*', channel,
                str(starttime.year), starttime.julday)
        else:
            dir = find_dir(top_level_dir,starttime.year,station,channel)
            filename = '%s.%s.%s.%s.%s.%03d.gz' % ('XA',station, '*', channel,
                str(starttime.year), starttime.julday)
        filename = os.path.join(dir,filename)
        try:
            stream += read(filename)
        except Exception as e:
            print(str(e))

        if dir_type=='processed_dir':
            dir = find_processed_dir(top_level_dir,starttime.year,station)
            filename = '*%s.%s.%s.%s.%s.%03d*.gz' % ('XA',station, '*', channel,
                str(starttime.year), starttime.julday-1)
        filename = os.path.join(dir,filename)
        try:
            stream += read(filename)
        except Exception as e:
            print(str(e))

        if starttime.julday != endtime.julday:
            if dir_type=='processed_dir':
                dir = find_processed_dir(top_level_dir,endtime.year,station)
                filename = '*%s.%s.%s.%s.%s.%03d*.gz' % ('XA',station, '*', channel,
                    str(endtime.year), endtime.julday)

            else: 
                dir = find_dir(top_level_dir,endtime.year,station,channel)
                filename = '*%s.%s.%s.%s.%s.%03d*.gz' % ('XA',station, '*', channel,
                    str(endtime.year), endtime.julday)
            filename = os.path.join(dir,filename)
            try: 
                stream += read(filename)
            except Exception as e:
                print(str(e))

        stream = stream.trim(starttime=starttime,endtime=endtime)

        if stream is not None and len(stream) > 0: 
            for tr in stream: 
                tr.stats.location = ''
                if tr.stats.channel not in channels:
                    stream.remove(tr)

            stream.merge()

    return stream

# is this not returning the amended trace properly???
def remove_response(trace, inventory=None, output="VEL", water_level=60,
                    pre_filt=None, zero_mean=True, taper=True,
                    taper_fraction=0.05, plot=False, fig=None, 
                    return_spectra=False, **kwargs):

    """
    Deconvolve instrument response. Based on a method on trace.
    It adds units, which are obviously dependent on the problem.
    axis ax5b is modified to show the pre-filtered trace instead of the
    time domain pre filter.

    Uses the adequate :class:`obspy.core.inventory.response.Response`
    from the provided
    :class:`obspy.core.inventory.inventory.Inventory` data. Raises an
    exception if the response is not present.

    Note that there are two ways to prevent overamplification
    while convolving the inverted instrument spectrum: One possibility is
    to specify a water level which represents a clipping of the inverse
    spectrum and limits amplification to a certain maximum cut-off value
    (`water_level` in dB). The other possibility is to taper the waveform
    data in the frequency domain prior to multiplying with the inverse
    spectrum, i.e. perform a pre-filtering in the frequency domain
    (specifying the four corner frequencies of the frequency taper as a
    tuple in `pre_filt`).

    .. note::

        Any additional kwargs will be passed on to
        :meth:`obspy.core.inventory.response.Response.get_evalresp_response`,
        see documentation of that method for further customization (e.g.
        start/stop stage).

    .. note::

        Using :meth:`~Trace.remove_response` is equivalent to using
        :meth:`~Trace.simulate` with the identical response provided as
        a (dataless) SEED or RESP file and when using the same
        `water_level` and `pre_filt` (and options `sacsim=True` and
        `pitsasim=False` which influence very minor details in detrending
        and tapering).

    .. rubric:: Example

    >>> from obspy import read, read_inventory
    >>> st = read()
    >>> tr = st[0].copy()
    >>> inv = read_inventory()
    >>> tr.remove_response(inventory=inv)  # doctest: +ELLIPSIS
    <...Trace object at 0x...>
    >>> tr.plot()  # doctest: +SKIP

    .. plot::

        from obspy import read, read_inventory
        st = read()
        tr = st[0]
        inv = read_inventory()
        tr.remove_response(inventory=inv)
        tr.plot()

    Using the `plot` option it is possible to visualize the individual
    steps during response removal in the frequency domain to check the
    chosen `pre_filt` and `water_level` options to stabilize the
    deconvolution of the inverted instrument response spectrum:

    >>> from obspy import read, read_inventory
    >>> st = read("/path/to/IU_ULN_00_LH1_2015-07-18T02.mseed")
    >>> tr = st[0]
    >>> inv = read_inventory("/path/to/IU_ULN_00_LH1.xml")
    >>> pre_filt = [0.001, 0.005, 45, 50]
    >>> tr.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP",
    ...                    water_level=60, plot=True)  # doctest: +SKIP
    <...Trace object at 0x...>

    .. plot::

        from obspy import read, read_inventory
        st = read("/path/to/IU_ULN_00_LH1_2015-07-18T02.mseed", "MSEED")
        tr = st[0]
        inv = read_inventory("/path/to/IU_ULN_00_LH1.xml", "STATIONXML")
        pre_filt = [0.001, 0.005, 45, 50]
        output = "DISP"
        tr.remove_response(inventory=inv, pre_filt=pre_filt, output=output,
                           water_level=60, plot=True)

    :type inventory: :class:`~obspy.core.inventory.inventory.Inventory`
        or None.
    :param inventory: Station metadata to use in search for adequate
        response. If inventory parameter is not supplied, the response
        has to be attached to the trace with :meth:`Trace.attach_response`
        beforehand.
    :type output: str
    :param output: Output units. One of:

        ``"DISP"``
            displacement, output unit is meters
        ``"VEL"``
            velocity, output unit is meters/second
        ``"ACC"``
            acceleration, output unit is meters/second**2

    :type water_level: float
    :param water_level: Water level for deconvolution.
    :type pre_filt: list or tuple of four float
    :param pre_filt: Apply a bandpass filter in frequency domain to the
        data before deconvolution. The list or tuple defines
        the four corner frequencies `(f1, f2, f3, f4)` of a cosine taper
        which is one between `f2` and `f3` and tapers to zero for
        `f1 < f < f2` and `f3 < f < f4`.
    :type zero_mean: bool
    :param zero_mean: If `True`, the mean of the waveform data is
        subtracted in time domain prior to deconvolution.
    :type taper: bool
    :param taper: If `True`, a cosine taper is applied to the waveform data
        in time domain prior to deconvolution.
    :type taper_fraction: float
    :param taper_fraction: Taper fraction of cosine taper to use.
    :type plot: bool or str
    :param plot: If `True`, brings up a plot that illustrates how the
        data are processed in the frequency domain in three steps. First by
        `pre_filt` frequency domain tapering, then by inverting the
        instrument response spectrum with or without `water_level` and
        finally showing data with inverted instrument response multiplied
        on it in frequency domain. It also shows the comparison of
        raw/corrected data in time domain. If a `str` is provided then the
        plot is saved to file (filename must have a valid image suffix
        recognizable by matplotlib e.g. '.png').
    """
    limit_numpy_fft_cache()

    from obspy.core.inventory import PolynomialResponseStage
    from obspy.signal.invsim import (cosine_taper, cosine_sac_taper,
                                     invert_spectrum)
    if plot:
        import matplotlib.pyplot as plt

    if (isinstance(inventory, (str, native_str)) and
            inventory.upper() in ("DISP", "VEL", "ACC")):
        from obspy.core.util.deprecation_helpers import \
            ObsPyDeprecationWarning
        output = inventory
        inventory = None
        msg = ("The order of optional parameters in method "
               "remove_response has changed. 'output' is not accepted "
               "as first positional argument in the next release.")
        warnings.warn(msg, category=ObsPyDeprecationWarning, stacklevel=3)
    response = trace._get_response(inventory)
    # polynomial response using blockette 62 stage 0
    if not response.response_stages and response.instrument_polynomial:
        coefficients = response.instrument_polynomial.coefficients
        trace.data = np.poly1d(coefficients[::-1])(trace.data)
        return trace

    # polynomial response using blockette 62 stage 1 and no other stages
    if len(response.response_stages) == 1 and \
       isinstance(response.response_stages[0], PolynomialResponseStage):
        # check for gain
        if response.response_stages[0].stage_gain is None:
            msg = 'Stage gain not defined for %s - setting it to 1.0'
            warnings.warn(msg % trace.id)
            gain = 1
        else:
            gain = response.response_stages[0].stage_gain
        coefficients = response.response_stages[0].coefficients[:]
        for i in range(len(coefficients)):
            coefficients[i] /= math.pow(gain, i)
        trace.data = np.poly1d(coefficients[::-1])(trace.data)
        return trace

    # use evalresp
    data = trace.data.astype(np.float64)
    npts = len(data)
    # time domain pre-processing
    if zero_mean:
        data -= data.mean()
    if taper:
        data *= cosine_taper(npts, taper_fraction,
                             sactaper=True, halfcosine=False)

    if plot:
        color1 = "blue"
        color2 = "red"
        bbox = dict(boxstyle="round", fc="w", alpha=1, ec="w")
        bbox1 = dict(boxstyle="round", fc="blue", alpha=0.15)
        bbox2 = dict(boxstyle="round", fc="red", alpha=0.15)
        if fig is None:
            fig = plt.figure(figsize=(14, 10))
        fig.suptitle("{}:\nRemoving the Instrument Response".format(str(trace)),
            fontsize=16)
        ax1 = fig.add_subplot(321)
        ax1b = ax1.twinx()
        ax2 = fig.add_subplot(323, sharex=ax1)
        ax2b = ax2.twinx()
        ax3 = fig.add_subplot(325, sharex=ax1)
        ax3b = ax3.twinx()
        ax4 = fig.add_subplot(322)
        ax5 = fig.add_subplot(324, sharex=ax4)
        ax6 = fig.add_subplot(326, sharex=ax4)
        for ax_ in (ax1, ax2, ax3, ax4, ax5, ax6):
            ax_.grid(zorder=-10)
        ax4.set_title("Waveform")
        # if pre_filt is not None:
        # #     text = 'pre_filt: None'
        # # else:

        output_dict = {'DISP':'Displacement', 'VEL':'Velocity',
             'ACC':'Acceleration'}

        # labels - specific to Apollo data
        unit_dict = {'DISP':'DU/m', 'VEL':'DU/(m/s)',
             'ACC':r'DU/(m/s$\mathrm{^2}$)'}
        unit_inverse_dict = {'DISP':'m/DU', 'VEL':'(m/s)/DU',
             'ACC':r'(m/s$\mathrm{^2}$)/DU'}
        unit_final_dict = {'DISP':'m', 'VEL':'m/s',
             'ACC':r'm/s$\mathrm{^2}$'}
        raw_unit = 'DU'


        ax1.set_ylabel("Raw [{}]".format(raw_unit), bbox=bbox1, fontsize=16)
        ax1.yaxis.set_label_coords(-0.13, 0.5)
        if pre_filt is not None:
            text = 'Pre-filter parameters:\n [{:.3g}, {:.3g}, {:.3g}, {:.3g}]'.format(
                *pre_filt)
            ax1.text(0.05, 0.1, text, ha="left", va="bottom",
                     transform=ax1.transAxes, fontsize="large", bbox=bbox,
                     zorder=5)
            ax1b.set_ylabel("Pre-filter", bbox=bbox2,
                     fontsize=16)
            ax1b.yaxis.set_label_coords(1.14, 0.5)
        ax1.set_title("Spectrum",color='b')
        output_dict = {'DISP':'Displacement', 'VEL':'Velocity',
             'ACC':'Acceleration'}
        evalresp_info = "\n".join(
            ['Output: %s' % output_dict[output]] +
            ['%s: %s' % (key, value) for key, value in kwargs.items()])
        ax2.text(0.05, 0.1, evalresp_info, ha="left",
                 va="bottom", transform=ax2.transAxes,
                 fontsize="large", zorder=5, bbox=bbox)
        ax2.set_ylabel("Pre-processed [{}]".format(raw_unit), bbox=bbox1,
                 fontsize=16)
        ax2.yaxis.set_label_coords(-0.13, 0.5)
        ax2b.set_ylabel("Instrument\nresponse [{}]".format(unit_dict[output]), 
                 bbox=bbox2, fontsize=16)

        ax2b.yaxis.set_label_coords(1.14, 0.5)
        if water_level is not None:
            ax3.text(0.05, 0.1, "Water Level: %s" % water_level,
                     ha="left", va="bottom", transform=ax3.transAxes,
                     fontsize="large", zorder=5, bbox=bbox)
        ax3.set_ylabel("Multiplied with inverted\n"
                       "instrument response [{}]".format(unit_final_dict[output]),
                       bbox=bbox1, fontsize=16)
        ax3.yaxis.set_label_coords(-0.13, 0.5)
        if water_level is not None:
            ax3b.set_ylabel("Inverted instrument response,\n"
                            "water level applied [{}]".format(
                            unit_inverse_dict[output]), bbox=bbox2, fontsize=16)
        else:
            ax3b.set_ylabel("Inverted instrument\nresponse [{}]".format(
                            unit_inverse_dict[output]), bbox=bbox2, fontsize=16)
        ax3b.yaxis.set_label_coords(1.14, 0.5)


        # xlbl = ax3b.yaxis.get_label()
        # print(xlbl.get_position())
        #
        #

        # xlbl = ax3b.yaxis.get_label()
        # print(xlbl.get_position())
        ax3.set_xlabel("Frequency [Hz]")
        times = trace.times()
        ax4.plot(times, trace.data, color="k")
        ax4.set_ylabel("Raw [{}]".format(raw_unit),fontsize=16)
        ax4.yaxis.set_ticks_position("right")
        ax4.yaxis.set_label_position("right")
        ax4.yaxis.set_label_coords(1.14, 0.5)
        # Ceri
        # ax5.plot(times, data, color="k")

        ax5.set_ylabel("Pre-processed [{}]".format(raw_unit), fontsize=16)
        ax5.yaxis.set_ticks_position("right")
        ax5.yaxis.set_label_position("right")
        ax5.yaxis.set_label_coords(1.14, 0.5)
        ax6.set_ylabel("Response removed [{}]".format(unit_final_dict[output]),
             fontsize=16)
        ax6.set_xlabel("Time [s]")
        ax6.yaxis.set_ticks_position("right")
        ax6.yaxis.set_label_position("right")
        ax6.yaxis.set_label_coords(1.14, 0.5)

    # smart calculation of nfft dodging large primes
    from obspy.signal.util import _npts2nfft
    nfft = _npts2nfft(npts)
    # Transform data to Frequency domain
    data = np.fft.rfft(data, n=nfft)
    # calculate and apply frequency response,
    # optionally prefilter in frequency domain and/or apply water level
    freq_response, freqs = \
        response.get_evalresp_response(trace.stats.delta, nfft,
                                       output=output, **kwargs)

    if plot:
        ax1.loglog(freqs, np.abs(data), color=color1, zorder=9)

    # frequency domain pre-filtering of data spectrum
    # (apply cosine taper in frequency domain)
    if pre_filt:
        freq_domain_taper = cosine_sac_taper(freqs, flimit=pre_filt)
        data *= freq_domain_taper

    if plot:
        try:
            freq_domain_taper
        except NameError:
            freq_domain_taper = np.ones(len(freqs))
        if pre_filt is not None:
            ax1b.semilogx(freqs, freq_domain_taper, color=color2, zorder=10)
            ax1b.set_ylim(-0.05, 1.05)
        ax2.loglog(freqs, np.abs(data), color=color1, zorder=9)
        ax2b.loglog(freqs, np.abs(freq_response), color=color2, zorder=10)
        # Ceri - transform data back into the time domain - just to view it
        filt_data = np.fft.irfft(data)[0:npts]
        ax5.plot(times, filt_data, color="k")


    if water_level is None:
        # No water level used, so just directly invert the response.
        # First entry is at zero frequency and value is zero, too.
        # Just do not invert the first value (and set to 0 to make sure).
        freq_response[0] = 0.0
        freq_response[1:] = 1.0 / freq_response[1:]
    else:
        # Invert spectrum with specified water level.
        invert_spectrum(freq_response, water_level)

    data *= freq_response
    data[-1] = abs(data[-1]) + 0.0j

    spectra = np.abs(data)
    if plot:
        ax3.loglog(freqs, np.abs(data), color=color1, zorder=9)
        ax3b.loglog(freqs, np.abs(freq_response), color=color2, zorder=10)

    # transform data back into the time domain
    data = np.fft.irfft(data)[0:npts]

    if plot:
        # Oftentimes raises NumPy warnings which we don't want to see.
        with np.errstate(all="ignore"):
            ax6.plot(times, data, color="k")
            plt.subplots_adjust(wspace=0.4)
            if plot is True and fig is None:
                plt.show()
            elif plot is True and fig is not None:
                # Ceri - changed this so that the graph does actually show
                plt.show()
            else:
                plt.savefig(plot)
                plt.close(fig)

    # assign processed data and store processing information
    trace.data = data
    if return_spectra:
        return trace, freqs, spectra
    else: 
        return trace 


def process_remove_response_Apollo(stream, inv, pre_filt=None, output='VEL'):
    # see also remove_response_from_seismogram() in extra_plots
    

    stream.attach_response(inv)

    # remove the response when we need to remove the mean and deal with 
    # gaps in the seismogram
    # stream = stream.merge()


    # remove location (ground station)
    for tr in stream:
        tr.stats.location = ''

    stream = stream.split()
    # detrend
    stream.detrend('linear')

    # merge the streams
    stream.merge()
    # if stream.count() > 1:
    #     print('Too many streams - exiting')


    for tr in stream:

        # find the gaps in the trace
        if isinstance(tr.data,np.ma.MaskedArray):
            mask = np.ma.getmask(tr.data)
            # refill the gaps in the trace with zeros 
            tr.data = np.ma.filled(tr.data,fill_value=0)
        else:
            mask = None

    stream = stream.merge()
    #     # remove the mean
    #     mean = tr.data.mean()
    #     tr.data = tr.data - mean
        

        # # find the gaps in the trace
        # if isinstance(tr.data,np.ma.MaskedArray):
        #     mask = np.ma.getmask(tr.data)
        #     # refill the gaps in the trace with zeros 
        #     tr.data = np.ma.filled(tr.data,fill_value=0)
        #     print('about to plot')
        #     tr.plot()
        # else:
        #     mask = None

        # zero_mean=False - because the trace can be asymmetric - remove the mean ourselves
        # do not taper here - it doesn't work well with the masked arrays - often required
        # when there are gaps - if necessary taper first
        # water level - this probably doesn't have much impact - because we are pre filtering
        # stream.remove_response(pre_filt=pre_filt,output="DISP",water_level=30,zero_mean=False,taper=False,plot=True,fig=outfile)
        # remove_response(tr, pre_filt=pre_filt,output=output,
        #         water_level=None,zero_mean=False,taper=False,plot=False)

    # is there an issue where this is not returning the updated trace? 
    trace, freqs, spectra =remove_response(stream[0], pre_filt=pre_filt,output=output,
            water_level=None,zero_mean=False,taper=False,plot=False,return_spectra=True)
    # remove_response(stream[0], pre_filt=pre_filt,output=output,
    #         water_level=None,zero_mean=False,taper=False,plot=True,return_spectra=False)
    stream[0] = trace

        # plt.loglog(freqs, spectra, alpha=0.5)
        # plt.show()
        # tr.spectrogram(log=False,wlen=100,dbscale=False,
        # clip=[0.85,1.0])
        #   title='{} S12 '.format(str(starttime)))
        # exit()

    # now put the masked areas back 
    if mask is not None:
        stream[0].data = np.ma.array(stream[0].data, mask = mask)  

    return stream 





def process_envelope(stream,inv,pre_filt,output='VEL',smooth_periods=10,square=True):

    stream = process_remove_response_Apollo(stream, inv, pre_filt=pre_filt,output=output)
    

    # stream.attach_response(inv)
    # print('len ', len(stream), square)
    # 
    # # remove the response when we need to remove the mean and deal with 
    # # gaps in the seismogram
    # stream = stream.merge()
    # 
    # for tr in stream:
    #     # remove the mean
    #     mean = tr.data.mean()
    #     tr.data = tr.data - mean
    # 
    #     # find the gaps in the trace
    #     if isinstance(tr.data,np.ma.MaskedArray):
    #         mask = np.ma.getmask(tr.data)
    #         # refill the gaps in the trace with zeros 
    #         tr.data = np.ma.MaskedArray(tr.data,fill_value=0)
    #     else:
    #         mask = None
    # 
    #     # zero_mean=False - because the trace can be asymmetric - remove the mean ourselves
    #     # do not taper here - it doesn't work well with the masked arrays - often required
    #     # when there are gaps - if necessary taper first
    #     # water level - this probably doesn't have much impact - because we are pre filtering
    #     # stream.remove_response(pre_filt=pre_filt,output="DISP",water_level=30,zero_mean=False,taper=False,plot=True,fig=outfile)
    #     remove_response(tr, pre_filt=pre_filt,output=output,
    #             water_level=None,zero_mean=False,taper=False,plot=False)
    for tr in stream:
        # find the gaps in the trace
        if isinstance(tr.data,np.ma.MaskedArray):
            mask = np.ma.getmask(tr.data)
            # refill the gaps in the trace with zeros 
            tr.data = np.ma.MaskedArray(tr.data,fill_value=0)
        else:
            mask = None

        if square: 
            tr.data = np.square(tr.data)

        # find the envelope 
        tr.data=envelope(tr.data)
        # calculate the smoothing length (no. of samples)
        if pre_filt is not None: 
            mid_period=1/(pre_filt[1] + (pre_filt[2] - pre_filt[1])/2)
            smooth_length = math.ceil(tr.stats.sampling_rate * mid_period * smooth_periods)
        else: 
            smooth_length = 100
    
        y_avg = np.zeros((1, len(tr)))
        avg_mask = np.ones(smooth_length)/smooth_length
        y_avg = np.convolve(tr.data, avg_mask, 'same')
        
        tr.data = y_avg

        # now put the masked areas back 
        if mask is not None:
            tr.data = np.ma.array(tr.data, mask = mask) 

    return stream 

#             SECONDS_PER_DAY = 86400
#             for i, ax in enumerate(fig.get_axes()):
#                 for ii, in enumerate([(0.3,0.9)])
# 
#                 y_avg = np.zeros((len(window_lst) , len(stream[i])))
#                 for ii, window in enumerate(window_lst):
#                     avg_mask = np.ones(window) / window
#                     x_time = ((stream[ii].times() / SECONDS_PER_DAY)) + date2num(stream[ii].stats.starttime.datetime)
#                     y_avg = np.convolve(stream[i].data, avg_mask, 'same')
#                     # y_avg[i, :] = np.convolve(stream[ii].data, avg_mask, 'same')
# 
#                 	# Plot each running average with an offset of 50
#                 	# in order to be able to distinguish them
#                     # print(y_avg[i, :].shape)
#                     # print(x_time.shape)
#                     ax.plot(x_time, y_avg, label=window, color='r')





    # for tr in stream:
    #     if tr.stats.channel == 'MH1':
    #         tr_envelope.stats.channel = 'EN1'
    #         if mask_1 is not None:
    #             tr_envelope.data = np.ma.array(tr_envelope.data, mask = mask_1)

        

        # # find the gaps in the trace
        # if isinstance(tr.data,np.ma.MaskedArray):
        #     mask = np.ma.getmask(tr.data)
        #     # refill the gaps in the trace with zeros 
        #     tr.data = np.ma.MaskedArray(tr.data,fill_value=0)
        # else:
        #     mask = None
        # 
        # # now put the masked areas back 
        # if mask is not None:
        #     tr.data = np.ma.array(tr.data, mask = mask)  

def plot_seismogram_with_envelopes(stream, inv, station_code, catalog, file_out, event=None, origin=None, picks=None, 
    begintime=None, output='VEL', pre_filt_main=None, 
    pre_filt_env=None,smooth_periods=10,time_before=600,time_after=3600,plot_type='with_seismogram',title=None):

    if len(stream) > 0: 
        station = stream[0].stats.station
        channel = stream[0].stats.channel

        if title:
            title = '{} - {} - {}'.format(title, channel, station)

        if origin is not None: 
            latitude = LATITUDES[STATION_LIST.index(station_code)]
            longitude = LONGITUDES[STATION_LIST.index(station_code)]
            MH1_azimuth = MH1_AZIMUTHS[STATION_LIST.index(station_code)]
            MH2_azimuth = MH2_AZIMUTHS[STATION_LIST.index(station_code)]

            distance, azimuth_A_B, azimuth_B_A =  gps2dist_azimuth(
              origin.latitude, origin.longitude, latitude, longitude, MOON_RADIUS, MOON_FLATTENING)
            distance = distance/1000.
            print(event.event_type)
            for desc in event.event_descriptions:
                print(desc.text)
            print('{}-{}'.format(station_code, channel))
            print('Distance: {:.1f} km, Azimuth: {:.1f} deg, Back Azimuth: {:.1f} deg'.format(distance, azimuth_A_B, azimuth_B_A))


        stream_orig = stream.copy()
        stream = process_remove_response_Apollo(stream, inv, pre_filt=pre_filt_main,output=output)
    
        if plot_type=='with_seismogram':
    
            fig = plt.figure(figsize=(8,3))
            for i, tr in enumerate(stream): 
                x_time = tr.times()-time_before
                plt.plot(x_time,tr.data,color='k')

                if pre_filt_env is not None:
                    for ii, pre_filt in enumerate(pre_filt_env):
                        stream_env = stream_orig.copy()
                        stream_env = process_envelope(stream_env, inv, pre_filt=pre_filt,output=output,smooth_periods=smooth_periods,square=False)
                        # print('before')
                        # stream_env.plot()
                        # plot the envelope  
                        plt.plot(x_time,stream_env[0].data,label='{} - {} Hz'.format(pre_filt[1],pre_filt[2]),color=colors[ii])
                        # print('after')
                        # stream_env.plot()

            if title:
                plt.title(title)
            plt.xlabel('Seconds',fontsize='12')
            plt.ylabel(r'Acceleration [m/s$^{-2}$]',fontsize='12')
            plt.legend(loc='lower right',framealpha=0)
            plt.xlim(-time_before,time_after)
            plt.subplots_adjust(bottom=0.15)

            if picks is not None: 
                for pick in picks:
                    if pick.phase_hint == 'P' and pick.waveform_id.station_code==station_code:
                        pick_markP = pick.time - begintime - time_before
                        plt.gca().axvline(x=pick_markP, 
                                              color=pick_colors[0], linewidth=2)


                    elif pick.phase_hint == 'S' and pick.waveform_id.station_code==station_code:
                        pick_markS = pick.time - begintime - time_before
                        plt.gca().axvline(x=pick_markS, 
                                              color=pick_colors[1], linewidth=2)

                    elif (pick.filter_id is not None and 
                      pick.phase_hint in ('t_max','t_onset') and 
                      pick.waveform_id.station_code==station_code): 

                        pick_mark = pick.time - begintime - time_before
                        if pick.filter_id == filter_ids[0]:
                            color=colors[0]
                        elif pick.filter_id == filter_ids[1]:
                            color=colors[1]
                        elif pick.filter_id == filter_ids[2]:
                            color=colors[2]

                        plt.gca().axvline(x=pick_mark, 
                            color=color, linewidth=2)

                    elif (pick.filter_id is not None and 
                      pick.phase_hint in ('t_max_lower','t_max_upper') and 
                      pick.waveform_id.station_code==station_code): 

                        pick_mark = pick.time - begintime - time_before
                        if pick.filter_id == filter_ids[0]:
                            color=colors[0]
                        elif pick.filter_id == filter_ids[1]:
                            color=colors[1]
                        elif pick.filter_id == filter_ids[2]:
                            color=colors[2]

                        plt.gca().axvline(x=pick_mark, 
                            color=color, linewidth=1, linestyle='dashed')


            # XXXXXX

            plt.show(block=False)
            while True:
                print('Left click on t_onset or t_max. Or press return.')
                cid = fig.canvas.mpl_connect('button_press_event', onclick)
                input_string = '''Left click on t_onset or t_max.
Set t_onset and t_max: (to1), (tm1), (to2), (tm2), (to3), (tm3)
Set t_max_lower or t_max_upper: (tl1), (tl2), (tl3), (tu1), (tu2), (tu3)
Or press return.'''
                output_str = input(input_string)
                if output_str in ('to1', 'tm1', 'to2', 'tm2', 'to3', 'tm3',
                      'tl1', 'tl2', 'tl3', 'tu1', 'tu2', 'tu3'):
                    measurement_type = output_str[0:2]
                    number = int(output_str[2])
                    freq = pre_filt_env[number-1]
                    if output_str in ('tl1', 'tl2', 'tl3', 'tu1', 'tu2', 'tu3'):
                        plt.gca().axvline(x=xdata,
                          color=colors[number-1], linewidth=1, linestyle='dashed')
                    else: 
                        plt.gca().axvline(x=xdata,
                          color=colors[number-1], linewidth=3)
                    plt.draw()
                    save_pick_to_event(event,output_str,xdata,station_code,begintime,
                      time_before,pre_filt_env)
                    # plt.pause(0.01)

                elif output_str in ('\n'):
                    break

                # this is not written yet, but it's to save out the changes
                # to the xml file 
                # if quit == False:
                #     res_id = str(pick.resource_id).replace('pick','amplitude_est')
                #     res_id = ResourceIdentifier(id=res_id)
                #     res_id.convert_id_to_quakeml_uri(authority_id=AUTHORITY_ID)
                #     amplitude = Amplitude(resource_id=res_id)
                # 
                #     if output_str not in ('n', 'na'):
                #         amplitude.unit = "m/(s*s)"
                #         amplitude.generic_amplitude = ampl_dict[output_str]
                #         amplitude.type = 'A'
                #     else: 
                #         amplitude.type = ampl_dict[output_str]
                #     amplitude.pick_id = pick.resource_id
                #     ev.amplitudes.append(amplitude)
                # 
                #     catalog.write(filename=file_out, format='QUAKEML') 
                #     print(i+start_no)
                # else:
                #     break

        
        elif plot_type=='normalized':
            # TODO check whether it is too curved. 
            print("""This is a squared version of the envelope (for intensity).
                But it is not clear whether the start is too curved.""")
            fig, axes = plt.subplots(1, 1, sharex=True, figsize=(7, 4))
            for i, tr in enumerate(stream): 
                x_time = tr.times()
                if pre_filt_env is not None:
                    for ii, pre_filt in enumerate(pre_filt_env):
                        stream_env = stream_orig.copy()
                        # calculate normalization factor 
                        # smooth this strongly to help with the normalization 
                        process_envelope(stream_env, inv, pre_filt=pre_filt,output=output,smooth_periods=400,square=True)
                        trace_max = stream_env[0].max()
                        stream_env = stream_orig.copy()
                        # process again 
                        process_envelope(stream_env, inv, pre_filt=pre_filt,output=output,smooth_periods=smooth_periods,square=True)
                        # don't plot them on top of each other, but shift by 1 each time
                        # plt.plot(x_time,stream_env[0].data/trace_max+ ii,label=str(pre_filt))
                        plt.plot(x_time,stream_env[0].data/trace_max,label='{} - {} Hz'.format(pre_filt[1],pre_filt[2]))


            if title:
                plt.title(title)
            plt.xlabel('Seconds',fontsize='12')
            plt.ylabel(r'Intensity [m$^2$/s$^{-4}$]',fontsize='12')
            plt.legend()
            plt.subplots_adjust(bottom=0.12)
            plt.ylim(0,2)

# yyyyyy
def save_pick_to_event(event,output_str,pick_seconds,station_code,begintime,
  time_before,pre_filt_env):

    measurement_type = output_str[0:2]
    if measurement_type == 'to':
        phase_hint = 't_onset'
    elif measurement_type == 'tm':
        phase_hint = 't_max'
    elif measurement_type == 'tl':
        phase_hint = 't_max_lower'
    elif measurement_type == 'tu':
        phase_hint = 't_max_upper'

    number = int(output_str[2])
    freq = pre_filt_env[number-1]
    filter_id = filter_ids[number-1] 

    found = False
    for pick in event.picks:
        if (pick.phase_hint == phase_hint and 
          pick.waveform_id.station_code==station_code and 
          pick.filter_id==filter_id): 
            found = True
            pick.time = begintime+time_before+pick_seconds

    if found is False: 
        pick = Pick()
        pick.phase_hint = phase_hint
        pick.time = begintime+time_before+pick_seconds
        pick.filter_id = filter_id
        waveformStreamID = WaveformStreamID()
        waveformStreamID.station_code = station_code
        waveformStreamID.network_code = 'XA'
        pick.waveform_id = waveformStreamID
        event.picks.append(pick)
    # 
    # stream += stream_orig


    # # automerge required to avoid changing the order of the plot, which we need
    # fig = stream.plot(handle=True,automerge=False,show=False,size=(1200,600))
    # 
    # window_lst = [20]
    # 
    # 
    # for i, ax in enumerate(fig.get_axes()):
    #     # for ii, in enumerate([(0.3,0.9)]):
    #         y_avg = np.zeros((len(window_lst) , len(stream[i])))
    #         for iii, window in enumerate(window_lst):
    #             avg_mask = np.ones(window) / window
    #             x_time = ((stream[iii].times() / SECONDS_PER_DAY)) + date2num(stream[iii].stats.starttime.datetime)
    #             y_avg = np.convolve(stream[i].data, avg_mask, 'same')
    #             # y_avg[i, :] = np.convolve(stream[iii].data, avg_mask, 'same')
    # 
    #         	# Plot each running average with an offset of 50
    #         	# in order to be able to distinguish them
    #             # print(y_avg[i, :].shape)
    #             # print(x_time.shape)
    #             ax.plot(x_time, y_avg, label=window, color='r')
    # plt.show()

    # for tr in stream_orig:
    #     tr_envelope = tr.copy()
    #     tr_envelope.data=envelope(tr.data)
    #             if tr.stats.channel == 'MH1':
    #                 tr_envelope.stats.channel = 'EN1'
    #                 if mask_1 is not None:
    #                     tr_envelope.data = np.ma.array(tr_envelope.data, mask = mask_1)
    #             elif tr.stats.channel == 'MH2':
    #                 tr_envelope.stats.channel = 'EN2'
    #                 if mask_2 is not None:
    #                     tr_envelope.data = np.ma.array(tr_envelope.data, mask = mask_2)
    #             elif tr.stats.channel == 'MHZ':
    #                 tr_envelope.stats.channel = 'ENZ'
    #                 if mask_Z is not None:
    #                     tr_envelope.data = np.ma.array(tr_envelope.data, mask = mask_Z)


    

# def plot_remove_response(stream, inv, output='ACC'):
# 
#             stream.attach_response(inv)
# 
#             # pre_filt = [0.07, 0.1,0.3,0.4]
#             pre_filt = None
#             pre_filt = [0.2,0.3,0.9,1.3]
# 
# 
#             for tr in stream:
#                 # remove the mean
#                 mean = tr.data.mean()
#                 tr.data = tr.data - mean
# 
#                 # find the gaps in the trace
#                 if isinstance(tr.data,np.ma.MaskedArray):
#                     if tr.stats.channel == 'MH1':
#                         mask_1 = np.ma.getmask(tr.data)
#                     elif tr.stats.channel == 'MH2':
#                         mask_2 = np.ma.getmask(tr.data)
#                     elif tr.stats.channel == 'MHZ':
#                         mask_Z = np.ma.getmask(tr.data)
#                 else:
#                     if tr.stats.channel == 'MH1':
#                         mask_1 = None
#                     elif tr.stats.channel == 'MH2':
#                         mask_2 = None
#                     elif tr.stats.channel == 'MHZ':
#                         mask_Z = None
# 
#             # split the stream, then refill it with zeros on the gaps
#             stream = stream.split()
#             stream = stream.merge(fill_value=0)
# 
#             # zero_mean=False - because the trace can be asymmetric - remove the mean ourselves
#             # do not taper here - it doesn't work well with the masked arrays - often required
#             # when there are gaps - if necessary taper first
#             # water level - this probably doesn't have much impact - because we are pre filtering
#             # stream.remove_response(pre_filt=pre_filt,output="DISP",water_level=30,zero_mean=False,taper=False,plot=True,fig=outfile)
#             for tr in stream:
#                 if tr.stats.channel == 'MHZ':
#                     trace, freqs, spectra = remove_response(tr, pre_filt=pre_filt,output=output,
#                             water_level=None,zero_mean=False,taper=False,plot=False,return_spectra=True)
#                     # tr.spectrogram(log=True,wlen=100,dbscale=False,
#                     #   # clip=[0.85,1.0],
#                     #   title='{} S12 '.format(str(starttime)))
#                     # plt.loglog(freqs, spectra, color='blue')
#                     # plt.show()
# 
#                 else:
#                     remove_response(tr, pre_filt=pre_filt,output=output,
#                             water_level=None,zero_mean=False,taper=False,plot=False)
# 
# 
#                 # from scipy.signal import hilbert, chirp
#                 # fs = 400.0
#                 # duration = 1.0
#                 # samples = int(fs*duration)
#                 # t = np.arange(samples) / fs
#                 # 
#                 # tr.data = chirp(t, 20.0, t[-1], 100.0)
#                 # tr.data *= (1.0 + 0.5 * np.sin(2.0*np.pi*3.0*t) )
# 
#                 # stream.append(tr_envelope)
# 
#                 if tr.stats.channel == 'MH1':
#                     if mask_1 is not None:
#                         tr.data = np.ma.array(tr.data, mask = mask_1)
#                 elif tr.stats.channel == 'MH2':
#                     if mask_2 is not None:
#                         tr.data = np.ma.array(tr.data, mask = mask_2)
#                 elif tr.stats.channel == 'MHZ':
#                     if mask_Z is not None:
#                         tr.data = np.ma.array(tr.data, mask = mask_Z)
# 
# 
# 
# 
#             # for pick in event.picks:
#             #     # print(pick.time)
#             #     pick_mark = (((pick.time - tr.stats.starttime) / SECONDS_PER_DAY) +
#             #       date2num(tr.stats.starttime.datetime))
#             #     plt.axvline(x=pick_mark, color='b')
# 
# 
#             # Compute moving averages using different window sizes
#             # window_lst = [3, 6, 10, 16, 22, 35]
#             window_lst = [3]
# 
#             fig = stream.plot(handle=True,automerge=False,show=False,size=(1200,600))
#                       #     pick_mark = (((pick.time - tr.stats.starttime) / SECONDS_PER_DAY) +
#                         #       date2num(tr.stats.starttime.datetime))
# 
#                 # add the envelope function (before the mask?)
#                 tr_envelope = tr.copy()
#                 tr_envelope.data=envelope(tr.data)
#                 if tr.stats.channel == 'MH1':
#                     tr_envelope.stats.channel = 'EN1'
#                     if mask_1 is not None:
#                         tr_envelope.data = np.ma.array(tr_envelope.data, mask = mask_1)
#                 elif tr.stats.channel == 'MH2':
#                     tr_envelope.stats.channel = 'EN2'
#                     if mask_2 is not None:
#                         tr_envelope.data = np.ma.array(tr_envelope.data, mask = mask_2)
#                 elif tr.stats.channel == 'MHZ':
#                     tr_envelope.stats.channel = 'ENZ'
#                     if mask_Z is not None:
#                         tr_envelope.data = np.ma.array(tr_envelope.data, mask = mask_Z)
# 
#             SECONDS_PER_DAY = 86400
#             for i, ax in enumerate(fig.get_axes()):
#                 for ii, in enumerate([(0.3,0.9)])
# 
#                 y_avg = np.zeros((len(window_lst) , len(stream[i])))
#                 for ii, window in enumerate(window_lst):
#                     avg_mask = np.ones(window) / window
#                     x_time = ((stream[ii].times() / SECONDS_PER_DAY)) + date2num(stream[ii].stats.starttime.datetime)
#                     y_avg = np.convolve(stream[i].data, avg_mask, 'same')
#                     # y_avg[i, :] = np.convolve(stream[ii].data, avg_mask, 'same')
# 
#                 	# Plot each running average with an offset of 50
#                 	# in order to be able to distinguish them
#                     # print(y_avg[i, :].shape)
#                     # print(x_time.shape)
#                     ax.plot(x_time, y_avg, label=window, color='r')
# 
#                 # Add legend to plot
#                 # ax.legend()
#                 # ax.axhline(-ampl_dict['purple'], color='purple')
#                 # ax.axhline(ampl_dict['purple'], color='purple')
#                 # ax.axhline(-ampl_dict['pink'], color='pink')
#                 # ax.axhline(ampl_dict['pink'], color='pink')
#                 # ax.axhline(-ampl_dict['darkblue'], color='darkblue')
#                 # ax.axhline(ampl_dict['darkblue'], color='darkblue')
#                 # ax.axhline(-ampl_dict['lightblue'], color='lightblue')
#                 # ax.axhline(ampl_dict['lightblue'], color='lightblue')
#                 # ax.axhline(-ampl_dict['green'], color='green')
#                 # ax.axhline(ampl_dict['green'], color='green')
#                 # ax.axhline(-ampl_dict['yellow'], color='yellow')
#                 # ax.axhline(ampl_dict['yellow'], color='yellow')
#                 # ax.axhline(-ampl_dict['orange'], color='orange')
#                 # ax.axhline(ampl_dict['orange'], color='orange')
#                 # ax.axhline(-ampl_dict['red'], color='red')
#                 # ax.axhline(ampl_dict['red'], color='red')
#                 # ax.set_ylim(-17e-9,17e-9)
#                 # ax.ylabel=''
#             plt.show()

def plot_example():
    np.random.seed(17)
    n = 10
    N = 500
    x = np.linspace(0, n, N)
    y0 = -0.05*x**4 + 5*x**2 + 7*x - 6
    yn = 4.5*np.random.standard_normal(N) * np.log10(y0**2 + 0.1)/2
    y = y0 + yn


    # Create a figure canvas and plot the original, noisy data
    fig, ax = plt.subplots()
    ax.plot(x, y, label="Original")
    # Compute moving averages using different window sizes
    window_lst = [3, 6, 10, 16, 22, 35]
    y_avg = np.zeros((len(window_lst) , N))
    for i, window in enumerate(window_lst):
    	avg_mask = np.ones(window) / window
    	y_avg[i, :] = np.convolve(y, avg_mask, 'same')
    	# Plot each running average with an offset of 50
    	# in order to be able to distinguish them
    	ax.plot(x, y_avg[i, :] + (i+1)*50, label=window)
    # Add legend to plot
    ax.legend()
    plt.show()

###########HERE
def view_catalog_with_envelopes(top_level_dir,file,inv_name,
  dir_type='processed_dir',pre_filt_main=None,pre_filt_env=None,output='VEL',
  smooth_periods=10,time_before=600,time_after=3600,plot_type='with_seismogram', 
  stations=STATION_LIST,channel='MHZ',start_no=0,end_no=None,title=None):
    """
    :param top_level_dir: string
        full path the directory 
    :param file: 
        
    :param inv_name: 
    :param dir_type: 
    :param pre_filt_main: 
    :param pre_filt_env: 
    :param output: default='VEL'
    :param top_level_dir: 
    :param top_level_dir: 
    """

    file_out = file.replace('.xml', '_out.xml')

    # read the response file
    inv = read_inventory(inv_name)
    
    # print(inv)
    # for network in inv:
    #     for station in network[0]:
    #         print(station)
    # print(inv[0][0][0].response)
    
    catalog = read_events(file)
    quit = False
    
    if end_no is None:
        end_no = len(catalog) + 1
    
    for i, ev in enumerate(catalog[start_no:end_no]):
        print('Event: ', start_no+i, ev)
        quit=True

        for station_code in stations: 
            channel = 'MHZ'
            while True: 
            
                found = display_envelopes(top_level_dir=top_level_dir, event=ev, 
                    event_no=start_no+i, 
                    inv=inv, catalog=catalog, file_out=file_out,
                    station_code=station_code, channel=channel, dir_type=dir_type,
                    pre_filt_main=pre_filt_main,output=output,
                    pre_filt_env=pre_filt_env, 
                    smooth_periods=smooth_periods,time_before=time_before,time_after=time_after,
                    plot_type=plot_type,title=title)



                if found == False: #next
                    if i+start_no+1 == end_no: 
                        quit=True
                    else: 
                        set_measurement = None
                        quit=False
                    break 

                plt.close()
                input_string = '''View again (v), next (n), save (s), MH1, MH2, MHZ, quit(q)'''
# \nOr set t_onset and t_max: (to1), (tm1), (to2), (tm2), (to3), (tm3)\n'''
                output_str = input(input_string)
                if output_str == 'v': #view again 
                    set_measurement = None
                elif output_str in ('MH1', 'MH2', 'MHZ'):
                    channel=output_str
                    set_measurement = None
                elif output_str == 's':
                    catalog.write(filename=file_out, format='QUAKEML') 
                    if i+start_no+1 == end_no: 
                        quit=True
                    else: 
                        set_measurement = None
                        quit=False
                    break 
                elif output_str == 'q': #quit
                    quit=True
                    break
                elif output_str == 'n': #next
                    if i+start_no+1 == end_no: 
                        quit=True
                    else: 
                        set_measurement = None
                        quit=False
                    break 

            # break from the station loop 
            if quit:
                break
            else: 
                continue 
        # break from the event loop 
        if quit:
            break

def display_envelopes(top_level_dir, event, event_no, inv, catalog, file_out, station_code, channel,
      dir_type='pdart_dir',pre_filt_main=None,pre_filt_env=None,
      output='VEL', 
      smooth_periods=10,time_before=600,time_after=3600,plot_type='with_seismogram',title=None):

    # print(top_level_dir, 
    # print(event_no)
    #      # event, event_no, inv, station_code, channel,
    #      #  dir_type,pre_filt_main,pre_filt_env,
    #      #  output, 
    #      #  smooth_periods,time_before,time_after,plot_type)
    # exit()

    # use the preferred origin if it exists
    origin = event.preferred_origin()
    if origin is None: 
        try: 
            origin = event.origins[0]
        except: 
            origin = None

    endtime = None
    starttime = None
    picks = event.picks
    if origin is not None and origin.time is not None: 
        starttime  = origin.time
    else: 
        
        for pick in picks:
            if pick.waveform_id.station_code == station_code: 
                print('Event in catalog number {}, Pick Time {}, Phase hint {}'.format(i+start_no, pick.time, pick.phase_hint))
                # continue
                if not starttime:
                    starttime = pick.time
                else:
                    if pick.time < starttime:
                        starttime = pick.time

        # if starttime is not None, then a pick was found for this event and station
        if starttime is None:
            print('Pick not found at station {} \n for Event: {}'.format(station_code,ev))
            return False


    # print('Temporarily changed to 1/2 hour')
    # starttime -= 1800.
    print('T ', time_before, time_after, starttime) 
    begintime = starttime - time_before
    endtime = starttime + time_after

    stream = find_seismogram(top_level_dir,begintime,endtime,
      stations=[station_code],dir_type=dir_type,channels=[channel])

    if stream is None or len(stream) == 0: 
        print('Event at {} not found for station {}'.format(str(starttime),station_code))
        return False
    # else: 
    #     print('Event found, ', event)

    plot_seismogram_with_envelopes(stream, inv, station_code=station_code,
         catalog=catalog,
         file_out=file_out,
         event=event,
         origin=origin,
         picks=picks,
         begintime=begintime,
         output=output, 
         pre_filt_main=pre_filt_main, pre_filt_env=pre_filt_env,
         smooth_periods=smooth_periods,
         time_before=time_before,time_after=time_after,
         plot_type=plot_type,title=str(starttime))


# 
# 
#                 # this is not written yet, but it's to save out the changes
#                 # to the xml file 
#                 # if quit == False:
#                 #     res_id = str(pick.resource_id).replace('pick','amplitude_est')
#                 #     res_id = ResourceIdentifier(id=res_id)
#                 #     res_id.convert_id_to_quakeml_uri(authority_id=AUTHORITY_ID)
#                 #     amplitude = Amplitude(resource_id=res_id)
#                 # 
#                 #     if output_str not in ('n', 'na'):
#                 #         amplitude.unit = "m/(s*s)"
#                 #         amplitude.generic_amplitude = ampl_dict[output_str]
#                 #         amplitude.type = 'A'
#                 #     else: 
#                 #         amplitude.type = ampl_dict[output_str]
#                 #     amplitude.pick_id = pick.resource_id
#                 #     ev.amplitudes.append(amplitude)
#                 # 
#                 #     catalog.write(filename=file_out, format='QUAKEML') 
#                 #     print(i+start_no)
#                 # else:
#                 #     break
# 
# 
#             # station level 
# 
#     # finally: 
#     # 
#     #     catalog.write(filename=file_out, format='QUAKEML') 
#     #     print(i+start_no)
# 
# # category = other
# # unit = m/(s*s)
# # genericAmplitude = number
# ###########HERE1

def plot_seismograms_with_envelopes(stream, inv, output='VEL', pre_filt_main=None, 
    pre_filt_env=None,smooth_periods=10,plot_type='section',title=None):

    if len(stream) > 0:

        stream_orig = stream.copy()
        stream = process_remove_response_Apollo(stream, inv, pre_filt=pre_filt_main,output=output)
    
        if plot_type=='with_seismogram':
            fig = plt.figure(figsize=(8,3))
            for i, tr in enumerate(stream): 
                x_time = tr.times()
                plt.plot(x_time,tr.data,color='k')
                if pre_filt_env is not None:
                    for pre_filt in pre_filt_env:
                        stream_env = stream_orig.copy()
                        stream_env = process_envelope(stream_env, inv, pre_filt=pre_filt,output=output,smooth_periods=smooth_periods,square=False)
                        # print('before')
                        # stream_env.plot()
                        # plot the envelope  
                        plt.plot(x_time,stream_env[0].data,label='{} - {} Hz'.format(pre_filt[1],pre_filt[2]))
                        # print('after')
                        # stream_env.plot()

            if title:
                plt.title(title)
            plt.xlabel('Seconds',fontsize='12')
            plt.ylabel(r'Acceleration [m/s$^{-2}$]',fontsize='12')
            plt.legend(loc='lower right',framealpha=0)
            plt.subplots_adjust(bottom=0.15)
            plt.show()


P_degrees =  [0.0, 0.52419354838709675, 1.0483870967741935, 1.5725806451612903, 2.096774193548387, 2.620967741935484, 3.1451612903225805, 3.669354838709677, 4.193548387096774, 4.717741935483871, 5.241935483870968, 5.7661290322580641, 6.290322580645161, 6.814516129032258, 7.3387096774193541, 7.8629032258064511, 8.387096774193548, 8.9112903225806441, 9.435483870967742, 9.9596774193548381, 10.483870967741936, 11.008064516129032, 11.532258064516128, 12.056451612903226, 12.580645161290322, 13.104838709677418, 13.629032258064516, 14.153225806451612, 14.677419354838708, 15.201612903225806, 15.725806451612902, 16.25, 16.25, 16.25, 16.774193548387096, 16.774193548387096, 16.774193548387096, 17.298387096774192, 17.298387096774192, 17.298387096774192, 17.822580645161288, 17.822580645161288, 17.822580645161288, 18.346774193548388, 18.346774193548388, 18.346774193548388, 18.870967741935484, 18.870967741935484, 18.870967741935484, 19.39516129032258, 19.39516129032258, 19.39516129032258, 19.919354838709676, 19.919354838709676, 19.919354838709676, 20.443548387096772, 20.443548387096772, 20.443548387096772, 20.967741935483872, 20.967741935483872, 20.967741935483872, 21.491935483870968, 21.491935483870968, 21.491935483870968, 22.016129032258064, 22.016129032258064, 22.016129032258064, 22.54032258064516, 22.54032258064516, 22.54032258064516, 23.064516129032256, 23.064516129032256, 23.064516129032256, 23.588709677419352, 23.588709677419352, 23.588709677419352, 24.112903225806452, 24.112903225806452, 24.112903225806452, 24.637096774193548, 24.637096774193548, 24.637096774193548, 25.161290322580644, 25.161290322580644, 25.161290322580644, 25.68548387096774, 25.68548387096774, 25.68548387096774, 26.209677419354836, 26.209677419354836, 26.209677419354836, 26.733870967741936, 26.733870967741936, 26.733870967741936, 26.733870967741936, 26.733870967741936, 27.258064516129032, 27.782258064516128, 28.306451612903224, 28.83064516129032, 29.354838709677416, 29.879032258064516, 30.403225806451612, 30.927419354838708, 31.451612903225804, 31.9758064516129, 32.5, 33.024193548387096, 33.548387096774192, 34.072580645161288, 34.596774193548384, 35.12096774193548, 35.645161290322577, 36.169354838709673, 36.693548387096776, 37.217741935483872, 37.741935483870968, 38.266129032258064, 38.79032258064516, 39.314516129032256, 39.838709677419352, 40.362903225806448, 40.887096774193544, 41.411290322580641, 41.935483870967744, 42.45967741935484, 42.983870967741936, 43.508064516129032, 44.032258064516128, 44.556451612903224, 45.08064516129032, 45.604838709677416, 46.129032258064512, 46.653225806451609, 47.177419354838705, 47.701612903225808, 48.225806451612904, 48.75, 49.274193548387096, 49.798387096774192, 50.322580645161288, 50.846774193548384, 51.37096774193548, 51.895161290322577, 52.419354838709673, 52.943548387096769, 53.467741935483872, 53.991935483870968, 54.516129032258064, 55.04032258064516, 55.564516129032256, 56.088709677419352, 56.612903225806448, 57.137096774193544, 57.661290322580641, 58.185483870967737, 58.709677419354833, 59.233870967741936, 59.758064516129032, 60.282258064516128, 60.806451612903224, 61.33064516129032, 61.854838709677416, 62.379032258064512, 62.903225806451609, 63.427419354838705, 63.951612903225801, 64.475806451612897, 65.0]
P_arrivals =  [0.0, 5.2104307435726067, 9.0113759055115565, 12.40855326763106, 15.582958391073161, 18.598421279913364, 21.506314564747264, 24.335240036283619, 27.098534531456746, 29.792860301552015, 32.438745946494699, 35.022899966449906, 37.579035381960757, 40.115878891421204, 42.617665642278546, 45.086005012536724, 47.496660385648482, 49.847376042401898, 52.139920503699962, 54.372177042828831, 56.565161581803373, 58.731461094664624, 60.877709158231575, 63.009910126452873, 65.134078114851533, 67.245948452472177, 69.350240733922291, 71.448642919980799, 73.543084558268433, 75.634450993930557, 77.725626255314978, 79.816465397288141, 80.463270761338052, 80.475781971084103, 81.906834175304439, 82.419322706100161, 82.470226230725828, 83.996702660374297, 84.375219569600503, 84.481446777208191, 86.086237877509831, 86.330882284340944, 86.504163904199302, 88.175316552327004, 88.286233264996739, 88.535358161780039, 90.241195909594325, 90.263250145867147, 90.573120715036836, 92.195685338822699, 92.349834344978817, 92.615890744717234, 94.149581165107719, 94.435384940457226, 94.662726323954985, 96.102814963934634, 96.520452056893774, 96.71277782779282, 98.055318572605628, 98.604249139752639, 98.765470338092769, 100.00698907925388, 100.68593351379361, 100.82032282535017, 101.95774440117115, 102.76594852306434, 102.87692051447968, 103.90749128222271, 104.84525492819786, 104.93497579215786, 105.85617149547005, 106.92318565571422, 106.99419650560002, 107.8036645717745, 108.99970348080942, 109.05437894102371, 109.74990177621025, 111.0754478921578, 111.11533424438443, 111.69479420155892, 113.14937383596144, 113.17689721317051, 113.63824870127065, 115.22166588889469, 115.23893331078929, 115.5801864975126, 117.29156270709669, 117.30129420598615, 117.52051417972703, 119.35982205757462, 119.3638803824136, 119.45914380332265, 121.42636008486512, 121.42646637372907, 121.42655378361536, 121.42659614083379, 121.39599700504812, 123.33099985916378, 125.26403830301368, 127.19503757863623, 129.12391790049128, 131.05093380700407, 132.97672023173845, 134.90092549273982, 136.82336635666115, 138.74393158147487, 140.66250198491906, 142.57897956385929, 144.49327056763784, 146.40534707796928, 148.31507673903644, 150.22240280272879, 152.12726052226787, 154.02955445947919, 155.92652644814805, 157.81930173209136, 159.7082247200228, 161.59334775478985, 163.47470363727058, 165.3522412880717, 167.22591129974884, 169.09568577144739, 170.96146460041734, 172.82329628857499, 174.68102615941447, 176.53461141683903, 178.38401250910957, 180.22914958400267, 182.07157618882164, 183.91190290268682, 185.74952954549656, 187.58413494244743, 189.41552438256596, 191.24355031557377, 193.06808130863826, 194.88901078047431, 196.70530422815315, 198.51544094528114, 200.32052657901977, 202.12068120134691, 203.91597470272512, 205.70639741011161, 207.49188718556334, 209.27309548714521, 211.05022413337724, 212.82296396594677, 214.59111021605204, 216.35438777663282, 218.1117766072885, 219.86366608649573, 221.61013552611507, 223.35125076029107, 225.08695264251779, 226.81719712290615, 228.54283064408381, 230.26372919085509, 231.97955195772536, 233.68984414036851, 235.39060250279124, 237.0838847792393, 238.77023374481377, 240.45091071067884, 242.12905816078427, 243.80376235587647, 245.47445082909744, 247.14079132468859, 248.80249900640999, 250.45739682983688, 252.10012481137173]
S_degrees =  [0.0, 0.52419354838709675, 1.0483870967741935, 1.5725806451612903, 2.096774193548387, 2.620967741935484, 3.1451612903225805, 3.669354838709677, 4.193548387096774, 4.717741935483871, 5.241935483870968, 5.7661290322580641, 6.290322580645161, 6.814516129032258, 7.3387096774193541, 7.8629032258064511, 8.387096774193548, 8.387096774193548, 8.387096774193548, 8.9112903225806441, 9.435483870967742, 9.9596774193548381, 10.483870967741936, 11.008064516129032, 11.532258064516128, 12.056451612903226, 12.580645161290322, 12.580645161290322, 12.580645161290322, 13.104838709677418, 13.629032258064516, 14.153225806451612, 14.677419354838708, 15.201612903225806, 15.725806451612902, 16.25, 16.774193548387096, 16.774193548387096, 16.774193548387096, 17.298387096774192, 17.298387096774192, 17.298387096774192, 17.822580645161288, 17.822580645161288, 17.822580645161288, 18.346774193548388, 18.346774193548388, 18.346774193548388, 18.870967741935484, 18.870967741935484, 18.870967741935484, 19.39516129032258, 19.39516129032258, 19.39516129032258, 19.39516129032258, 19.39516129032258, 19.919354838709676, 19.919354838709676, 19.919354838709676, 20.443548387096772, 20.443548387096772, 20.443548387096772, 20.967741935483872, 20.967741935483872, 20.967741935483872, 21.491935483870968, 21.491935483870968, 21.491935483870968, 21.491935483870968, 21.491935483870968, 22.016129032258064, 22.016129032258064, 22.016129032258064, 22.54032258064516, 22.54032258064516, 22.54032258064516, 23.064516129032256, 23.064516129032256, 23.064516129032256, 23.588709677419352, 23.588709677419352, 23.588709677419352, 24.112903225806452, 24.112903225806452, 24.112903225806452, 24.112903225806452, 24.112903225806452, 24.637096774193548, 24.637096774193548, 24.637096774193548, 24.637096774193548, 24.637096774193548, 25.161290322580644, 25.161290322580644, 25.161290322580644, 25.68548387096774, 25.68548387096774, 25.68548387096774, 25.68548387096774, 25.68548387096774, 26.209677419354836, 26.209677419354836, 26.209677419354836, 26.733870967741936, 26.733870967741936, 26.733870967741936, 27.258064516129032, 27.258064516129032, 27.258064516129032, 27.782258064516128, 28.306451612903224, 28.83064516129032, 29.354838709677416, 29.879032258064516, 30.403225806451612, 30.927419354838708, 31.451612903225804, 31.9758064516129, 32.5, 33.024193548387096, 33.548387096774192, 34.072580645161288, 34.596774193548384, 35.12096774193548, 35.645161290322577, 36.169354838709673, 36.693548387096776, 37.217741935483872, 37.741935483870968, 38.266129032258064, 38.79032258064516, 39.314516129032256, 39.838709677419352, 40.362903225806448, 40.887096774193544, 41.411290322580641, 41.935483870967744, 42.45967741935484, 42.983870967741936, 43.508064516129032, 44.032258064516128, 44.556451612903224, 45.08064516129032, 45.604838709677416, 46.129032258064512, 46.129032258064512, 46.129032258064512, 46.653225806451609, 47.177419354838705, 47.701612903225808, 48.225806451612904, 48.75, 49.274193548387096, 49.798387096774192, 50.322580645161288, 50.846774193548384, 51.37096774193548, 51.895161290322577, 52.419354838709673, 52.943548387096769, 53.467741935483872, 53.991935483870968, 54.516129032258064, 55.04032258064516, 55.564516129032256, 56.088709677419352, 56.612903225806448, 57.137096774193544, 57.661290322580641, 58.185483870967737, 58.709677419354833, 59.233870967741936, 59.758064516129032, 60.282258064516128, 60.806451612903224, 61.33064516129032, 61.854838709677416, 62.379032258064512, 62.903225806451609, 63.427419354838705, 63.951612903225801, 64.475806451612897, 65.0]
S_arrivals =  [0.0, 9.0258103928685394, 15.60936813770744, 21.494318075818921, 26.991825139833175, 32.214725623641975, 37.249637684317889, 42.149990246754015, 46.936593767902821, 51.603246753428394, 56.18580135464908, 60.664163986342245, 65.08975833108812, 69.482438745453493, 73.815809602460334, 78.090304585788516, 82.266716808578963, 82.269869844774917, 82.269984997457072, 86.33857556370269, 90.309873819258456, 94.177832461302756, 97.975346773217396, 101.7275623707406, 105.44485952517289, 109.13773077733833, 112.81600794101465, 112.81673385156955, 112.81697284165203, 116.47363752266297, 120.11793470230116, 123.75190394536691, 127.37988295846604, 131.00253488233068, 134.62467271889398, 138.245981453779, 141.86611227101537, 142.88554593392112, 142.89947287242782, 145.4855228798053, 146.29188582085283, 146.36334839882988, 149.10440323233885, 149.69798939013313, 149.85521851608786, 152.72251852347273, 153.10372551616643, 153.3658139106285, 156.33970406960395, 156.50893760558134, 156.89037509992619, 159.91348247651831, 159.95493635644829, 159.95580279723734, 159.9558547565376, 160.4250736835933, 163.31721254099514, 163.56690879727014, 163.96781680689026, 166.71999376802594, 167.17816601404792, 167.51693586425898, 170.12167433750329, 170.78831415232639, 171.07105497093758, 173.52211043896662, 174.39411403737572, 174.39708815299863, 174.3971100349637, 174.62927018566967, 176.92114014122063, 177.99561996554789, 178.19081005283579, 180.31865400300913, 181.59636540620903, 181.75503706307873, 183.71447034169572, 185.19588590125551, 185.32145821985799, 187.1084478576139, 188.79384925689672, 188.88965772814203, 190.50044271430031, 192.39002138280435, 192.39154221517938, 192.39173869568714, 192.4592857336537, 193.89031329511207, 195.98123984913204, 195.98421394695012, 195.98422665316315, 196.03003512383887, 197.27789657682661, 199.56998268788095, 199.60167211980016, 200.66304689721267, 203.15495280946124, 203.15711904244762, 203.15722260475417, 203.17392243529457, 204.04563306397884, 206.73651983970447, 206.74658644339337, 207.42550251243449, 210.31651940448353, 210.31946744683773, 210.80248935141699, 213.8944673660574, 213.89462913895468, 214.176462472358, 217.54727533172394, 220.91477055652464, 224.27880625373763, 227.63925477901631, 230.99596790396825, 234.34878096055434, 237.69756432022578, 241.04218180037529, 244.38248762705959, 247.71834235915875, 251.04961184080446, 254.37615691333718, 257.6978372628609, 261.01453107188843, 264.3260500272176, 267.63239388022947, 270.93330721581304, 274.22868639374184, 277.51843586486279, 280.8024012747311, 284.0804539725259, 287.3524899986445, 290.61837080991955, 293.87827939648503, 297.13593666196647, 300.39153459076863, 303.64435203662134, 306.89391914653049, 310.13987004172037, 313.38195025533884, 316.61994789885699, 319.8536428832989, 323.08286996769039, 326.30748528982105, 329.51990564917162, 329.52729284026043, 329.52740133202451, 332.71558480869822, 335.90135930474833, 339.07795722103737, 342.24621673342608, 345.40808217288316, 348.56270486176322, 351.70968309890259, 354.84872591813519, 357.97999098705702, 361.10308947154448, 364.21792212974913, 367.31332035851852, 370.3953490475231, 373.46706540817081, 376.53408180067811, 379.59480253147387, 382.64831712820097, 385.69407686911126, 388.73173538286647, 391.76101746088494, 394.78171890565915, 397.79467306076953, 400.80491535635957, 403.81184689010632, 406.81448237112994, 409.81237328681823, 412.80490323729708, 415.79182634076352, 418.73855025749634, 421.67958207802781, 424.61323499975219, 427.54472291391289, 430.47266957202191, 433.39613557440657, 436.31439474111039, 439.2269928498871]
ONE_DEGREE_IN_KM = 30.32335042414645


def view_section_with_envelopes(top_level_dir,file,inv_name,
  dir_type='pdart_dir',pre_filt_main=None,pre_filt_env=None,output='VEL',
  smooth_periods=10,time_before=100,time_after=2000,
  stations=STATION_LIST,channel='MHZ',start_no=0,end_no=None,exclude_list=None,scale=10,
  xlim=None,ylim=(0,1400),title=None):
    """
    View several traces at once for all the artificial impacts.

    :param top_level_dir: string
        full path the directory 
    :param file: 
        
    :param inv_name: 
    :param dir_type: 
    :param pre_filt_main: 
    :param pre_filt_env: 
    :param output: default='VEL'
    :param top_level_dir: 
    :param top_level_dir: 
    """

    file_out = file.replace('.xml', '_out.xml')

    # read the response file
    inv = read_inventory(inv_name)
    
    # print(inv)
    # for network in inv:
    #     for station in network[0]:
    #         print(station)
    # print(inv[0][0][0].response)
    
    catalog = read_events(file)
    quit = False
    
    fig, ax = plt.subplots(figsize=(7.2, 10))

    # estimates for the P and S waves? 
    P_kilometers = [i*ONE_DEGREE_IN_KM for i in P_degrees]
    S_kilometers = [i*ONE_DEGREE_IN_KM for i in S_degrees]

    plt.plot(P_arrivals, P_kilometers, color='#4B0082')
    plt.plot(S_arrivals, S_kilometers, color='#FF1493', linestyle='dashed')

    if end_no is None:
        end_no = len(catalog) + 1
    for i, ev in enumerate(catalog[start_no:end_no]):
        
        if ev.event_type != 'crash':
            print('Event is not a crash: ', start_no+i, ev)
            continue
        else: 
            print('stuff found ')
    
        # use the preferred origin if it exists
        origin = ev.preferred_origin()
        if origin is None: 
            try: 
                origin = ev.origins[0]
            except: 
                print('No origin found for: \m{}'.format(str(ev)))
                break 

        starttime = origin.time - time_before - 60
        endtime = origin.time + time_after + 60
        

        # if origin.time is not None: 
        #     starttime  = origin.time
        #     print('No origin set for Event {}'.format(str(starttime)))
        #     continue # continue event loop 

        for station_code in stations: 

            # else: 
            #     picks = ev.picks
            #     for pick in picks:
            #         print('now doing some more picks')
            #         if pick.waveform_id.station_code == station_code: 
            #             print('Event in catalog number {}, Pick Time {}, Phase hint {}'.format(i+start_no, pick.time, pick.phase_hint))
            #             # continue
            #             if not starttime:
            #                 starttime = pick.time
            #             else:
            #                 if pick.time < starttime:
            #                     starttime = pick.time


            latitude = LATITUDES[STATION_LIST.index(station_code)]
            longitude = LONGITUDES[STATION_LIST.index(station_code)]
            MH1_azimuth = MH1_AZIMUTHS[STATION_LIST.index(station_code)]
            MH2_azimuth = MH2_AZIMUTHS[STATION_LIST.index(station_code)]

            distance, azimuth_A_B, azimuth_B_A =  gps2dist_azimuth(
              origin.latitude, origin.longitude, latitude, longitude, MOON_RADIUS, MOON_FLATTENING)
            distance = distance/1000.
            print(ev.event_type)
            for desc in ev.event_descriptions:
                print(desc.text)
            print('Event: ', start_no+i, ev)
            print('{}-{}'.format(station_code, channel))
            print('Distance: {:.1f} km, Azimuth: {:.1f} deg, Back Azimuth: {:.1f} deg'.format(distance, azimuth_A_B, azimuth_B_A))

            stream = find_seismogram(top_level_dir,starttime,endtime,
              stations=[station_code],dir_type=dir_type,channels=['MHZ','MH1','MH2'])

            if stream is None or len(stream) == 0: 
                print('Event at {} not found for station {}'.format(str(starttime),station_code))
                continue # continue station loop 
         
            working_stream = process_remove_response_Apollo(stream.select(channel=channel).copy(), inv, pre_filt=pre_filt_main,output=output)
            working_stream.trim(starttime=starttime+60,endtime=endtime-60)
            # calculate normalization factor 
            # smooth this strongly to help with the normalization 
            stream_norm = process_envelope(stream.select(channel=channel).copy(), inv, pre_filt=pre_filt_main,output=output,smooth_periods=400,square=False)
            trace_max = stream_norm[0].max()
            ax.plot(working_stream[0].times()-time_before, (working_stream[0].data/trace_max)*scale+distance, color='gray',zorder=1)
            for i, pre_filt in enumerate(pre_filt_env): 
                stream_env = process_envelope(stream.select(
                  channel=channel).copy(), inv, pre_filt=pre_filt,
                  output=output,smooth_periods=smooth_periods,square=False)
                stream_env.trim(starttime=starttime+60,endtime=endtime-60)
                
                ax.plot(stream_env[0].times()-time_before, (stream_env[0].data/trace_max)*scale+distance, color=colors[i],linewidth=2,zorder=10)
   
            
            picks = ev.picks
            for pick in picks:
                if pick.phase_hint == 'P' and pick.waveform_id.station_code==station_code:
                    pick_markP = pick.time - origin.time
                    plt.vlines(x=pick_markP, ymin=distance-scale*2, ymax=distance+scale*2, color='#4B0082', linewidth=3)

                elif pick.phase_hint == 'S' and pick.waveform_id.station_code==station_code:
                    pick_markS = pick.time - origin.time
                    plt.vlines(x=pick_markS, ymin=distance-scale*2, ymax=distance+scale*2, color='#FF1493', linewidth=3)

                elif (pick.phase_hint in ('t_max','t_onset') and 
                  pick.filter_id is not None and  
                  pick.waveform_id.station_code==station_code): 

                    pick_mark = pick.time - origin.time
                    if pick.filter_id == filter_ids[0]:
                        color=colors[0]
                    elif pick.filter_id == filter_ids[1]:
                        color=colors[1]
                    elif pick.filter_id == filter_ids[2]:
                        color=colors[2]

                    plt.vlines(x=pick_mark, ymin=distance-scale*2, ymax=distance+scale*2, color=color, linewidth=3)

                elif (pick.phase_hint in ('t_max_lower','t_max_lower') and 
                  pick.filter_id is not None and  
                  pick.waveform_id.station_code==station_code): 

                    pick_mark = pick.time - origin.time
                    if pick.filter_id == filter_ids[0]:
                        color=colors[0]
                    elif pick.filter_id == filter_ids[1]:
                        color=colors[1]
                    elif pick.filter_id == filter_ids[2]:
                        color=colors[2]

                    plt.vlines(x=pick_mark, ymin=distance-scale*2, 
                      ymax=distance+scale*2, color=color, linewidth=1,
                      linestyle='dashed')


            # stream_env = stream.select(channel=channel).copy()


            
            
                        # trace_max = stream_env[0].max()
                        # stream_env = stream_orig.copy()
                # plot_seismogram_with_envelopes(working_stream, inv, output=output, 
                #      pre_filt_main=pre_filt_main, pre_filt_env=pre_filt_env,
                #      smooth_periods=smooth_periods,
                #      plot_type=plot_type,title=str(starttime))

            # manual_scale = [1, 9, 2, 2, 1.5, 5, 1]



    plt.xlabel('Time [s]',labelpad=-5)
    plt.ylabel('Offset [km]')
    if xlim is None: 
        xlim = (-time_before,time_after)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if title is not None: 
        plt.title(title)

    legend_lines = []
    legend_titles = []
    for i, pre_filt in enumerate(pre_filt_env):
        legend_lines.append(Line2D([0], [0], color=colors[i], lw=1))
        legend_titles.append('{} - {} Hz'.format(pre_filt[1],pre_filt[2]))
    # cut
    # custom_lines = [Line2D([0], [0], color=cmap(0.), lw=4),
    #                 Line2D([0], [0], color=cmap(.5), lw=4),
    #                 Line2D([0], [0], color=cmap(1.), lw=4)]

    ax.legend(legend_lines, legend_titles, loc='upper right', framealpha=0)
    # plt.savefig(filename,bbox_inches='tight')
    ax.axvline(x=0, color='gray', linewidth=1, linestyle='dashed')


    plt.show()
            # for i, tr in enumerate(stream): 
            #     x_time = tr.times()
            #     if pre_filt_env is not None:
            #         for ii, pre_filt in enumerate(pre_filt_env):
            #             stream_env = stream_orig.copy()
            #             # calculate normalization factor 
            #             # smooth this strongly to help with the normalization 
            #             process_envelope(stream_env, inv, pre_filt=pre_filt,output=output,smooth_periods=400,square=True)
            #             trace_max = stream_env[0].max()
            #             stream_env = stream_orig.copy()
            #             # process again 
            #             process_envelope(stream_env, inv, pre_filt=pre_filt,output=output,smooth_periods=smooth_periods,square=True)
            #             # don't plot them on top of each other, but shift by 1 each time
            #             # plt.plot(x_time,stream_env[0].data/trace_max+ ii,label=str(pre_filt))
            #             plt.plot(x_time,stream_env[0].data/trace_max,label='{} - {} Hz'.format(pre_filt[1],pre_filt[2]))

def onclick(event):
    global xdata
    xdata = round(event.xdata,1)
    print('Clicked {} s'.format(xdata))
