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
from obspy.core.util.misc import limit_numpy_fft_cache
from obspy.signal.filter import envelope
from future.utils import native_str
import os.path
import numpy as np
from matplotlib.dates import date2num
import math


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})

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

        # print('Before ', stream)
        stream = stream.trim(starttime=starttime,endtime=endtime)
        # print('After ', stream)

        if stream is not None and len(stream) > 0: 
            for tr in stream: 
                tr.stats.location = ''
                if tr.stats.channel not in channels:
                    stream.remove(tr)

            stream.merge()

    return stream

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


    stream.attach_response(inv)

    # remove the response when we need to remove the mean and deal with 
    # gaps in the seismogram
    stream = stream.merge()

    for tr in stream:
        # remove the mean
        mean = tr.data.mean()
        tr.data = tr.data - mean

        # find the gaps in the trace
        if isinstance(tr.data,np.ma.MaskedArray):
            mask = np.ma.getmask(tr.data)
            # refill the gaps in the trace with zeros 
            tr.data = np.ma.MaskedArray(tr.data,fill_value=0)
        else:
            mask = None

        # zero_mean=False - because the trace can be asymmetric - remove the mean ourselves
        # do not taper here - it doesn't work well with the masked arrays - often required
        # when there are gaps - if necessary taper first
        # water level - this probably doesn't have much impact - because we are pre filtering
        # stream.remove_response(pre_filt=pre_filt,output="DISP",water_level=30,zero_mean=False,taper=False,plot=True,fig=outfile)
        # remove_response(tr, pre_filt=pre_filt,output=output,
        #         water_level=None,zero_mean=False,taper=False,plot=False)

        trace, freqs, spectra =remove_response(tr, pre_filt=pre_filt,output=output,
                water_level=None,zero_mean=False,taper=False,plot=False,return_spectra=True)

        # plt.loglog(freqs, spectra, alpha=0.5)
        # plt.show()
        # tr.spectrogram(log=False,wlen=100,dbscale=False,
        # clip=[0.85,1.0])
        #   title='{} S12 '.format(str(starttime)))
        # exit()

        # now put the masked areas back 
        if mask is not None:
            tr.data = np.ma.array(tr.data, mask = mask)  



def process_envelope(stream,inv,pre_filt,output='VEL',smooth_periods=10,square=True):

    stream.attach_response(inv)
    
    # remove the response when we need to remove the mean and deal with 
    # gaps in the seismogram
    stream = stream.merge()

    for tr in stream:
        # remove the mean
        mean = tr.data.mean()
        tr.data = tr.data - mean

        # find the gaps in the trace
        if isinstance(tr.data,np.ma.MaskedArray):
            mask = np.ma.getmask(tr.data)
            # refill the gaps in the trace with zeros 
            tr.data = np.ma.MaskedArray(tr.data,fill_value=0)
        else:
            mask = None

        # zero_mean=False - because the trace can be asymmetric - remove the mean ourselves
        # do not taper here - it doesn't work well with the masked arrays - often required
        # when there are gaps - if necessary taper first
        # water level - this probably doesn't have much impact - because we are pre filtering
        # stream.remove_response(pre_filt=pre_filt,output="DISP",water_level=30,zero_mean=False,taper=False,plot=True,fig=outfile)
        remove_response(tr, pre_filt=pre_filt,output=output,
                water_level=None,zero_mean=False,taper=False,plot=False)

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

        # now put the masked areas back 
        if mask is not None:
            tr.data = np.ma.array(tr.data, mask = mask) 



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

def plot_seismogram_with_envelopes(stream, inv, output='VEL', pre_filt_main=None, 
    pre_filt_env=None,smooth_periods=10,plot_type='with_seismogram',title=None):



    # stream = stream.merge()
    stream = stream.select(channel='MHZ')

    if title:
        title = '{} - {}'.format(title, 'MHZ')

    stream_orig = stream.copy()
    process_remove_response_Apollo(stream, inv, pre_filt=pre_filt_main, output=output)

    if plot_type=='with_seismogram':
        fig = plt.figure(figsize=(8,3))
        for i, tr in enumerate(stream): 
            x_time = tr.times()
            plt.plot(x_time,tr.data,color='k')
            if pre_filt_env is not None:
                for pre_filt in pre_filt_env:
                    stream_env = stream_orig.copy()
                    process_envelope(stream_env, inv, pre_filt=pre_filt,smooth_periods=smooth_periods,square=False)
                    # plot the square 
                    plt.plot(x_time,stream_env[0].data,label='{} - {} Hz'.format(pre_filt[1],pre_filt[2]))

        if title:
            plt.title(title)
        plt.xlabel('Seconds',fontsize='12')
        plt.ylabel(r'Acceleration [m/s$^{-2}$]',fontsize='12')
        plt.legend()
        plt.subplots_adjust(bottom=0.15)
        plt.show()
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
                    process_envelope(stream_env, inv, pre_filt=pre_filt,smooth_periods=400,square=True)
                    trace_max = stream_env[0].max()
                    stream_env = stream_orig.copy()
                    # process again 
                    process_envelope(stream_env, inv, pre_filt=pre_filt,smooth_periods=smooth_periods,square=True)
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
        plt.show()

        

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


def view_catalog_with_envelopes(top_level_dir,file,inv_name,
  dir_type='pdart_dir',pre_filt_main=None,pre_filt_env=None,output='VEL',
  smooth_periods=10,time_before=600,time_after=3600,plot_type='with_seismogram'): 

    file_out = file.replace('.xml', '_out.xml')

    # read the response file
    inv = read_inventory(inv_name)
    
    # print(inv)
    # for network in inv:
    #     for station in network[0]:
    #         print(station)
    # print(inv[0][0][0].response)
    # exit()
    station_code = 'S12'
    
    catalog = read_events(file)
    quit = False
    
    # put the whole thing in a try, because we need to save if something goes 
    # wrong
    # try: 
    # start_no = 129
    start_no = 0
    for i, ev in enumerate(catalog[start_no:]):
            starttime = None
            endtime = None
            picks = ev.picks
            for pick in picks:
                print(i+start_no, pick.time)
                # continue
                if not starttime:
                    starttime = pick.time
                else:
                    if pick.time < starttime:
                        starttime = pick.time
    
            if starttime:
            
                # print('Temporarily changed to 1/2 hour')
                # starttime -= 1800.
                starttime -= time_before
                endtime = starttime + time_after
            
                while True:
                    stream = find_seismogram(top_level_dir,starttime,endtime,
                      stations=[station_code],dir_type=dir_type)
            
                    plot_seismogram_with_envelopes(stream, inv, output=output, 
                         pre_filt_main=pre_filt_main, pre_filt_env=pre_filt_env,
                         smooth_periods=smooth_periods,
                         plot_type=plot_type,title=str(starttime))
            
                    input_string = '''view again (v), next (n), quit(q)\n'''
                    output_str = input(input_string)
                    if output_str == 'v': #view again 
                        pass 
                    elif output_str == 'q': #quit
                        quit=True
                        break
                    elif output_str == 'n': #next
                        quit=False
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
                if quit:
                    break

    # finally: 
    # 
    #     catalog.write(filename=file_out, format='QUAKEML') 
    #     print(i+start_no)

# category = other
# unit = m/(s*s)
# genericAmplitude = number
        

if __name__ == "__main__":
    # top_level_dir = '/Users/nunn/lunar_data/PDART'
    top_level_processed_dir = '/Users/nunn/lunar_data/PDART_PROCESSED'
    inv_name = "/Users/nunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    # these catalog are just some shorter versions - e.g. for one year
    # the full catalog takes a long time to run
    # view_catalog(top_level_processed_dir,'catalogs/1973_S12_deep.xml',inv_name,dir_type='processed_dir',output='VEL')          
    # view_catalog(top_level_processed_dir,'catalogs/1973_S12_shallow.xml',inv_name,dir_type='processed_dir',output='VEL') 
    # view_catalog(top_level_processed_dir,'catalogs/1973_S12_met.xml',inv_name,dir_type='processed_dir',output='VEL')
    # view_catalog(top_level_processed_dir,'catalogs/S12_shallow.xml',inv_name,dir_type='processed_dir',output='VEL') 
    pre_filt_main = [0.2,0.3,0.9,1.3]
    pre_filt_env = [[0.1,0.25,0.75,1],[0.5,0.75,1.25,1.5],[1,1.25,1.75,2.25]]

    # pre_filt_main = [0.2,0.3,0.9,1.3]
    # pre_filt_env = [[0.05,0.08,0.15,0.25],[0.08,0.15,0.35,0.6],[0.3,0.35,0.55,0.7],[0.35,0.55,0.75,0.9]]



    # view_catalog_with_envelopes(top_level_processed_dir,'../lunar_proposal/catalogs/S12_impact.xml',
    #   inv_name,dir_type='processed_dir',pre_filt_main=pre_filt_main,pre_filt_env=pre_filt_env,output='ACC',smooth_periods=20)

    # view_catalog_with_envelopes(top_level_processed_dir,
    #   '../lunar_proposal/catalogs/S12_shallow.xml',inv_name,
    #  dir_type='processed_dir',pre_filt_main=pre_filt_main,
    #  pre_filt_env=pre_filt_env,output='ACC',smooth_periods=10,
    #  plot_type='with_seismogram')

    view_catalog_with_envelopes.view_catalog_with_envelopes(top_level_processed_dir,
      '../lunar_proposal/catalogs/S12_shallow.xml',inv_name,
     dir_type='processed_dir',pre_filt_main=pre_filt_main,
     pre_filt_env=pre_filt_env,output='ACC',smooth_periods=10,
     plot_type='normalized')



    # plot_example()
