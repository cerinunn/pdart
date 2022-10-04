#!/usr/bin/env python

from __future__ import print_function
from scipy import interpolate, signal
from time import *
import numpy as np
from obspy import *
import matplotlib
import matplotlib.pyplot as plt
import os
import glob
from datetime import datetime, timedelta
from obspy.core.utcdatetime import UTCDateTime
import urllib.request
from obspy.core.util.misc import (limit_numpy_fft_cache)
from future.utils import native_str

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

# test the different types of filtering
def filter_test(
    inv_name,
    start_time,
    station,
    channel,
    pre_filt,
    end_time=None,
    outfile=None,
    output='DISP',
    water_level=None):

    # read the response file
    inv = read_inventory(inv_name)

    if end_time is None:
        time_interval = timedelta(hours=3)
        end_time = start_time + time_interval


    filename = '%s.%s.*.%s.%s.%03d.gz' % ('XA',station, channel,
       str(start_time.year), start_time.julday)
    filename = os.path.join('/Users/cnunn/lunar_data/PDART',str(start_time.year),'XA',station,channel,filename)

    stream = read(filename)
    stream = stream.select(channel=channel)
    stream.trim(starttime=start_time, endtime=end_time)

    # remove location (ground station)
    for tr in stream:
        tr.stats.location = ''

    stream = stream.merge()

    # print(pre_filt,1/pre_filt[1],1/pre_filt[2])

    # detrend
    stream = stream.split()
    stream.detrend('linear')

    # filter in time domain
    stream_time_filter = stream.copy()
    stream_time_filter.filter("bandpass",freqmin=pre_filt[1],freqmax=pre_filt[2],zerophase=True)

    for tr in stream_time_filter:
        tr.stats.location = 'time_filter'

    stream_freq_filter = stream.copy()
    # smart calculation of nfft dodging large primes
    from obspy.signal.util import _npts2nfft

    for tr in stream_freq_filter:
        nfft = _npts2nfft(len(tr.data))
        npts = len(tr.data)
        response = tr._get_response(inv)
        freq_response, freqs = \
            response.get_evalresp_response(tr.stats.delta, nfft,
                                           output='DISP')
        # print(freq_response)
        # print(freqs)
        # exit()
        # Transform data to Frequency domain
        tr.data = np.fft.rfft(tr.data, n=nfft)
        from obspy.signal.invsim import (cosine_taper, cosine_sac_taper,
                                         invert_spectrum)
        freq_domain_taper = cosine_sac_taper(freqs, flimit=pre_filt)
        tr.data *= freq_domain_taper
        # transform data back into the time domain
        tr.data = np.fft.irfft(tr.data)[0:npts]
        tr.stats.location = 'freq_filter'

    stream = stream + stream_time_filter + stream_freq_filter
    # stream = stream.trim(UTCDateTime('1976-01-13T07:16'),UTCDateTime('1976-01-13T07:20'))
    stream.plot(equal_scale=False,size=(1200,600))


def remove_response_from_seismogram(
    inv_name,
    start_time,
    station,
    channel,
    pre_filt,
    end_time=None,
    outfile=None,
    output='DISP',
    water_level=None):

    # read the response file
    inv = read_inventory(inv_name)

    if end_time is None:
        time_interval = timedelta(hours=3)
        end_time = start_time + time_interval
# xa.s12..att.1969.324.0.mseed

    filename = '%s.%s.*.%s.%s.%03d.0.mseed' % ('xa',station.lower(), channel.lower(),
       str(start_time.year), start_time.julday)
    filename = os.path.join('/Users/cnunn/lunar_data/PDART_CONTINUOUS_MAIN_TAPES',station.lower(),str(start_time.year),str(start_time.julday),filename)

    stream = read(filename)
    stream = stream.select(channel=channel)
    stream.trim(starttime=start_time, endtime=end_time)

    # remove location (ground station)
    for tr in stream:
        tr.stats.location = ''

    # detrend
    stream.detrend('linear')

    # taper the edges
    # if there are gaps in the seismogram - EVERY short trace will be tapered
    # this is required to remove the response later
    # stream.taper(max_percentage=0.05, type='cosine')

    # experiment with tapering? not tapering preserves the overall shape better
    # but it may required

    # merge the streams
    stream.merge()
    if stream.count() > 1:
        print('Too many streams - exiting')

    # find the gaps in the trace
    if isinstance(stream[0].data,np.ma.MaskedArray):
        mask = np.ma.getmask(stream[0].data)
    else:
        mask = None

    # split the stream, then refill it with zeros on the gaps
    stream = stream.split()
    stream = stream.merge(fill_value=0)

    # for i, n in enumerate(stream[0].times()):
    #     # print(n)
    #     stream[0].data[i]=np.sin(2*np.pi*(1/25)*n)

    stream.attach_response(inv)
    # print('here')

    # zero_mean=False - because the trace can be asymmetric - remove the mean ourselves
    # do not taper here - it doesn't work well with the masked arrays - often required
    # when there are gaps - if necessary taper first
    # water level - this probably doesn't have much impact - because we are pre filtering
    # stream.remove_response(pre_filt=pre_filt,output="DISP",water_level=30,zero_mean=False,taper=False,plot=True,fig=outfile)
    for tr in stream:
        remove_response(tr, pre_filt=pre_filt,output=output,water_level=water_level,zero_mean=False,taper=False,plot=True,fig=outfile)

    for tr in stream:
        tr.stats.location = 'changed'

    if mask is not None:
        stream[0].data = np.ma.array(stream[0].data, mask = mask)

    return stream


# test the different types of filtering with a sinusoid
def filter_test_sinusoid(
    inv_name,
    start_time,
    station,
    channel,
    pre_filt,
    end_time=None,
    outfile=None,
    output='DISP',
    water_level=None):

    # read the response file
    inv = read_inventory(inv_name)

    if end_time is None:
        time_interval = timedelta(hours=3)
        end_time = start_time + time_interval


    filename = '%s.%s.*.%s.%s.%03d.gz' % ('XA',station, channel,
       str(start_time.year), start_time.julday)
    filename = os.path.join('/Users/cnunn/lunar_data/PDART',str(start_time.year),'XA',station,channel,filename)

    stream = read(filename)
    stream = stream.select(channel=channel)
    stream.trim(starttime=start_time, endtime=end_time)

    # remove location (ground station)
    for tr in stream:
        tr.stats.location = ''

    stream = stream.merge()
    for i, n in enumerate(stream[0].times()):
        stream[0].data[i]=(1000*np.sin(2*np.pi*(1/25)*n))

    # print(pre_filt,1/pre_filt[1],1/pre_filt[2])

    # detrend
    stream = stream.split()
    stream.detrend('linear')

    # filter in time domain
    stream_time_filter = stream.copy()
    stream_time_filter.filter("bandpass",freqmin=pre_filt[1],freqmax=pre_filt[2],zerophase=True)

    for tr in stream_time_filter:
        tr.stats.location = 'time_filter'

    stream_freq_filter = stream.copy()
    # smart calculation of nfft dodging large primes
    from obspy.signal.util import _npts2nfft

    for tr in stream_freq_filter:
        nfft = _npts2nfft(len(tr.data))
        npts = len(tr.data)
        response = tr._get_response(inv)
        freq_response, freqs = \
            response.get_evalresp_response(tr.stats.delta, nfft,
                                           output='DISP')
        # print(freq_response)
        # print(freqs)
        # exit()
        # Transform data to Frequency domain
        tr.data = np.fft.rfft(tr.data, n=nfft)
        from obspy.signal.invsim import (cosine_taper, cosine_sac_taper,
                                         invert_spectrum)
        freq_domain_taper = cosine_sac_taper(freqs, flimit=pre_filt)
        tr.data *= freq_domain_taper
        # transform data back into the time domain
        tr.data = np.fft.irfft(tr.data)[0:npts]
        tr.stats.location = 'freq_filter'

    stream = stream + stream_time_filter + stream_freq_filter
    # stream = stream.trim(UTCDateTime('1976-01-13T07:16'),UTCDateTime('1976-01-13T07:20'))
    stream.plot(equal_scale=False,size=(1200,600))

    # original_stream = stream.copy()
    #
    # # exit()
    #
    # stream.split()
    #
    # f1 = stream.copy()
    # # stream.filter("bandpass",freqmin=0.05,freqmax=0.2)
    # f1.filter("bandpass",freqmin=1/15,freqmax=1/5)
    # for tr in f1:
    #     tr.stats.location = 'f1'
    # f2 = stream.copy()
    # f2.filter("bandpass",freqmin=1/20,freqmax=1/10)
    # for tr in f2:
    #     tr.stats.location = 'f2'
    # f3 = stream.copy()
    # f3.filter("bandpass",freqmin=1/25,freqmax=1/15)
    # for tr in f3:
    #     tr.stats.location = 'f3'
    # f4 = stream.copy()
    # f4.filter("bandpass",freqmin=1/35,freqmax=1/20)
    # # for tr in f4:
    # #     tr.stats.location = 'f4'
    #
    # f5 = stream.copy()
    # f5.filter("bandpass",freqmin=1/40,freqmax=1/25)
    # for tr in f5:
    #     tr.stats.location = 'f5'
    #
    # f6 = stream.copy()
    # f6.filter("bandpass",freqmin=1/45,freqmax=1/30)
    # for tr in f6:
    #     tr.stats.location = 'f6'
    # # stream.merge()
    #
    # f7 = stream.copy()
    # f7.filter("bandpass",freqmin=1/50,freqmax=1/35)
    # for tr in f7:
    #     tr.stats.location = 'f7'
    # # stream.merge()
    #
    # f8 = stream.copy()
    # f8.filter("bandpass",freqmin=1/55,freqmax=1/40)
    # for tr in f8:
    #     tr.stats.location = 'f8'
    # # stream.merge()
    #
    # f9 = stream.copy()
    # f9.filter("bandpass",freqmin=1/55,freqmax=1/40)
    # for tr in f9:
    #     tr.stats.location = 'f9'
    # # stream.merge()
    #
    # f10 = stream.copy()
    # f10.filter("bandpass",freqmin=1/55,freqmax=1/40)
    # for tr in f10:
    #     tr.stats.location = 'f10'
    # # stream.merge()
    #
    # stream = f4
    #

    # stream = stream + f1 + f2 + f3 + f4 + f5 + f6 + f7 +f8 +f9 +f10
    # stream = stream.trim(UTCDateTime('1976-01-13T07:16'),UTCDateTime('1976-01-13T07:20'))
    # stream.plot(equal_scale=False,size=(1200,600))






    # print('temp')
    # stream.cutout(start_time+600,start_time+667)
    # stream[1].stats.starttime = stream[0].stats.endtime + 0.15094
    # stream = stream.merge()

    # stream.plot()
    # exit()
    # stream.filter("lowpass", freq=1/20)
    # stream.plot()
    # exit()


def remove_response_from_seismogram_method2(
    inv_name,
    start_time,
    station,
    channel,
    pre_filt,
    end_time=None,
    outfile=None,
    output='DISP',
    water_level=None):

    # read the response file
    inv = read_inventory(inv_name)

    if end_time is None:
        time_interval = timedelta(hours=3)
        end_time = start_time + time_interval


    filename = '%s.%s.*.%s.%s.%03d.gz' % ('XA',station, channel,
       str(start_time.year), start_time.julday)
    filename = os.path.join('/Users/cnunn/lunar_data/PDART',str(start_time.year),'XA',station,channel,filename)

    stream = read(filename)
    stream = stream.select(channel=channel)
    stream.trim(starttime=start_time, endtime=end_time)

    # remove location (ground station)
    for tr in stream:
        tr.stats.location = ''

    # detrend
    stream.detrend('linear')

    # taper the edges
    # if there are gaps in the seismogram - EVERY short trace will be tapered
    # this is required to remove the response later
    # stream.taper(max_percentage=0.05, type='cosine')

    # experiment with tapering? not tapering preserves the overall shape better
    # but it may required

    # merge the streams
    stream.merge()
    if stream.count() > 1:
        print('Too many streams - exiting')

    # find the gaps in the trace
    if isinstance(stream[0].data,np.ma.MaskedArray):
        mask = np.ma.getmask(stream[0].data)
    else:
        mask = None

    # split the stream, then refill it with zeros on the gaps
    stream = stream.split()
    # stream.merge(fill_value=0)


    stream.attach_response(inv)

    # zero_mean=False - because the trace can be asymmetric - remove the mean ourselves
    # do not taper here - it doesn't work well with the masked arrays - often required
    # when there are gaps - if necessary taper first
    # water level - this probably doesn't have much impact - because we are pre filtering
    # stream.remove_response(pre_filt=pre_filt,output="DISP",water_level=30,zero_mean=False,taper=False,plot=True,fig=outfile)
    stream.remove_response(pre_filt=pre_filt,output=output,water_level=water_level,zero_mean=False,taper=False,plot=True)

    stream = stream.merge()


    stream = stream.trim(UTCDateTime('1971-02-07T00:43'),UTCDateTime('1971-02-07T00:58'))


    stream.plot()
    # if mask is not None:
    #     stream[0].data = np.ma.array(stream[0].data, mask = mask)






def some_example_times():
    pass
    # station='S12'
    # time = UTCDateTime('1973-11-10T') #peak
    # time = UTCDateTime('1974-10-17T') #flat
    # time = UTCDateTime('1975-05-01T') #peak
    # time = UTCDateTime('1977-03-27T') #flat
    # time = UTCDateTime('1977-09-15T') #peak
    #
    # station='S14'
    # time = UTCDateTime('1973-11-10T') #peak OK
    # time = UTCDateTime('1976-11-10T') #flat
    # time = UTCDateTime('1977-09-15T') #peak
    #
    # station='S15'
    # time = UTCDateTime('1973-11-10T') #peak
    # time = UTCDateTime('1975-06-29T') #flat
    # time = UTCDateTime('1977-09-15T') #peak
    #
    # station='S16'
    # time = UTCDateTime('1973-11-10T') #peak
    # time = UTCDateTime('1975-06-29T') #flat
    # time = UTCDateTime('1977-09-15T') #peak

def example_without_noise():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [0.2, 0.3,0.9,1.3]
    end_time = UTCDateTime('1971:02:07T01:28')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      end_time=end_time)

    # saved as response_1971-02-07T00T00-35.25.pdf


def large_pulse_flat_mode():
    # this is a large calibration in flat mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime(year=1976,julday=149,hour=16,minute=8)
    station = 'S15'
    channel = 'MH1'
    pre_filt = [0.2, 0.3,0.9,1.3]
    pre_filt = None
    # pre_filt = [0.001, 0.005,0.02,0.05]
    pre_filt = [0.001, 0.005,2,3]
    water_level = 30
    # water_level = None
    end_time = UTCDateTime(year=1976,julday=149,hour=16,minute=14)
    end_time = UTCDateTime(year=1976,julday=149,hour=16,minute=8) + 302.4
    output='ACC'

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time+60,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      end_time=end_time,
      output=output,
      water_level=water_level)



def example_flat_mode():
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1976:01:13T06:00')
    station = 'S12'
    channel = 'MHZ'
    pre_filt = [0.2, 0.3,0.9,1.3]
    end_time = UTCDateTime('1976:01:13T10')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      end_time=end_time)

    # saved as response_1976:01:13T06:00_flat1.pdf

def example_flat_mode_lower_freq():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1976:01:13T06:00')
    station = 'S12'
    channel = 'MHZ'
    pre_filt = [0.07, 0.1,0.3,0.4]
    end_time = UTCDateTime('1976:01:13T10')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      end_time=end_time)

    # saved as response_1976:01:13T06:00_flat2.pdf



    # saved as response_1976:01:13T06:00_flat_10s_period_zoomed.png

# this isn't giving the correct answer. Something wrong
def simulate():
    from obspy.core import Trace
    trace = Trace()
    trace.stats.delta = 0.15094
    trace.data = np.zeros(3000)
    trace.data[1500:] = 1
    trace.stats.network = 'XA'
    trace.stats.station = 'S15'
    trace.stats.channel = 'MHZ'
    print(trace.stats)
    # exit()
    date = UTCDateTime('1972-12-01')


    trace.integrate()
    trace.integrate()
    # trace.plot()
    #
    #
    # # trace.plot()
    #
    # exit()

    #     # define a filter band to prevent amplifying noise during the deconvolution
    #     pre_filt = (0.005, 0.006, 30.0, 35.0)
    #
    #     # this can be the date of your raw data or any date for which the
    #     # SEED RESP-file is valid
    #     date = t1
    #
    seedresp = {'filename': '/Users/cnunn/lunar_data/IRIS_dataless_seed/RESP.XA.S15..MHZ',
                # RESP filename
                # when using Trace/Stream.simulate() the "date" parameter can
                # also be omitted, and the starttime of the trace is then used.
                # Units to return response in ('DIS', 'VEL' or ACC)
               'date': date,
                'units': 'DIS'
                }

    paz_sts2 = {
        'poles':
        [-3.560470E-01  + 2.837020E+00j,
        -3.560470E-01  - 2.837020E+00j,
        -6.283190E-02  + 0.000000E+00j,
        -8.062370E+00  + 3.339540E+00j,
        -8.062370E+00  - 3.339540E+00j,
        -8.062370E+00  + 3.339540E+00j,
        -8.062370E+00  - 3.339540E+00j,
        -3.339540E+00  + 8.062370E+00j,
        -3.339540E+00  - 8.062370E+00j,
        -3.339540E+00  + 8.062370E+00j,
        -3.339540E+00  - 8.062370E+00j],
        'zeros':[0.000000E+00  + 0.000000E+00j,
        0.000000E+00  + 0.000000E+00j,
        0.000000E+00  + 0.000000E+00j],
        'gain': 60077000.0,
        'sensitivity': 2516778400.0}

# gain and senstivity not set
#		Complex poles:
#		  i  real          imag          real_error    imag_error



    #
    #     # Remove instrument response using the information from the given RESP file
    trace.simulate(paz_simulate=paz_sts2)
    trace.plot()
    exit()
    #
    # # plot original and simulated data
    # tr = st[0]
    # tr_orig = st_orig[0]
    # time = tr.times()
    #
    # plt.subplot(211)
    # plt.plot(time, tr_orig.data, 'k')
    # plt.ylabel('STS-2 [counts]')
    # plt.subplot(212)
    # plt.plot(time, tr.data, 'k')
    # plt.ylabel('Displacement [m]')
    # plt.xlabel('Time [s]')
    # plt.show()

def calibration_pulse_remove():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    time = UTCDateTime('1972-12-31T07:45:05.407000Z')
    start_time = time - 10 * 60
    end_time = time + 10 * 60
    # start_time = UTCDateTime('1976:01:13T06:00')
    station = 'S12'
    channel = 'MHZ'
    # pre_filt = [0.2, 0.3,0.9,1.3]
    pre_filt = None
    # pre_filt = [0.1, 0.2,0.9,1.3]
    # pre_filt = [0.0005, 0.001,0.1,0.4]
    # pre_filt = [0.005, 0.01,0.02,0.04]
    #
    # end_time = UTCDateTime('1976:01:13T10')
    #  stream.filter("lowpass", freq=1/20)
    # pre_filt = [0.001, 0.002, 1/20, 1/18]

    # strongly filtered
    pre_filt = [0.001, 0.002, 1/40, 1/20]
    water_level = None

    # not filtered, but with a water level to see the pulse
    # pre_filt = None
    # water_level=30
    output='ACC'

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      end_time=end_time,
      water_level=water_level,
      output=output)




def example_with_noise_my_filtering():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [0.1,0.2,0.9,1.3]
    end_time = UTCDateTime('1971:02:07T02:35.25')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=None,
      water_level=60,
      end_time=end_time)

def example_with_noise():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [0.1,0.2,0.9,1.3]
    end_time = UTCDateTime('1971:02:07T02:35.25')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      end_time=end_time)

    # saved as response_1971-02-07T00T00-35.25_withnoise_zoomed.pdf

def example_with_noise_method2():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [0.2,0.3,0.9,1.3]
    end_time = UTCDateTime('1971:02:07T02:35.25')

    remove_response_from_seismogram_method2(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      end_time=end_time)

# this tests the filtering with a sinusoid
def test_sinusoid():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [1/45, 1/35,1/20,1/10]
    end_time = UTCDateTime('1971:02:07T02:35.25')

    filter_test_sinusoid(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      water_level=60,
      end_time=end_time)

# this tests the filtering
def run_test_filter():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [1/45, 1/35,1/20,1/10]
    end_time = UTCDateTime('1971:02:07T02:35.25')

    filter_test(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      water_level=60,
      end_time=end_time)

def example1_peaked_mode_4point5_Hz():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [0.1, 0.3,0.7,1]
    end_time = UTCDateTime('1971:02:07T02:35.25')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      water_level=None,
      end_time=end_time)
    # saved as response_example1_peaked_mode_4point5_Hz.png

def example2_peaked_mode_10s_period():
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [0.03, 0.07,0.13,0.18]
    end_time = UTCDateTime('1971:02:07T02:35.25')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      water_level=None,
      end_time=end_time)
    # saved as response_example2_peaked_mode_10s_period_zoom.png

def example3_flat_mode_20s_period():
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1976:01:13T06:00')
    station = 'S12'
    channel = 'MHZ'
    # pre_filt = [0.02,0.03,0.07,0.1]
    pre_filt = [0.03,0.04,0.07,0.1]
    end_time = UTCDateTime('1976:01:13T10')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      end_time=end_time,
      water_level=None)
    # saved as response_example3_flat_mode_20s_period.png



def example4_peaked_mode_1971_02_07_4point5_Hz():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S14'
    channel = 'MHZ'
    pre_filt = [0.1, 0.3,0.7,1]
    end_time = UTCDateTime('1971:02:07T02:35.25')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      water_level=None,
      end_time=end_time)
    # saved as response_example4_peaked_mode_1971_02_07_4point5_Hz.png

def Apollo12_LM_impact_spectra():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1969-11-20TT22:17:17.700000Z')
    station = 'S12'
    channel = 'MHZ'
    pre_filt = [0.1, 0.3,0.7,1]
    # end_time = UTCDateTime('1971:02:07T02:35.25')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      water_level=None,
      end_time=end_time)
    # saved as response_example4_peaked_mode_1971_02_07_4point5_Hz.png

def example_with_noise_test():
    # peaked mode
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    start_time = UTCDateTime('1971-02-07T00:00:35.250000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [1/45, 1/35,1/20,1/10]
    end_time = UTCDateTime('1971:02:07T02:35.25')

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      water_level=60,
      end_time=end_time)

def check_calibration():
    # just to show whether flat or peaked mode

    # onset = UTCDateTime('1970-03-26T02:04:48.295000Z')
    # onset2 = UTCDateTime(')
    inv_name = "/Users/cnunn/lunar_data/IRIS_dataless_seed/XA.1969-1977.xml"
    # example 1
    start_time =  UTCDateTime('1970-03-26T02:04:48.295000Z')
    # example 2
    start_time =  UTCDateTime('1976-05-18T04:29:04.511000Z')
    station = 'S12'
    channel = 'MH2'
    pre_filt = [1/45, 1/35,1/20,1/10]
    end_time = start_time + 60

    remove_response_from_seismogram(inv_name=inv_name,
      start_time=start_time,
      station=station,
      channel=channel,
      pre_filt=pre_filt,
      water_level=60,
      end_time=end_time)

# 0.075, 0.083, 0.124, 0.15


if __name__ == "__main__":

    Apollo12_LM_impact_spectra()

    # check_calibration()

    # waveform instrument response
    # example1_peaked_mode_4point5_Hz()
    # example2_peaked_mode_10s_period()
    # example3_flat_mode_20s_period()
    # example4_peaked_mode_1971_02_07_4point5_Hz()
    # this is a nice one

    # example
    # test_sinusoid()
    # run_test_filter()

    # example_with_noise_test()
    # example_with_noise_my_filtering()

    # example_with_noise_method2()

    # large_pulse_flat_mode()
    # calibration_pulse_remove()
    # simulate()


    # example_without_noise()
    # example_with_noise()
    # example_flat_mode()
    # example_flat_mode_lower_freq()



    # some examples
    # station='S12'
    # time = UTCDateTime('1973-11-10T') #peak
    # time = UTCDateTime('1974-10-17T') #flat
    # time = UTCDateTime('1975-05-01T') #peak
    # time = UTCDateTime('1977-03-27T') #flat
    # time = UTCDateTime('1977-09-15T') #peak
    #
    # station='S14'
    # time = UTCDateTime('1973-11-10T') #peak OK
    # time = UTCDateTime('1976-11-10T') #flat
    # time = UTCDateTime('1977-09-15T') #peak
    #
    # station='S15'
    # time = UTCDateTime('1973-11-10T') #peak
    # time = UTCDateTime('1975-06-29T') #flat
    # time = UTCDateTime('1977-09-15T') #peak
    #
    # station='S16'
    # time = UTCDateTime('1973-11-10T') #peak
    # time = UTCDateTime('1975-06-29T') #flat
    # time = UTCDateTime('1977-09-15T') #peak
