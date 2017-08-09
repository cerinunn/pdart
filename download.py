#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Import Apollo data to MINISEED

:copyright:
    The ObsPy Development Team (devs@obspy.org) & C. J. Ammon & J. MacCarthy
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from datetime import timedelta
import os
import io
import gzip

from obspy.core.utcdatetime import UTCDateTime
from urllib.request import Request, build_opener, HTTPError

SECONDS_PER_DAY = 3600.0 * 24.0

def download_from_jaxa(
    base_url='http://darts.jaxa.jp/planet/seismology/apollo/dump/dump/lp/START_TIME/',
    base_dir='', start_time=UTCDateTime('1969-07-21T03:00:00.000000Z'),
    end_time=UTCDateTime('1977-09-30T:21:00.000000Z')):

    time_interval = timedelta(hours=3)
    start = start_time
    end = start_time + time_interval

    while start < end_time:

        # in the last time inteval, download only to the end time
        if end > end_time:
            end = end_time

        # name csv outputfile with the start time
        csv_file = '%s.csv.gz' % start.strftime("%Y-%m-%dT%H:%M:%S")

        # build url for the Jaxa website
        url = '{}{}/STOP_TIME/{}'.format(base_url,
          start.strftime("%Y-%m-%dT%H:%M:%S"),
          end.strftime("%Y-%m-%dT%H:%M:%S"))
        print (url)

        opener = build_opener()
        code, data = download_url(url, opener,timeout=180, debug=False,
          return_string=True, use_gzip=True)

        with gzip.open(os.path.join(base_dir, csv_file), "wb") as gzip_file:
            gzip_file.write(data)

        # increment the time interval
        start += time_interval
        end += time_interval

    print('Done')

def download_url(url, opener, timeout=10, debug=False,
                 return_string=True, use_gzip=True):
    """
    Returns a pair of tuples.

    The first one is the returned HTTP code and the second the data as
    string.

    Will return a tuple of Nones if the service could not be found.
    All encountered exceptions will get raised unless `debug=True` is
    specified.

    Performs an http GET
    """
    if debug is True:
        print("Downloading %s %s requesting gzip compression" % (
            url, "with" if use_gzip else "without"))

    try:
        request = Request(url=url)
        # Request gzip encoding if desired.
        if use_gzip:
            request.add_header("Accept-encoding", "gzip")

        url_obj = opener.open(request, timeout=timeout)
    # Catch HTTP errors.
    except HTTPError as e:
        if debug is True:
            msg = "HTTP error %i, reason %s, while downloading '%s': %s" % \
                  (e.code, str(e.reason), url, e.read())
            print(msg)
        return e.code, e

    code = url_obj.getcode()

    # Unpack gzip if necessary.
    if url_obj.info().get("Content-Encoding") == "gzip":
        if debug is True:
            print("Uncompressing gzipped response for %s" % url)
        # Cannot directly stream to gzip from urllib!
        # http://www.enricozini.org/2011/cazzeggio/python-gzip/
        buf = io.BytesIO(url_obj.read())
        buf.seek(0, 0)
        f = gzip.GzipFile(fileobj=buf)
    else:
        f = url_obj


    if return_string is False:
        data = io.BytesIO(f.read())
    else:
        data = f.read()

    if debug is True:
        print("Downloaded %s with HTTP code: %i" % (url, code))

    return code, data
