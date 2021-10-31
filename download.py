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
import os
import io
import gzip

# from obspy.core.utcdatetime import UTCDateTime
from urllib.request import Request, build_opener, HTTPError
import httplib2
from bs4 import BeautifulSoup, SoupStrainer

# use environment env_obspy_downgrade

def create_links(url, match,base_dir):

    http = httplib2.Http()
    status, response = http.request(url)

    links = []
    for link in BeautifulSoup(response, parseOnlyThese=SoupStrainer('a')):

        if link.has_attr('href'):
            if link['href'][0:len(match)] == match:
                links.append(link['href'])

    out_filename = os.path.join(base_dir,'links.txt')
    with open(out_filename, 'w') as f:
        for link in links:
            f.write("%s\n" % link)

def read_links(base_url,match,base_dir):

    links = []
    in_filename = os.path.join(base_dir,'links.txt')

    # if file exists already, we can use it
    if not os.path.isfile(in_filename):
        create_links(base_url,match,base_dir)

    with open(in_filename, 'rt') as f:
        for line in f:
            links.append(line.rstrip('\n'))

    return links


def download_from_jaxa(base_url,match,base_dir):

    links = read_links(base_url,match,base_dir)
    for filename in links:

        url = "{}/{}".format(base_url,filename)

        opener = build_opener()
        code, data = download_url(url, opener,timeout=180, debug=False,
          return_string=True, use_gzip=True)

        filename = '{}.gz'.format(filename)
        with gzip.open(os.path.join(base_dir, filename), "wb") as gzip_file:
            gzip_file.write(data)

        # with open(os.path.join(base_dir, filename), "wb") as filewrite:
        #     filewrite.write(data)

        # exit()


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


# if __name__ == "__main__":
    # url='http://darts.jaxa.jp/pub/apollo/pse/p11s/'
    # match='pse.a11'
    # base_dir='/Users/nunn/lunar_data/PDART_TAPES/S11'
    #
    # download_from_jaxa(url,match,base_dir)

    # url='http://darts.jaxa.jp/pub/apollo/pse/p12s/'
    # match='pse.a12'
    # base_dir='/Users/nunn/lunar_data/PDART_TAPES/S12'
    #
    # download_from_jaxa(url,match,base_dir)
    #
    # url='http://darts.jaxa.jp/pub/apollo/pse/p14s/'
    # match='pse.a14'
    # base_dir='/Users/nunn/lunar_data/PDART_TAPES/S14'
    #
    # download_from_jaxa(url,match,base_dir)
    #
    # url='http://darts.jaxa.jp/pub/apollo/pse/p15s/'
    # match='pse.a15'
    # base_dir='/Users/nunn/lunar_data/PDART_TAPES/S15'
    #
    # download_from_jaxa(url,match,base_dir)
    #
    # url='http://darts.jaxa.jp/pub/apollo/pse/p16s/'
    # match='pse.a16'
    # base_dir='/Users/nunn/lunar_data/PDART_TAPES/S16'
    #
    # url='http://darts.jaxa.jp/pub/apollo/pse/wtns/'
    # match='wtn'
    # base_dir='/Users/nunn/lunar_data/PDART_TAPES/work_tapes'
    #
    # download_from_jaxa(url,match,base_dir)
