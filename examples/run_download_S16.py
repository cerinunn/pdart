#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example file to run the download.

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from pdart.download import download_from_jaxa, create_links, read_links

# use environment env_obspy_downgrade

if __name__ == "__main__":

    url='http://darts.jaxa.jp/pub/apollo/pse/p16s/'
    match='pse.a16'
    base_dir='/Users/nunn/lunar_data/PDART_TAPES/S16'

    download_from_jaxa(url,match,base_dir)
