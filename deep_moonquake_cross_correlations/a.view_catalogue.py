#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

View the Catalog

:copyright:
    The PDART Development Team & Katja Heger & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

from obspy import read_events



def view_catalog(file):
    catalog = read_events(file)
    for ev in catalog:
        picks = ev.picks
        for pick in picks:
           # print(pick.time)
            print(pick)


if __name__ == "__main__":
    view_catalog('../LunarCatalog_Nakamura_1981_and_updates_v1/LunarCatalog_Nakamura_1981_and_updates_v1_A01.xml')
