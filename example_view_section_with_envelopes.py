#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
View the catalog with envelopes 

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)



"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA
from pdart.diffusion.view_catalog_with_envelopes import (
  view_catalog_with_envelopes, view_section_with_envelopes)


if __name__ == "__main__":

    top_level_processed_dir = '/Users/nunn/Google Drive/for_Britt/PDART_PROCESSED'
    inv_name = "/Users/nunn/Google Drive/for_Britt/inventory/XA.1969-1977.xml"
    # Garcia 2011 filters between 0.3 and 0.9 
    # pre_filt_main = [0.2,0.3,0.9,1.3]
    # trying these filtering params for the seismogram behind the envelopes 
    pre_filt_main = [0.1,0.25,1.75,2.25]
    pre_filt_env = [[0.1,0.25,0.75,1],[0.5,0.75,1.25,1.5],[1,1.25,1.75,2.25]]


    view_section_with_envelopes(top_level_processed_dir,
      '/Users/nunn/lunar_data/lunar_catalogs/output/Lognonne_2003/Lognonne_2003.xml',inv_name,
     dir_type='processed_dir',pre_filt_main=pre_filt_main,
     pre_filt_env=pre_filt_env,output='ACC',smooth_periods=5,
     scale=10,xlim=(0,150),ylim=(0,1400),
     title='Artificial Impacts',start_no=0,end_no=None)
