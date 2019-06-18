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
from pdart.diffusion.view_catalog_with_envelopes import view_catalog_with_envelopes


if __name__ == "__main__":

    top_level_processed_dir = '/Users/nunn/Google Drive/for_Britt/PDART_PROCESSED'
    inv_name = "/Users/nunn/Google Drive/for_Britt/inventory/XA.1969-1977.xml"
    pre_filt_main = [0.2,0.3,0.9,1.3]
    pre_filt_env = [[0.1,0.25,0.75,1],[0.5,0.75,1.25,1.5],[1,1.25,1.75,2.25]]


    print('''The Lognonne_2003.xml catalog contains measured longitude/latitudes 
        for the artificial impacts.
        The catolog contains estimates for most other events.''')
    view_catalog_with_envelopes(top_level_processed_dir,
      '/Users/nunn/lunar_data/lunar_catalogs/output/Lognonne_2003/Lognonne_2003.xml',inv_name,
     dir_type='processed_dir',pre_filt_main=pre_filt_main,
     pre_filt_env=pre_filt_env,output='ACC',smooth_periods=30,
     plot_type='with_seismogram')

    # original catalog
    # view_catalog_with_envelopes(top_level_processed_dir,
    #   '/Users/nunn/Google Drive/for_Britt/catalogs/S12_shallow.xml',inv_name,
    #  dir_type='processed_dir',pre_filt_main=pre_filt_main,
    #  pre_filt_env=pre_filt_env,output='ACC',smooth_periods=30,
    #  plot_type='with_seismogram')




