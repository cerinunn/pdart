#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot the Moon.

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
# from __future__ import (absolute_import, division, print_function,
#                         unicode_literals,absolute_import, division, print_function)
# from future.builtins import *  # NOQA
import matplotlib
import matplotlib.pyplot as plt

def plot_moon(size='normal'):

    from netCDF4 import Dataset as netcdf_dataset
    import cartopy.crs as ccrs
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    if size == 'normal':
        # this is really irritating - you can't specify figure size if
        # you change the dpi
        mpl.rcParams['figure.dpi'] = 80
        mpl.rcParams['savefig.dpi'] = 200
        mpl.rcParams['figure.figsize'] = [3.3, 2.8]


        mpl.rcParams['font.size'] = 12
        mpl.rcParams['legend.fontsize'] = 'large'
        mpl.rcParams['figure.titlesize'] = 'medium'

        filename='../extra_plots_output/moon_topo_networkX.png'

    elif size == 'large':
        mpl.rcParams['figure.figsize'] = [11.69,8.27]
        mpl.rcParams['figure.dpi'] = 80
        mpl.rcParams['savefig.dpi'] = 200

        mpl.rcParams['font.size'] = 12
        mpl.rcParams['legend.fontsize'] = 'large'
        mpl.rcParams['figure.titlesize'] = 'medium'

        filename='../extra_plots_output/moon_topo_network_largeX.png'

    # dataset = netcdf_dataset('../extra_plots_input/lata_a.grd')
    dataset = netcdf_dataset('../extra_plots_input/lalt_topo_ver3.grd')

    elev = dataset.variables['z'][:].astype('float32')
    lats = dataset.variables['lat'][:].astype('float32')
    lons = dataset.variables['lon'][:].astype('float32')

    # print(elev.shape)
    # print(lons[15:-15])
    # print(elev[:,15:-15].max())
    # print(elev[:,15:-15].min())
    # exit()

    # max / min elevation (approximate)
    # 8.51977
    # -7.24387

    ax = plt.axes(projection=ccrs.Orthographic())
    im = ax.pcolormesh(lons, lats, elev, vmin=-7.25, vmax=8.5, cmap='jet', transform=ccrs.PlateCarree())
    # im = ax.pcolormesh(lons, lats, elev, vmin=-, vmax=-4, cmap='viridis', transform=ccrs.PlateCarree())


    # station_names1 = ['S11', 'S17']
    # station_lats1 = [ 0.67409, 20.18935]
    # station_lons1 = [23.47298, 30.76796]


    # station_names1 = ['S17']
    # station_lats1 = [ 20.18935]
    # station_lons1 = [30.76796]

    station_names = ['S11','S12', 'S14','S15','S16']
    station_lats = [  0.67409,-3.01084, -3.64450, 26.13407,-8.97577]
    station_lons = [23.47298,-23.42456, -17.47753, 3.62981,15.49649]


    # plot black inverted triangles
    # ax.plot(station_lons1, station_lats1, 'kv', markersize=10, transform=ccrs.PlateCarree())

    transform = ccrs.PlateCarree()._as_mpl_transform(ax)
    # for (i, txt) in enumerate(station_names1):
        # ax.annotate(txt, xy=(station_lons1[i]+1, station_lats1[i]+1), xycoords=transform,
        #             ha='left', va='bottom',size=18, color='#606060')

    # plot inverted triangles
    ax.plot(station_lons, station_lats, 'kv', markersize=10, transform=ccrs.PlateCarree())
    for (i, txt) in enumerate(station_names):

        if txt == 'S14' or txt =='S16':
            ax.annotate(txt, xy=(station_lons[i], station_lats[i]-5), xycoords=transform,
                        ha='center', va='top',size=18, color='k')
        elif txt == 'S15':
            ax.annotate(txt, xy=(station_lons[i], station_lats[i]+2), xycoords=transform,
                        ha='right', va='bottom',size=18, color='k')
        else:
            ax.annotate(txt, xy=(station_lons[i]-4, station_lats[i]), xycoords=transform,
                        ha='right', va='center',size=18, color='k')


    ax.set_global()

    cbar = plt.colorbar(im,orientation='vertical',label='km')
    cbar.ax.tick_params(labelsize=12)
    # need a png file - pdfs are just too slow

    plt.savefig(filename)
    print(filename)
    # plt.show()
    plt.close()

# #############################
# RUN this plot_moon
###############################

if __name__ == "__main__":
    plot_moon()
    # plot_moon('large')
