#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot the response files for S14

:copyright:
    The PDART Development Team & Ceri Nunn
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import print_function
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math import pi

from obspy.core.inventory.inventory import read_inventory
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.inventory.response import _pitick2latex

def plot_seismometer_response():

    inv = read_inventory("/Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0.xml")
    # select times, so that we only get one copy of each of the response curves
    # peaked mode
    matching1 = inv.select(time=UTCDateTime('1973-01-01'), channel='MHZ')
    # flat mode
    matching2 = inv.select(time=UTCDateTime('1976-09-19'), channel='MHZ')
    # SHZ
    matching3 = inv.select(time=UTCDateTime('1974-10-17'), channel='SHZ')

    inv = matching1 + matching2 + matching3

    station='S14'
    # inv.plot_response(0.001, station=station,  output="DISP", outfile='S14_response1.png')
    fig = inv.plot_response(0.01, station=station,  output="DISP", unwrap_phase=True, show=False)
    # make some adjustments to the plot afterwards
    ax1, ax2 = fig.axes[:2]
    ax1.set_ylabel('Amplitude [DU/m]')
    for child in ax1.get_children():
        # remove the arrow annotations
        if isinstance(child, matplotlib.text.Annotation):
            child.remove()
    # adjust the size - now the annotations are not used
    fig.subplots_adjust(top=0.95, right=0.95)
    # make the labels a bit clearer for this plot
    handles, labels = ax1.get_legend_handles_labels()
    legend = ax1.legend()
    legend.remove()
    labels = ['Peaked Mode', 'Flat Mode', 'Short Period']
    print('BE CAREFUL - check that the colours correspond to the legend !!!!')
    print('Also note - the scale might need changing slightly - this used to work, but for some reason it is missing 1pi now')
    ax1.legend(handles,labels,loc="lower center", ncol=3, fontsize=14)
    plt.xlim(xmax=30,xmin=0.01)

    minmax2 = ax2.yaxis.get_data_interval()
    yticks2 = np.arange(minmax2[0] - minmax2[0] % (pi),
                        minmax2[1] - minmax2[1] % (pi) + pi, pi)
    ax2.set_yticks(yticks2)
    yticks2_minor = np.arange(minmax2[0] - minmax2[0] % (pi / 2),
                        minmax2[1] - minmax2[1] % (pi / 2) + pi, pi / 2)
    ax2.set_yticks(yticks2_minor, minor=True)

    ytick_labels2 = [_pitick2latex(x) for x in yticks2]
    ax2.set_yticklabels(ytick_labels2, fontsize=16)
    ax2.grid(True)
    ax2.set_ylim(-16.56067583,   1.57)

    plt.savefig('../extra_plots_output/S14_response_dispXX.pdf')
    plt.show()

# def plot_seismometer_response_disp():
# 
#     output='DISP'
# 
#     inv = read_inventory("/Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0.xml")
#     # select times, so that we only get one copy of each of the response curves
#     # peaked mode
#     matching1 = inv.select(time=UTCDateTime('1973-01-01'), channel='MHZ', location='00')
#     # flat mode
#     matching2 = inv.select(time=UTCDateTime('1976-09-19'), channel='MHZ', location='01')
#     # SHZ
#     matching3 = inv.select(time=UTCDateTime('1974-10-17'), channel='SHZ')
# 
#     inv = matching1 + matching2 + matching3
# 
#     # inv = matching1
# 
#     station='S14'
#     # inv.plot_response(0.001, station=station,  output="DISP", outfile='S14_response1.png')
#     fig = inv.plot_response(0.01, station=station,  output=output, unwrap_phase=True, show=False)
#     # make some adjustments to the plot afterwards
#     ax1, ax2 = fig.axes[:2]
#     if output == 'DISP':
#         ax1.set_ylabel('Amplitude in Displacement\n[DU/m]')
#     elif output ==  'VEL':
#         ax1.set_ylabel('Amplitude in Velocity\n[DU/m/s]')
#     elif output ==  'ACC':
#         ax1.set_ylabel('Amplitude in Acceleration\n[DU/m/s^2]')
#     for child in ax1.get_children():
#         # remove the arrow annotations
#         if isinstance(child, matplotlib.text.Annotation):
#             child.remove()
#     # adjust the size - now the annotations are not used
#     fig.subplots_adjust(top=0.95, right=0.95)
#     # make the labels a bit clearer for this plot
#     handles, labels = ax1.get_legend_handles_labels()
#     legend = ax1.legend()
#     legend.remove()
#     labels = ['Peaked Mode', 'Flat Mode', 'Short Period']
#     print('BE CAREFUL - check that the colours correspond to the legend !!!!')
#     ax1.legend(handles,labels,loc="lower center", ncol=3, fontsize=14)
#     plt.xlim(xmax=30, xmin=0.01)
# 
#     # plt.savefig('S14_response_acc.png')
#     # plt.savefig('S14_response_acc.pdf')
#     plt.show()


# def plot_seismometer_response_acc():
# 
#     inv = read_inventory("/Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0.xml")
#     # select times, so that we only get one copy of each of the response curves
#     # peaked mode
#     matching1 = inv.select(time=UTCDateTime('1973-01-01'), channel='MHZ')
#     # flat mode
#     matching2 = inv.select(time=UTCDateTime('1976-09-19'), channel='MHZ')
#     # SHZ
#     matching3 = inv.select(time=UTCDateTime('1974-10-17'), channel='SHZ')
# 
#     inv = matching1 + matching2 + matching3
# 
#     inv = matching1
# 
#     station='S14'
#     # inv.plot_response(0.001, station=station,  output="DISP", outfile='S14_response1.png')
#     fig = inv.plot_response(0.01, station=station,  output="ACC", unwrap_phase=True, show=False)
#     # make some adjustments to the plot afterwards
#     ax1, ax2 = fig.axes[:2]
#     ax1.set_ylabel('Amplitude [DU/m/s^2]')
#     for child in ax1.get_children():
#         # remove the arrow annotations
#         if isinstance(child, matplotlib.text.Annotation):
#             child.remove()
#     # adjust the size - now the annotations are not used
#     fig.subplots_adjust(top=0.95, right=0.95)
#     # make the labels a bit clearer for this plot
#     handles, labels = ax1.get_legend_handles_labels()
#     legend = ax1.legend()
#     legend.remove()
#     labels = ['Peaked Mode', 'Flat Mode', 'Short Period']
#     print('BE CAREFUL - check that the colours correspond to the legend !!!!')
#     ax1.legend(handles,labels,loc="lower center", ncol=3, fontsize=14)
#     plt.xlim(xmax=30)
# 
#     plt.savefig('S14_response_acc.png')
#     # plt.savefig('S14_response_acc.pdf')
#     # plt.show()

# def test_plot_response():
#     # from obspy import read_inventory
#     # resp = read_inventory()[0][0][0].response
#     # resp.plot(0.001, output="DISP")
#     inv = read_inventory("/Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0.xml")
#     resp = inv.select(time=UTCDateTime('1973-01-01'), channel='MHZ')[0][0][0].response
# 
#     resp.plot(0.001, output="DISP",unwrap_phase=True)
# 
# def plot_seismometer_response_all_test():
# 
#     inv = read_inventory("/Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0.xml")
#     # select times, so that we only get one copy of each of the response curves
# 
# 
#     # check the following in turn 
# 
#     # matching1 = inv.select( channel='MH1',station='S11')
#     # matching1 = inv.select( channel='MH2',station='S11')
#     # matching1 = inv.select( channel='MHZ',station='S11')
#     # matching1 = inv.select( channel='SHZ',station='S11')
#     # matching1 = inv.select(station='S11')
# 
#     # matching1 = inv.select( channel='MH1',station='S12')
#     # matching1 = inv.select( channel='MH2',station='S12')
#     # matching1 = inv.select( channel='MHZ',station='S12')
#     # matching1 = inv.select( channel='SHZ',station='S12')
#     matching1 = inv.select(station='S12')
# 
#     # matching1 = inv.select( channel='MH1',station='S14')
#     # matching1 = inv.select( channel='MH2',station='S14')
#     # matching1 = inv.select( channel='MHZ',station='S14')
#     # matching1 = inv.select( channel='SHZ',station='S14')
#     # matching1 = inv.select(station='S14')
# 
#     # matching1 = inv.select( channel='MH1',station='S15')
#     # matching1 = inv.select( channel='MH2',station='S15')
#     # matching1 = inv.select( channel='MHZ',station='S15')
#     # matching1 = inv.select( channel='SHZ',station='S15')
#     # matching1 = inv.select(station='S15')
# 
#     # matching1 = inv.select( channel='MH1',station='S16')
#     # matching1 = inv.select( channel='MH2',station='S16')
#     # matching1 = inv.select( channel='MHZ',station='S16')
#     # matching1 = inv.select( channel='SHZ',station='S16')
#     # matching1 = inv.select(station='S16')
# 
# 
#     # peaked mode
#     # matching1 = inv.select(time=UTCDateTime('1973-01-01'), channel='MHZ', location='00')
#     # matching2 = inv.select(time=UTCDateTime('1973-01-01'), channel='MH1', location='00')
#     # matching3 = inv.select(time=UTCDateTime('1973-01-01'), channel='MH2', location='00')
#     # # flat mode
#     # matching2 = inv.select(time=UTCDateTime('1976-09-19'), channel='MHZ')
#     # # SHZ
#     # matching3 = inv.select(time=UTCDateTime('1974-10-17'), channel='SHZ')
# 
# 
#     # inv = matching1 + matching2 + matching3
#     inv = matching1
#     for net in inv:
#         for sta in net:
#             for cha in net:
#                 for x in cha:
#                     print(x)
#                 # print(cha)
#     # print(inv)
#     # exit()
# 
#     # inv.plot_response(0.001, station=station,  output="DISP", outfile='S14_response1.png')
#     fig = inv.plot_response(0.01, output="DISP", unwrap_phase=True, show=False)
#     # make some adjustments to the plot afterwards
#     ax1, ax2 = fig.axes[:2]
#     ax1.set_ylabel('Amplitude [DU/m]')
#     for child in ax1.get_children():
#         # remove the arrow annotations
#         if isinstance(child, matplotlib.text.Annotation):
#             child.remove()
#     # adjust the size - now the annotations are not used
#     fig.subplots_adjust(top=0.95, right=0.95)
#     # make the labels a bit clearer for this plot
#     handles, labels = ax1.get_legend_handles_labels()
#     legend = ax1.legend()
#     # legend.remove()
#     # labels = ['Peaked Mode', 'Flat Mode', 'Short Period']
#     # print('BE CAREFUL - check that the colours correspond to the legend !!!!')
#     # ax1.legend(handles,labels,loc="lower center", ncol=3, fontsize=14)
#     plt.xlim(xmax=30)
# 
#     minmax2 = ax2.yaxis.get_data_interval()
#     yticks2 = np.arange(minmax2[0] - minmax2[0] % (pi),
#                         minmax2[1] - minmax2[1] % (pi) + pi, pi)
#     ax2.set_yticks(yticks2)
#     yticks2_minor = np.arange(minmax2[0] - minmax2[0] % (pi / 2),
#                         minmax2[1] - minmax2[1] % (pi / 2) + pi, pi / 2)
#     ax2.set_yticks(yticks2_minor, minor=True)
# 
#     ytick_labels2 = [_pitick2latex(x) for x in yticks2]
#     ax2.set_yticklabels(ytick_labels2, fontsize=16)
#     ax2.grid(True)
#     ax2.set_ylim(-16.56067583,   1.57)
# 
#     # plt.savefig('S14_response_XX_XX.png')
#     plt.savefig('S14_responseXX_wrapped.pdf')
#     plt.show()
# 

if __name__ == "__main__":
    # plot_seismometer_response()
    # plot_seismometer_response_acc()
    # test_plot_response()
    # plot_seismometer_response_all_test()
    # plot_seismometer_response_disp()
    plot_seismometer_response()
