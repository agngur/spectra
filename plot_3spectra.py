#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Description:
# The following script displays a window containing simultaneously three consecutive spectral plots (SDSS) on the screen. 
# After closing this window, in each subsequent one, the second and third spectrum from the previous window 
# will be displayed, and the new one - the next spectrum in order in the working directory.
#
# Required packages:
# - python >=3.5
# - numpy >=1.0.3
# - astropy
# - matplotlib


import math
import glob
import os
import numpy as np
#import pylab
from matplotlib import pyplot as plt
from astropy.io import fits

# Absorbtion and emission lines
line =[[3704.906, 'H', 'darkgrey'],
       [3713.027, '', 'darkgrey'],
       [3722.997, '', 'darkgrey'],
       [3735.43, '', 'darkgrey'],
       [3751.217, '', 'darkgrey'],
       [3771.701, '', 'darkgrey'],
       [3798.976, '', 'darkgrey'],
       [3835.0, 'MgI', 'red'],
       [3836.472, '', 'darkgrey'],
       [3890.151, '', 'darkgrey'],
       [3934.777, 'CaII', 'navy'],
       [3969.591, '', 'navy'],
       [3971.195, '', 'darkgrey'],
       [4102.892, 'H-d', 'lime'],
       [4227.9179, 'CaI', 'navy'],
       [4308.9556, '', 'navy'],
       [4341.684, 'H-g', 'lime'],
       [4862.683, 'H-b', 'lime'],
       [5174.1251, 'MgI', 'red'],
       [5185.0479, '', 'red'],
       [5891.5833, 'NaI', 'mediumorchid'],
       [5897.5581, '', 'mediumorchid'],
       [6564.61, 'H-a', 'lime'],
       [6560.1, '', 'aqua'],
       [5411.5, 'HeI', 'aqua'],
       [4859.32, '', 'aqua'],
       [4541.59, 'HeI', 'aqua'],
       [4685.682, '', 'aqua'],
       [3819.607, '', 'orange'],
       [3964.729, '', 'orange'],
       [4026.191, '', 'orange'],
       [4120.821, '', 'orange'],
       [4143.76, 'HeII', 'orange'],
       [4387.929, '', 'orange'],
       [4471.479, '', 'orange'],
       [4713.146, '', 'orange'],
       [4921.931, '', 'orange'],
       [5015.678, '', 'orange'],
       [5047.74, '', 'orange'],
       [5875.7, '', 'orange'],
       [6678.15, 'HeII', 'orange'],
       [7065.19, '', 'orange'],
       [7281.35, 'HeII', 'orange']]

col = 20*["g","r","y"]

def plot_ref():
    ax.set_ylabel("Flux($\lambda$) [erg/cm$^2$/s/$\AA$]")
    ax.set_xlabel("Wavelength [$\AA$]")
    ax.locator_params(nbins=16)
    ax.yaxis.set_offset_position("right")
    x_bounds = ax.get_xlim()
    ax.set_xlim(xmin=3500, xmax=9600)
    for r in range(0,len(line)):
        plt.axvline(x=line[r][0], ymin=0.02, ymax=1, gid=line[r][1], lw=0.5, ls=":", c=line[r][2])
        plt.annotate(s=line[r][1], xy =(((line[r][0]-x_bounds[0])/(x_bounds[1]-x_bounds[0])),1.01), xycoords='axes fraction', verticalalignment='right', horizontalalignment='bottom', rotation = 90, color=line[r][2], fontsize='x-small')
    fig = plt.gcf()
    plt.legend(loc=1,  fontsize="x-small")
    fig.set_size_inches(10, 5)
    plt.show()
    
    
widma = glob.glob(os.getcwd()+'/'+'*.fits')
widma.sort()

if len(widma) <= 3:
    fig, ax = plt.subplots()

    for i in range(len(widma)):
        d, h = fits.getdata(widma[i], header=True)
        n = os.path.split(widma[i])[-1]
        lam = d['loglam']
        lamb = np.power(10,lam)
        flu = d['flux']*1e-17

        plt.plot(lamb,flu, label=str(n), color=col[i])
    plot_ref()

else:
    for i in range(len(widma)):
        spectra = widma[i:i+3]
        
        if len(spectra) < 3:
            break
        else:
            fig, ax = plt.subplots()

            for j in range(len(spectra)):
                d, h = fits.getdata(spectra[j], header=True)
                n = os.path.split(spectra[j])[-1]
                lam = d['loglam']
                lamb = np.power(10,lam)
                flu = d['flux']*1e-17        
              
                plt.plot(lamb, flu, label=str(n), color=col[i+j])
            plot_ref()

    
