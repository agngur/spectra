#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Description:
# Script which convert all .fits (binary) files to .asc (text files) in current directory.
# Extract main table from SDSS 1d spectra .fits files and save in asci text file.
# Column of asci files are same as rows in .fits:
# lambda[A], Flux_lambda[erg/cm2/s/A], continuum_substracted_flux, Flux_error, Bitmask, 
# Sky_template (undocumented - only in the DR6 or DR7)???

# Attention:
# output: lambda[A], Flux_lambda[erg/cm2/s/A], wdisp (Flux error), model fit
#
# Required packages:
# - python >=3.5
# - numpy >=1.0.3
# - astropy
#

import math
import glob
import os
import numpy as np
import pylab
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Column, Table
#from scipy.stats import binned_statistic

# To change binning factor You have to change parameters above:
bin1 = 3
bin2 = 5

# Absorbtion and emission lines
lines=(3704.906, 3713.027, 3722.997, 3735.430, 3751.217, 3771.701, 3798.976, 3835.000,
        3836.472, 3890.151, 3934.777, 3969.591, 3971.195, 4102.892, 4227.9179, 4308.9556, 
        4341.684, 4862.683, 5174.1251, 5185.0479, 5891.5833, 5897.5581, 6564.61,
        6560.10, 5411.50, 4859.32, 4541.59, 4685.682, 3819.607, 3964.729, 4026.191, 4120.821, 4143.76, 4387.929, 4471.479, 4713.146, 4921.931, 5015.678, 5047.74, 5875.70, 6678.15, 7065.19, 7281.35)
#names=("H-16","H-15","H-14","H-13","H-12","H-11","H-10","MgI","H-9","H-8","CaII_K","CaII_H",
#        "H-e","H-d_h","CaI_g","CaI_G","H-g_f","H-b_F","MgI_b2","MgI_b1","NaI_D2","NaI_D1","H-a_C","HeII","HeII","HeII","HeII","HeII","HeI","HeI","HeI","HeI","HeI","HeI","HeI","HeI","HeI","HeI","HeI","HeI","HeI","HeI","HeI")
names=("H","","","","","","","MgI","","","CaII","",
        "","H-d","CaI","","H-g","H-b","MgI","","NaI","","H-a","HeI","HeI","HeI","HeI","",
        "HeII","","","","HeII","","","","","","","HeII","","","HeII")

colors=("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","red","darkgrey","darkgrey","navy","navy",
"darkgrey","lime","navy","navy","lime","lime","red","red","mediumorchid","mediumorchid","lime","aqua","aqua","aqua","aqua","aqua",
"orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange")

# Definition of function that calculates the flux measure value and the wavelength, depending on the binning factor and the original values of flux and wavelengths:
def binning_spectra(wavelength, flux, binning):
    n_bins = math.ceil(len(wavelength)/binning)-1
    
    n_lamb = []
    m_flu = []
    
# CHECKPOINT whether entered data are correct   
    if len(wavelength) != len(flux):
        sys.exit("Lambda must be the same length as flux!")
    
    else:
        index = len(wavelength)-1
        for i in range(index):
        
            while (index-i) >= binning:
                n=np.mean(wavelength[i:i+(binning-1)])
                m=np.mean(flux[i:i+(binning-1)])
                n_lamb.append(n) 
                m_flu.append(m)

# We only leave wavelength values up to 9500 A:                
                if wavelength[i] > 9500:
                    wavelength = np.delete(wavelength,i)
                    flu = np.delete(flu,i)
                    ferr = np.delete(ferr,i)
                    mod = np.delete(mod,i)
                i+=binning
        
            if ((index-i) == 0) | ((index-i) < 0):
                break
            else:   
# To get the last bin, which is without full coverage of flux measures - uncomment the following lines:
#                n=np.mean(wavelength[i:-1])
#                m=np.mean(flux[i:-1])
#                n_lamb.append(n)
#                m_flu.append(m)
                break
                
    return n_lamb, m_flu  
        

# All found spectra (.fits) in this directory:
widma = glob.glob(os.getcwd()+'/'+'*.fits')
widma.sort()

for pliksp in widma:
# Opening fits file and conversion of required columns
    hdulist = fits.open(pliksp)
    dane = hdulist[1].data
    hdr = hdulist[0].header
    ra = str(hdr['PLUG_RA'])
    dec = str(hdr['PLUG_DEC'])
    clas = str(hdulist[2].data['class'])
    subcl = str(hdulist[2].data['subclass'])
    z = str(hdulist[2].data['z']) 
    z = float(z[1:-1])  
    zerr = str(hdulist[2].data['z_err'])
    zerr = float(zerr[1:-1])
    vel = z*3e+5       #[km/s]
    name = os.path.split(pliksp)[-1]
# SDSS SkyServer url:
    link = "http://skyserver.sdss.org/dr14/en/tools/explore/summary.aspx?plate="+str(hdr['PLATEID'])+"&mjd="+str(hdr['MJD'])+"&fiber="+str(hdr['FIBERID'])
    #x = urllib.request.urlopen(link)
    print(link+"\n")
    hdulist.close()
    
    lam = dane['loglam']        # Wavelengths [m]
    lamb = np.power(10,lam)     # Wavelengths [A]
    flu = dane['flux']*1e-17    # Flux measure [erg/cm2/s/A]
    ferr = dane['wdisp']*1e-17  # Flux error [erg/cm2/s/A] or 'w_dispersion'
    mod = dane['model']*1e-17   # Model fit
    
# Converting data from fits file and saving as ascii file
# (but we will operate on '.fits' spectra file)
    dd = np.vstack((lamb,flu,ferr,mod))
    np.savetxt(pliksp[:-4]+'asc',np.transpose(dd))

# Spectral plots  
    tab = Table(data=(lamb,flu,ferr,mod), names=("lamb","flu","ferr","mod"))
    mask = tab["lamb"]<=9500

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312, sharex = ax1)
    ax3 = fig.add_subplot(313, sharex = ax1)
    
# Plot - without binning
    #plt.subplot(311)
    ax1.plot(tab["lamb"][mask],tab["flu"][mask], label="no binning")
    ax1.locator_params(tight=True, nbins=12, bottom="on")
    ax1.legend(loc=1, fontsize="x-small")
    ax1.set_xlim(xmin=3700, xmax=9600)
    ax1.tick_params('y', colors='k', labelbottom='off')
    ax1.set_title(str(name)+"\nRA="+ra+" DEC="+dec+' '+clas+subcl+"\n\n", fontsize=8)
    x_bounds = ax1.get_xlim()
    for i in range(0,len(lines)-20):
        ax1.axvline(x=lines[i], ymin=0.01, ymax=1, gid=names[i], lw=0.5, ls=":", c=colors[i]) 
        ax1.annotate(s=names[i], xy =(((lines[i]-x_bounds[0])/(x_bounds[1]-x_bounds[0])),1.04), xycoords='axes fraction', verticalalignment='right', horizontalalignment='bottom', rotation = 90, color=colors[i], fontsize='x-small')    

# Plot - with binning 'bin1'
    #plt.subplot(312)
    nn_lamb, nn_flu = binning_spectra(tab["lamb"][mask],tab["flu"][mask],bin1)
    ax2.plot(nn_lamb,nn_flu,label="binning "+str(bin1)+"pts")
    ax2.locator_params(tight=True, nbins=12, bottom="on")
    ax2.legend(loc=1, fontsize='x-small')
    ax2.set_xlim(xmin=3700, xmax=9600)
    ax2.tick_params(labelbottom='off')
    ax2.set_ylabel("Flux($\lambda$) [erg/cm$^2$/s/$\AA$]")
    for i in range(0,len(lines)):
        ax2.axvline(x=lines[i], ymin=0.01, ymax=1, gid=names[i], lw=0.5, ls=":", c=colors[i])

# Plot - with binning 'bin2'
    #plt.subplot(313)
    nnn_lamb, nnn_flu = binning_spectra(tab["lamb"][mask],tab["flu"][mask],bin2)
    ax3.plot(nnn_lamb,nnn_flu,label="binning "+str(bin2)+"pts")
    ax3.locator_params(tight=True, nbins=12, bottom="on")
    ax3.legend(loc=1,  fontsize="x-small")
    ax3.set_xlim(xmin=3700, xmax=9600)
    ax3.set_xlabel("Wavelength [$\AA$]"+"\nz="+str(round(z,4))+"±"+str(round(zerr,4))+"    RV="+str(vel)+" km/s")
    for i in range(len(lines)-20,len(lines)):
        ax3.axvline(x=lines[i], ymin=0.01, ymax=1, gid=names[i], lw=0.5, ls=":", c=colors[i])
        ax3.annotate(s=names[i], xy =(((lines[i]-x_bounds[0])/(x_bounds[1]-x_bounds[0])),0.04), xycoords='axes fraction', verticalalignment='right', horizontalalignment='bottom', rotation = 90, color=colors[i], fontsize='x-small')

# Position and visibility of yaxis exponential notation
    #ax1.yaxis.get_offset_text().set_visible(False)
    ax1.yaxis.set_offset_position("right")
    #ax2.yaxis.get_offset_text().set_visible(False)
    ax2.yaxis.set_offset_position("right")
    #ax3.yaxis.get_offset_text().set_visible(False)
    ax3.yaxis.set_offset_position("right")
    
    fig = plt.gcf()
    fig.set_size_inches(12, 6) #size of the window x,y
    plt.subplots_adjust(top=0.90, bottom=0.12, left=0.10, right=0.95, hspace=0.00, wspace=0.25)

    plt.show()
    
    #input("Enter ")
  
