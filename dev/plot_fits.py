#!/usr/bin/env python
#-*- coding:utf-8 -*-

#
# Plots an AIA image
#
import sys
import pyfits
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.patches import Circle
import thesun
import datetime

def main():
    # Load fits file
    fits = pyfits.open('AIA20110319_105400_0171.fits')
    header = fits[0].header
    data = fits[0].data

    # Get useful header information
    date  = header['date-obs']
    fitsDatetime = datetime.datetime(date[0:4],date[5:7],date[8:10],date[11:13],date[14:16],date[17:19])
    instr = header['instrume']
    rSun  = header['r_sun']
    wavelength = header['wavelnth']
    
    centerX = header['crpix1']
    centerY = header['crpix1']
    scaleX  = header['cdelt1']
    scaleY  = header['cdelt2']
    
    # Create a figure and add title and axes
    fig = plt.figure()
    
    ax = fig.add_subplot(111)
    ax.set_title(instr + ' ' + date)
    ax.set_xlabel('X-postion (arcseconds)')
    ax.set_ylabel('Y-postion (arcseconds)')
    
    # Draw circle at solar limb
    circ = Circle([0,0], radius = thesun.radius(fitsDatetime), fill = False, color = 'white')
    ax.add_artist(circ)
        
    # Determine extent
    xmin = -(centerX -1) * scaleX
    xmax =  (centerX -1) * scaleX
    ymin = -(centerY -1) * scaleY
    ymax =  (centerY -1) * scaleY
    
    extent = [xmin, xmax, ymin, ymax]
    
    # Draw image
    imgplot = plt.imshow(data, cmap=cm.Greys_r, origin='lower', extent=extent)

    plt.colorbar()
    
    # Show figure
    plt.show()

if __name__ == '__main__':
    sys.exit(main())