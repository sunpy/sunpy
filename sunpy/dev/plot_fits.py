#-*- coding:utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
# Author: Keith Hughitt <keith.hughitt@nasa.gov>
# Written: 2011/03/29
#
# <License info will go here...>
#
"""Test function for plotting an AIA FITS image."""

import os
import datetime
import pyfits
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.patches import Circle
from sunpy import Sun
import sunpy.dev.cm

AIA_SAMPLE_IMAGE = 'sample-data/AIA20110319_105400_0171.fits'

def plot_fits(filepath=None):
    '''Plots an AIA image.'''

    if filepath is None:
        filepath = os.path.join(os.path.dirname(__file__), AIA_SAMPLE_IMAGE)

    # Load fits file
    fits = pyfits.open(filepath)
    header = fits[0].header
    data = fits[0].data

    # Get useful header information
    date  = header['date-obs']
    fitsDatetime = datetime.datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f")
    instr = header['instrume']
    rSun  = header['r_sun']
    wavelength = header['wavelnth']
    
    centerX = header['crpix1']
    centerY = header['crpix1']
    scaleX = header['cdelt1']
    scaleY = header['cdelt2']
    
    # Create a figure and add title and axes
    fig = plt.figure()
    
    ax = fig.add_subplot(111)
    ax.set_title(instr + ' ' + date)
    ax.set_xlabel('X-postion (arcseconds)')
    ax.set_ylabel('Y-postion (arcseconds)')
    
    # Draw circle at solar limb
    circ = Circle([0, 0], radius=Sun.radius(fitsDatetime), fill=False, color='white')
    ax.add_artist(circ)
        
    # Determine extent
    xmin = -(centerX - 1) * scaleX
    xmax = (centerX - 1) * scaleX
    ymin = -(centerY - 1) * scaleY
    ymax = (centerY - 1) * scaleY
    
    extent = [xmin, xmax, ymin, ymax]
    
    # Draw image
    cm = sunpy.dev.cm.log_adaptive(data)
    imgplot = plt.imshow(data, cmap=cm, origin='lower', extent=extent)

    plt.colorbar()
    
    # Show figure
    plt.show()

