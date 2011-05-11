#-*- coding:utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
# Author: Keith Hughitt <keith.hughitt@nasa.gov>
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
#from sunpy import Sun
import numpy as np

AIA_SAMPLE_IMAGE = 'doc/sample-data/AIA20110319_105400_0171.fits'

def plot_fits(filepath=None):
    '''Plots an AIA image.'''

    if filepath is None:
        filepath = os.path.join(os.path.dirname(__file__), AIA_SAMPLE_IMAGE)

    # Load fits file
    fits = pyfits.open(filepath)
    header = fits[0].header
    data = fits[0].data

    #show all of the header items
    # header.items

    # Get useful header information
    date  = header.get('date-obs')
    #fitsDatetime = datetime.datetime(date[0:4], date[5:7], date[8:10], date[11:13], date[14:16], date[17:19])
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
    ax.set_title(instr + ' ' + str(wavelength) + '$\AA$' + ' ' + date)
    ax.set_xlabel('X-postion (arcseconds)')
    ax.set_ylabel('Y-postion (arcseconds)')
    
    # Draw circle at solar limb
    #circ = Circle([0, 0], radius=Sun.radius(fitsDatetime), fill=False, color='white')
    #ax.add_artist(circ)
        
    # Determine extent
    xmin = -(centerX - 1) * scaleX
    xmax =  (centerX - 1) * scaleX
    ymin = -(centerY - 1) * scaleY
    ymax =  (centerY - 1) * scaleY
    
    extent = [xmin, xmax, ymin, ymax]
    
    color_map = aia_color_table(wavelength)
    
    # Draw image
    #imgplot = plt.imshow(data, cmap=cm.Greys_r, origin='lower', extent=extent)
    imgplot = plt.imshow(data, cmap=color_map, origin='lower', extent=extent)

    plt.colorbar()
    
    # Show figure
    plt.show()

def aia_color_table(wavelength = None):
    '''Returns a color table for SDO AIA images'''
    # Based on aia_lct.pro by Karl Schriver (2010/04/12)

    # The following values describe color table 3 for IDL (Red Temperature)
    r0 = np.array([0,1,2,4,5,7,8,10,11,13,14,15,17,18,20,21,23,24,26,27,28,30,31,33,34,36,37,39,40,42,43,44,46,47,49,50,52,53,55,56,57,59,60,62,63,65,66,68,69,70,72,73,75,76,78,79,81,82,84,85,86,88,89,91,92,94,95,97,98,99,101,102,104,105,107,108,110,111,113,114,115,117,118,120,121,123,124,126,127,128,130,131,133,134,136,137,139,140,141,143,144,146,147,149,150,152,153,155,156,157,159,160,162,163,165,166,168,169,170,172,173,175,176,178,179,181,182,184,185,186,188,189,191,192,194,195,197,198,199,201,202,204,205,207,208,210,211,212,214,215,217,218,220,221,223,224,226,227,228,230,231,233,234,236,237,239,240,241,243,244,246,247,249,250,252,253,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255])
    g0 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3,5,7,9,11,13,15,17,18,20,22,24,26,28,30,32,34,35,37,39,41,43,45,47,49,51,52,54,56,58,60,62,64,66,68,69,71,73,75,77,79,81,83,85,86,88,90,92,94,96,98,100,102,103,105,107,109,111,113,115,117,119,120,122,124,126,128,130,132,134,136,137,139,141,143,145,147,149,151,153,154,156,158,160,162,164,166,168,170,171,173,175,177,179,181,183,185,187,188,190,192,194,196,198,200,202,204,205,207,209,211,213,215,217,219,221,222,224,226,228,230,232,234,236,238,239,241,243,245,247,249,251,253,255])
    b0 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,7,11,15,19,23,27,31,35,39,43,47,51,54,58,62,66,70,74,78,82,86,90,94,98,102,105,109,113,117,121,125,129,133,137,141,145,149,153,156,160,164,168,172,176,180,184,188,192,196,200,204,207,211,215,219,223,227,231,235,239,243,247,251,255])
        
    c0 = np.arange(256, dtype = 'f')
    c1 = (np.sqrt(c0)*np.sqrt(255.0))
    c2 = (np.arange(256)**2/255.0)
    c3 = ((c1+c2/2.0)*255.0/(c1.max() + c2.max()/2.0))

    if wavelength == 1600: 
        r=c3
        g=c3
        b=c2
    elif wavelength == 1700:
        r=c1
        g=c0
        b=c0
    elif wavelength == 4500:
        r=c0
        g=c0
        b=b0/2.0
    elif wavelength == 94:
        r=c2
        g=c3
        b=c0
    elif wavelength == 131:
        r=g0
        g=r0
        b=r0
    elif wavelength == 171:
        r=r0
        g=c0
        b=b0
    elif wavelength == 193:
        r=c1
        g=c0
        b=c2
    elif wavelength == 211:
        r=c1
        g=c0
        b=c3
    elif wavelength == 304:
        r=r0
        g=g0
        b=b0 
    elif wavelength == 335:
        r=c2
        g=c0
        b=c1
    else:
        print("Please choose a valid AIA wavelength (1600,1700,4500,94,131,171,193,211,304,335).")
        return None
   
    # Now create the color tuples
    i = np.arange(r0.size,dtype = 'f')/r0.size
    
    rtuple = list(zip(i,r/255.0,r/255.0))
    rtuple.append((1.0,r[-1],r[-1]))
    gtuple = list(zip(i,g/255.0,g/255.0))
    gtuple.append((1.0,g[-1]/255.0,g[-1]/255.0))
    btuple = list(zip(i,b/255.0,b/255.0))
    btuple.append((1.0,b[-1]/255.0,b[-1]/255.0))
    
    cdict = {'red':rtuple, 'green':gtuple, 'blue': btuple}    
    color_table = colors.LinearSegmentedColormap('mytable', cdict)

    return color_table