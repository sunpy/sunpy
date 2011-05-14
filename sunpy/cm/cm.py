"""
This module provides a set of colormaps specific to solar data (e.g. SDO/AIA color maps), 
functions for getting a colormap by name.

"""

import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import pylab
import matplotlib.cbook as cbook

def __aia_color_table__(wavelength = None):
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

sdoaia94 = __aia_color_table__(wavelength = 94)
sdoaia131 = __aia_color_table__(wavelength = 131)
sdoaia171 = __aia_color_table__(wavelength = 171)
sdoaia193 = __aia_color_table__(wavelength = 193)
sdoaia211 = __aia_color_table__(wavelength = 211)
sdoaia304 = __aia_color_table__(wavelength = 304)
sdoaia335 = __aia_color_table__(wavelength = 335)
sdoaia1600 = __aia_color_table__(wavelength = 1600)
sdoaia1700 = __aia_color_table__(wavelength = 1700)
sdoaia4500 = __aia_color_table__(wavelength = 4500)

cmlist = {
          'sdoaia94': sdoaia94,
          'sdoaia131': sdoaia131,
          'sdoaia171': sdoaia171,
          'sdoaia193': sdoaia193,
          'sdoaia211': sdoaia211,
          'sdoaia304': sdoaia304,
          'sdoaia335': sdoaia335,
          'sdoaia1600': sdoaia1600,
          'sdoaia1700': sdoaia1700,
          'sdoaia4500': sdoaia4500,
          }

def get_cmap(name=None, lut=None):
    """
    Get a colormap instance.

    """
    if name is None:
        name = 'sdoaia94'

    if name in cmlist:
        return cmlist.get(name)
    else:
        raise ValueError("Colormap %s is not recognized" % name)

def show_colormaps(cm = None):

    maps = sorted(m for m in cmlist)
    nmaps = len(maps) + 1
    
    a = np.linspace(0, 1, 256).reshape(1,-1)
    a = np.vstack((a,a))
    
    fig = plt.figure(figsize=(5,10))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    for i,name in enumerate(maps):
        ax = plt.subplot(nmaps, 1, i+1)
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=get_cmap(name), origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], name, fontsize=10, horizontalalignment='right')

    plt.show()


def test_equalize(data):
    '''Test'''

    dfile = cbook.get_sample_data('s1045.ima', asfileobj=False)

    im = np.fromstring(file(dfile, 'rb').read(), np.uint16).astype(float)
    im.shape = 256, 256

    #imshow(im, ColormapJet(256))
    #imshow(im, cmap=cm.jet)
    
    imvals = np.sort(im.flatten())
    lo = imvals[0]
    hi = imvals[-1]
    steps = (imvals[::len(imvals)/256] - lo) / (hi - lo)
    num_steps = float(len(steps))
    interps = [(s, idx/num_steps, idx/num_steps) for idx, s in enumerate(steps)]
    interps.append((1, 1, 1))
    cdict = {'red' : interps,
             'green' : interps,
             'blue' : interps}
    histeq_cmap = colors.LinearSegmentedColormap('HistEq', cdict)
    pylab.figure()
    pylab.imshow(im, cmap=histeq_cmap)
    pylab.title('histeq')
    pylab.show()