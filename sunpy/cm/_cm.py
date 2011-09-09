from __future__ import absolute_import

"""
Nothing here but dictionaries for generating LinearSegmentedColormaps,
and a dictionary of these dictionaries.

"""

import numpy as np
import matplotlib.colors as colors

# FIXME: Give me a proper name.
def _mkx(i, steps, n):
    """ Generate list according to pattern of g0 and b0. """
    x = []
    for step in steps:
        x.extend(range(i, step + n, n))
        i = step + (n - 1)
    return x

def padfr(lst, len_, pad=0):
    """ Pad lst to contain at least len_ items by adding pad to the front. """
    diff = len_ - len(lst)
    diff = 0 if diff < 0 else diff
    return [pad] * diff + lst

def paden(lst, len_, pad=0):
    """ Pad lst to contain at least len_ items by adding pad to the end. """
    diff = len_ - len(lst)
    diff = 0 if diff < 0 else diff
    return lst + [pad] * diff


# The following values describe color table 3 for IDL (Red Temperature)
r0 = np.array(paden([0,1,2,4,5,7,8,10,11,13,14,15,17,18,20,21,23,24,26,27,28,30,31,33,34,36,37,39,40,42,43,44,46,47,49,50,52,53,55,56,57,59,60,62,63,65,66,68,69,70,72,73,75,76,78,79,81,82,84,85,86,88,89,91,92,94,95,97,98,99,101,102,104,105,107,108,110,111,113,114,115,117,118,120,121,123,124,126,127,128,130,131,133,134,136,137,139,140,141,143,144,146,147,149,150,152,153,155,156,157,159,160,162,163,165,166,168,169,170,172,173,175,176,178,179,181,182,184,185,186,188,189,191,192,194,195,197,198,199,201,202,204,205,207,208,210,211,212,214,215,217,218,220,221,223,224,226,227,228,230,231,233,234,236,237,239,240,241,243,244,246,247,249,250,252,253], 256, 255))
g0 = np.array(padfr(_mkx(1, xrange(17, 256, 17), 2), 256))
b0 = np.array(padfr(_mkx(3, xrange(51, 256, 51), 4), 256))
    
c0 = np.arange(256, dtype = 'f')
c1 = (np.sqrt(c0)*np.sqrt(255.0)).astype('f')
c2 = (np.arange(256)**2/255.0).astype('f')
c3 = ((c1+c2/2.0)*255.0/(c1.max() + c2.max()/2.0)).astype('f')

def aia_color_table(wavelength):
    '''Returns one of the fundamental color tables for SDO AIA images.
       Based on aia_lct.pro part of SDO/AIA on SSWIDL written by Karl Schriver (2010/04/12).
    '''
    try:
        r, g, b = {
            1600: (c3, c3, c2), 1700: (c1, c0, c0), 4500: (c0, c0, b0 / 2.0),
            94: (c2, c3, c0), 131: (g0, r0, r0), 171: (r0, c0, b0),
            193: (c1, c0, c2), 211: (c1, c0, c3), 304: (r0, g0, b0),
            335: (c2, c0, c1)
        }[wavelength]
    except KeyError:
        raise ValueError(
            "Invalid AIA wavelength. Valid values are "
            "1600,1700,4500,94,131,171,193,211,304,335."
        )
   
    # Now create the color tuples
    i = np.linspace(0, 1, r0.size)
    
    cdict = dict(
        (name, list(zip(i, el/255.0, el/255.0)))
        for el, name in [(r, 'red'),  (g, 'green'), (b, 'blue')]
    )
    
    return colors.LinearSegmentedColormap('mytable', cdict)
