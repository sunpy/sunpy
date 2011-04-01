#-*- coding:utf-8 -*-
#
# Author: Keith Hughitt <keith.hughitt@nasa.gov>
# Author: Jack Ireland <jack.ireland@nasa.gov>
#
# Written: 2011/04/01
#
# <License info will go here...>
#
import matplotlib.colors
import operator
import numpy

#
# Things to try:
#   use sub-sampling to speed up histogram generation
#
def log_adaptive(data, bits=12):
    """Creates a custom color map scaled to the specified image data."""
    
    # Number of entries in the color map
    cmap_size = 2 ** bits
    
    # Ignore negative values and outliers
    data = numpy.log(data.clip(1, 2000))
    
    # Sort bins by frequency
    sorted_counts = _get_pixel_counts(data, cmap_size)
    indices = sorted([x[0] for x in sorted_counts[0:cmap_size - 1]])
    
    # Take exponent of each indice
    indices = [numpy.exp(x) for x in indices]
    
    # Scale from 0 to 1
    maximum = numpy.exp(data.max())
    indices =  [x / maximum for x in indices]
    
    # Create a matplotlib-formatted color dictionary
    cdict = _generate_cdict_for_indices(indices, cmap_size)
    
    return matplotlib.colors.LinearSegmentedColormap('automap', cdict, cmap_size)


def _generate_cdict_for_indices(indices, cmap_size):
    """Converts a list of indice values to an RGB color dictionary needed 
       to generate a linear segmented colormap
       
       See: http://matplotlib.sourceforge.net/api/colors_api.html#matplotlib.colors.LinearSegmentedColormap
    """
    step = 1. / cmap_size
    cdict = {'red': [], 'green': [], 'blue': []}
    
    value = 0
    
    for i in indices:
        t = (i, value, value)
        cdict['red'].append(t)
        cdict['green'].append(t)
        cdict['blue'].append(t)
        value += step
        
    # cmap values must range from 0 to 1
    cdict['red'][0] = cdict['green'][0] = cdict['blue'][0] = (0, 0, 0)
    cdict['red'][-1] = cdict['green'][-1] = cdict['blue'][-1] = (1, 1, 1)
    
    # convert rgb lists to tuples
    cdict['red'] = tuple(cdict['red'])
    cdict['green'] = tuple(cdict['green'])
    cdict['blue'] = tuple(cdict['blue'])

    return cdict

def _get_pixel_counts(data, cmap_size):
    """Returns a list of the pixel values sorted by frequency"""
    # Bin the data
    hist, bins = numpy.histogram(data, bins=cmap_size)
    
    a = {}
    i = 0
    
    # Create a dictionary which maps from mean bin values to counts
    for freq in hist:
        a[(bins[i] + bins[i + 1]) / 2] = freq
        i += 1
    
    # Return as a list of tuples sorted by frequency
    return sorted(a.iteritems(), key=operator.itemgetter(1), reverse=True)

