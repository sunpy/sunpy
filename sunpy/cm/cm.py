from __future__ import absolute_import
"""
This module provides a set of colormaps specific to solar data (e.g. SDO/AIA 
color maps), functions for getting a colormap by name.
"""

import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib.cm as cm
from sunpy.cm import _cm

sdoaia94 = _cm.aia_color_table(94)
sdoaia131 = _cm.aia_color_table(131)
sdoaia171 = _cm.aia_color_table(171)
sdoaia193 = _cm.aia_color_table(193)
sdoaia211 = _cm.aia_color_table(211)
sdoaia304 = _cm.aia_color_table(304)
sdoaia335 = _cm.aia_color_table(335)
sdoaia1600 = _cm.aia_color_table(1600)
sdoaia1700 = _cm.aia_color_table(1700)
sdoaia4500 = _cm.aia_color_table(4500)

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
          'rhessi': cm.jet #pylint: disable=E1101
          }

def get_cmap(name='sdoaia94'):
    """Get a colormap instance."""
    if name in cmlist:
        return cmlist.get(name)
    else:
        raise ValueError("Colormap %s is not recognized" % name)

def show_colormaps():
    """Displays custom color maps supported in SunPy"""
    maps = sorted(cmlist)
    nmaps = len(maps) + 1
    
    a = np.linspace(0, 1, 256).reshape(1, -1) #pylint: disable=E1103
    a = np.vstack((a, a))
    
    fig = plt.figure(figsize=(5, 10))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    for i,name in enumerate(maps):
        ax = plt.subplot(nmaps, 1, i + 1)
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=get_cmap(name), origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], name, fontsize=10, 
                 horizontalalignment='right')

    plt.show()

def test_equalize():
    '''Test'''
    import pylab

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
    cdict = {'red': interps,
             'green': interps,
             'blue': interps}
    histeq_cmap = colors.LinearSegmentedColormap('HistEq', cdict)
    pylab.figure()
    pylab.imshow(im, cmap=histeq_cmap)
    pylab.title('histeq')
    pylab.show()

if __name__ == '__main__':
    test_equalize()
