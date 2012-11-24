"""
This module provides a set of colormaps specific for solar data.
"""
from __future__ import absolute_import

__all__ = ["get_cmap", "show_colormaps"]

import numpy as np
import matplotlib.pyplot as plt
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

sohoeit171 = _cm.eit_color_table(171)
sohoeit195 = _cm.eit_color_table(195)
sohoeit284 = _cm.eit_color_table(284)
sohoeit304 = _cm.eit_color_table(304)

soholasco2 = _cm.lasco_color_table(2)
soholasco3 = _cm.lasco_color_table(3)

yohkohsxtal = _cm.sxt_color_table('al')
yohkohsxtwh = _cm.sxt_color_table('wh')

hinodexrt = _cm.xrt_color_table()



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
          'sohoeit171': sohoeit171,
          'sohoeit195': sohoeit195,
          'sohoeit284': sohoeit284,
          'sohoeit304': sohoeit304,
          'soholasco2': soholasco2,
          'soholasco3': soholasco3,
          'rhessi': cm.jet,  # pylint: disable=E1101
          'yohkohsxtal': yohkohsxtal,
          'yohkohsxtwh': yohkohsxtwh,
          'hinodexrt': hinodexrt
          }


def get_cmap(name='sdoaia94'):
    """Get a colormap.

    Parameters
    ----------
    name : string
        The name of a color map.

    Returns
    -------
    value : matplotlib colormap

    See Also
    --------

    Examples
    --------
    >>> import sunpy.cm as cm
    >>> colormap = cm.get_cmap(name = 'sdoaia94')

    References
    ----------
    | http://matplotlib.sourceforge.net/api/cm_api.html

    """
    if name in cmlist:
        return cmlist.get(name)
    else:
        raise ValueError("Colormap %s is not recognized" % name)


def show_colormaps():
    """Displays a plot of the custom color maps supported in SunPy.

    Parameters
    ----------
    None : none

    Returns
    -------
    None : none

    See Also
    --------

    Examples
    --------
    >>> import sunpy.cm as cm
    >>> cm.show_colormaps()

    References
    ----------

    """
    maps = sorted(cmlist)
    nmaps = len(maps) + 1

    a = np.linspace(0, 1, 256).reshape(1, -1)  # pylint: disable=E1103
    a = np.vstack((a, a))

    fig = plt.figure(figsize=(5, 10))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    for i, name in enumerate(maps):
        ax = plt.subplot(nmaps, 1, i + 1)
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=get_cmap(name), origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], name, fontsize=10,
                 horizontalalignment='right')

    plt.show()

#def test_equalize(data):
#    """Returns a color map which performs histogram equalization on the data.
#
#    Parameters
#    ----------
#    data : ndarray
#
#    Returns
#    -------
#    value : matplotlib colormap
#
#    See Also
#    --------
#
#    Examples
#    --------
#    >>> import sunpy.cm as cm
#    >>> cm.test_equalize()
#
#    Reference
#    ---------
#    | http://matplotlib.sourceforge.net/api/cm_api.html
#
#    .. warning:: this function is under development
#
#    .. todo:: finish coding this function!
#
#    """
#    dfile = cbook.get_sample_data('s1045.ima', asfileobj=False)
#
#    im = np.fromstring(file(dfile, 'rb').read(), np.uint16).astype(float)
#    im.shape = 256, 256
#
#    #imshow(im, ColormapJet(256))
#    #imshow(im, cmap=cm.jet)
#
#    imvals = np.sort(im.flatten())
#    lo = imvals[0]
#    hi = imvals[-1]
#    steps = (imvals[::len(imvals)/256] - lo) / (hi - lo)
#    num_steps = float(len(steps))
#    interps = [(s, idx/num_steps, idx/num_steps) for idx,
#        s in enumerate(steps)]
#    interps.append((1, 1, 1))
#    cdict = {'red': interps,
#             'green': interps,
#             'blue': interps}
#    histeq_cmap = colors.LinearSegmentedColormap('HistEq', cdict)
#    pylab.figure()
#    pylab.imshow(im, cmap=histeq_cmap)
#    pylab.title('histeq')
#    pylab.show()
