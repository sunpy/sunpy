"""
This module provides a set of colormaps specific for solar data.
"""
from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from sunpy.cm import color_tables as ct

__all__ = ['get_cmap', 'show_colormaps', 'cmlist']

sdoaia94 = ct.aia_color_table(94)
sdoaia131 = ct.aia_color_table(131)
sdoaia171 = ct.aia_color_table(171)
sdoaia193 = ct.aia_color_table(193)
sdoaia211 = ct.aia_color_table(211)
sdoaia304 = ct.aia_color_table(304)
sdoaia335 = ct.aia_color_table(335)
sdoaia1600 = ct.aia_color_table(1600)
sdoaia1700 = ct.aia_color_table(1700)
sdoaia4500 = ct.aia_color_table(4500)

sohoeit171 = ct.eit_color_table(171)
sohoeit195 = ct.eit_color_table(195)
sohoeit284 = ct.eit_color_table(284)
sohoeit304 = ct.eit_color_table(304)

soholasco2 = ct.lasco_color_table(2)
soholasco3 = ct.lasco_color_table(3)

stereocor1 = ct.cor_color_table(1)
stereocor2 = ct.cor_color_table(2)

stereohi1 = ct.stereo_hi_color_table(1)
stereohi2 = ct.stereo_hi_color_table(2)

yohkohsxtal = ct.sxt_color_table('al')
yohkohsxtwh = ct.sxt_color_table('wh')

hinodexrt = ct.xrt_color_table()
hinodesotintensity = ct.sot_color_table('intensity')
#hinodesotstokesquv = ct.sot_color_table('stokesQUV')
#hinodesotmagneticf = ct.sot_color_table('magnetic field')
#hinodesotvelocity = ct.sot_color_table('velocity')
#hinodesotwidth =  ct.sot_color_table('width')

trace171 = ct.trace_color_table('171')
trace195 = ct.trace_color_table('195')
trace284 = ct.trace_color_table('284')
trace1216 = ct.trace_color_table('1216')
trace1550 = ct.trace_color_table('1550')
trace1600 = ct.trace_color_table('1600')
trace1700 = ct.trace_color_table('1700')
traceWL = ct.trace_color_table('WL')

hmimag = ct.hmi_mag_color_table()

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
          'stereocor1': stereocor1,
          'stereocor2': stereocor2,
          'stereohi1': stereohi1,
          'stereohi2': stereohi2,
          'rhessi': cm.jet,  # pylint: disable=E1101
          'yohkohsxtal': yohkohsxtal,
          'yohkohsxtwh': yohkohsxtwh,
          'hinodexrt': hinodexrt,
          'hinodesotintensity': hinodesotintensity,
          #'hinodesotstokesquv': hinodesotstokesquv,
          #'hinodesotmagneticf': hinodesotmagneticf,
          #'hinodesotvelocity': hinodesotvelocity,
          #'hinodesotwidth': hinodesotwidth,
          'trace171': trace171,
          'trace195': trace195,
          'trace284': trace284,
          'trace1216': trace1216,
          'trace1550': trace1550,
          'trace1600': trace1600,
          'trace1700': trace1700,
          'traceWL': traceWL,
          'hmimag': hmimag
          }

# Register the colormaps with matplotlib so plt.get_cmap('sdoaia171') works
for name, cmap in cmlist.items():
    cm.register_cmap(name=name, cmap=cmap)

def get_cmap(name):
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
        raise ValueError("Colormap {name!s} is not recognized".format(name=name))


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

# def test_equalize(data):
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
