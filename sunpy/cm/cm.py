"""
This module provides a set of colormaps specific for solar data.
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from sunpy.cm import color_tables as ct
from sunpy.util import deprecated

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
          'rhessi': cm.jet,
          'yohkohsxtal': yohkohsxtal,
          'yohkohsxtwh': yohkohsxtwh,
          'hinodexrt': hinodexrt,
          'hinodesotintensity': hinodesotintensity,
          'trace171': trace171,
          'trace195': trace195,
          'trace284': trace284,
          'trace1216': trace1216,
          'trace1550': trace1550,
          'trace1600': trace1600,
          'trace1700': trace1700,
          'traceWL': traceWL,
          'hmimag': hmimag,
          'irissji1330': ct.iris_sji_color_table('1330'),
          'irissji1400': ct.iris_sji_color_table('1400'),
          'irissji1600': ct.iris_sji_color_table('1600'),
          'irissji2796': ct.iris_sji_color_table('2796'),
          'irissji2832': ct.iris_sji_color_table('2832'),
          'irissji5000': ct.iris_sji_color_table('5000'),
          'irissjiFUV': ct.iris_sji_color_table('FUV'),
          'irissjiNUV': ct.iris_sji_color_table('NUV'),
          'irissjiSJI_NUV': ct.iris_sji_color_table('SJI_NUV')
}

# Register the colormaps with matplotlib so plt.get_cmap('sdoaia171') works
for name, cmap in cmlist.items():
    cm.register_cmap(name=name, cmap=cmap)

@deprecated("0.9", "Use Matplotlib to load the colormaps", alternative='plt.get_cmap')
def get_cmap(name):
    """
    Get a colormap.

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
    | https://matplotlib.org/api/cm_api.html

    """
    if name in cmlist:
        return cmlist.get(name)
    else:
        raise ValueError("Colormap {name!s} is not recognized".format(name=name))


def show_colormaps(filter=None):
    """Displays a plot of the custom color maps supported in SunPy.

    Parameters
    ----------
    filter : str
        A string to filter the color maps presented (e.g. aia, EIT, 171). Case
        insensitive.

    Returns
    -------
    None : none

    Examples
    --------
    >>> import sunpy.cm as cm
    >>> cm.show_colormaps()
    >>> cm.show_colormaps(filter='aia')
    >>> cm.show_colormaps(filter='171')

    References
    ----------

    """

    if filter:
        maps =  sorted({k:v for (k,v) in cmlist.items() if k.lower().count(filter.lower())})
        if len(maps) == 0:
            raise KeyError('No color maps found for key - ' + filter)
    else:
        maps = sorted(cmlist)

    nmaps = len(maps) + 1

    a = np.linspace(0, 1, 256).reshape(1, -1)  # pylint: disable=E1103
    a = np.vstack((a, a))

    fig = plt.figure(figsize=(5, 10),dpi=64)
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    for i, name in enumerate(maps):
        ax = plt.subplot(nmaps, 1, i + 1)
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=get_cmap(name), origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], name, fontsize=10,
                 horizontalalignment='right')
    plt.show()

