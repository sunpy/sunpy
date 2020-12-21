"""
This module provides a set of colormaps specific for solar data.
"""
from copy import deepcopy

import matplotlib.cm as mplcm
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u

from sunpy.visualization.colormaps import color_tables as ct

__all__ = ['show_colormaps', 'cmlist']

sdoaia94 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[94*u.angstrom], 'SDO AIA {:s}'.format(str(94*u.angstrom)))
sdoaia131 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[131*u.angstrom], 'SDO AIA {:s}'.format(str(131*u.angstrom)))
sdoaia171 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[171*u.angstrom], 'SDO AIA {:s}'.format(str(171*u.angstrom)))
sdoaia193 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[193*u.angstrom], 'SDO AIA {:s}'.format(str(193*u.angstrom)))
sdoaia211 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[211*u.angstrom], 'SDO AIA {:s}'.format(str(211*u.angstrom)))
sdoaia304 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[304*u.angstrom], 'SDO AIA {:s}'.format(str(304*u.angstrom)))
sdoaia335 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[335*u.angstrom], 'SDO AIA {:s}'.format(str(335*u.angstrom)))
sdoaia1600 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[1600*u.angstrom], 'SDO AIA {:s}'.format(str(1600*u.angstrom)))
sdoaia1700 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[1700*u.angstrom], 'SDO AIA {:s}'.format(str(1700*u.angstrom)))
sdoaia4500 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[4500*u.angstrom], 'SDO AIA {:s}'.format(str(4500*u.angstrom)))

sohoeit171 =  ct.cmap_from_rgb_file('SOHO EIT {:s}'.format(str(171*u.angstrom)), f'eit_{ct.Colors.color.value[171*u.angstrom]}.csv')
sohoeit195 =  ct.cmap_from_rgb_file('SOHO EIT {:s}'.format(str(195*u.angstrom)), f'eit_{ct.Colors.color.value[195*u.angstrom]}.csv')
sohoeit284 =  ct.cmap_from_rgb_file('SOHO EIT {:s}'.format(str(284*u.angstrom)), f'eit_{ct.Colors.color.value[284*u.angstrom]}.csv')
sohoeit304 =  ct.cmap_from_rgb_file('SOHO EIT {:s}'.format(str(304*u.angstrom)), f'eit_{ct.Colors.color.value[304*u.angstrom]}.csv')
    
goesrsuvi94 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[94*u.angstrom], 'GOES-R SUVI {:s}'.format(str(94*u.angstrom)))
goesrsuvi131 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[131*u.angstrom], 'GOES-R SUVI {:s}'.format(str(131*u.angstrom)))
goesrsuvi171 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[171*u.angstrom], 'GOES-R SUVI {:s}'.format(str(171*u.angstrom)))
goesrsuvi304 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[304*u.angstrom], 'GOES-R SUVI {:s}'.format(str(304*u.angstrom)))

goesrsuvi195 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[193*u.angstrom], 'GOES-R SUVI {:s}'.format(str(193*u.angstrom)))
goesrsuvi284 = ct._cmap_from_rgb(ct.Colors.aia_wave_dict.value[335*u.angstrom], 'GOES-R SUVI {:s}'.format(str(335*u.angstrom)))

# The color tables below returns one of the fundamental color tables for SOHO
# LASCO images. These are not the same as those used in SSWIDL.  This is
# because the SSWIDL color scaling for LASCO level 0.5 and 1.0 is highly
# compressed and does not display the data well.

soholasco2 = deepcopy(mplcm.get_cmap("gist_heat"))
soholasco2.name = 'SOHO LASCO C2'
soholasco3 = deepcopy(mplcm.get_cmap("bone"))
soholasco3.name = 'SOHO LASCO C3'

# These are the SSWIDL color tables.
sswidlsoholasco2 = ct.cmap_from_rgb_file(f'SOHO LASCO C{2}', f'lasco_c{2}.csv')
sswidlsoholasco3 = ct.cmap_from_rgb_file(f'SOHO LASCO C{3}', f'lasco_c{3}.csv')

stereocor1 = ct.cmap_from_rgb_file(f'STEREO COR{1}', f'stereo_cor{1}.csv')
stereocor2 = ct.cmap_from_rgb_file(f'STEREO COR{2}', f'stereo_cor{2}.csv')

stereohi1 = ct.cmap_from_rgb_file(f'STEREO HI{1}', f'hi{1}.csv')
stereohi2 = ct.cmap_from_rgb_file(f'STEREO HI{2}', f'hi{2}.csv')

yohkohsxtal = ct._cmap_from_rgb(ct.Colors.sxt_colors.value['al'], 'Yohkoh SXT {:s}'.format('al'.title()))
yohkohsxtwh = ct._cmap_from_rgb(ct.Colors.sxt_colors.value['wh'], 'Yohkoh SXT {:s}'.format('wh'.title()))

hinodexrt = ct._cmap_from_rgb(ct.Colors.rgb_colors.value['intensity'], 'Hinode XRT')
hinodesotintensity = ct._cmap_from_rgb(ct.Colors.rgb_colors.value['intensity'], 'Hinode SOT intensity')


trace171 = ct.cmap_from_rgb_file('TRACE 171', 'trace_171.csv')
trace195 = ct.cmap_from_rgb_file('TRACE 195', 'trace_195.csv')
trace284 = ct.cmap_from_rgb_file('TRACE 284', 'trace_284.csv')
trace1216 = ct.cmap_from_rgb_file('TRACE 1216', 'trace_1216.csv')
trace1550 = ct.cmap_from_rgb_file('TRACE 1550', 'trace_1550.csv')
trace1600 = ct.cmap_from_rgb_file('TRACE 1600', 'trace_1600.csv')
trace1700 = ct.cmap_from_rgb_file('TRACE 1700', 'trace_1700.csv')

traceWL = ct.cmap_from_rgb_file('TRACE WL', 'grayscale.csv')

hmimag = ct.cmap_from_rgb_file('SDO HMI magnetogram', 'hmi_mag.csv')

kcor = deepcopy(mplcm.get_cmap("gist_gray"))
kcor.name = 'MLSO KCor'

rhessi = ct.cmap_from_rgb_file("rhessi", "rhessi.csv")
std_gamma_2 = ct.cmap_from_rgb_file("std_gamma_2", "std_gamma_2.csv")

cmlist = {
    'goes-rsuvi94': goesrsuvi94,
    'goes-rsuvi131': goesrsuvi131,
    'goes-rsuvi171': goesrsuvi171,
    'goes-rsuvi195': goesrsuvi195,
    'goes-rsuvi284': goesrsuvi284,
    'goes-rsuvi304': goesrsuvi304,
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
    'sswidlsoholasco2': sswidlsoholasco2,
    'sswidlsoholasco3': sswidlsoholasco3,
    'stereocor1': stereocor1,
    'stereocor2': stereocor2,
    'stereohi1': stereohi1,
    'stereohi2': stereohi2,
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
    'irissjiSJI_NUV': ct.iris_sji_color_table('SJI_NUV'),
    'kcor': kcor,
    'rhessi': rhessi,
    'std_gamma_2': std_gamma_2,
}

# Register the colormaps with matplotlib so plt.get_cmap('sdoaia171') works
for name, cmap in cmlist.items():
    mplcm.register_cmap(name=name, cmap=cmap)


def show_colormaps(search=None):
    """
    Displays a plot of the custom color maps supported in SunPy.
    Parameters
    ----------
    search : str
        A string to search for in the names of the color maps (e.g. aia, EIT,
        171). Case insensitive.
    Examples
    --------
    >>> import sunpy.visualization.colormaps as cm
    >>> cm.show_colormaps()  # doctest: +IGNORE_WARNINGS
    >>> cm.show_colormaps(search='aia')  # doctest: +IGNORE_WARNINGS
    >>> cm.show_colormaps(search='171')  # doctest: +IGNORE_WARNINGS
    """

    if search is not None:
        maps = sorted({k: v for (k, v) in cmlist.items() if k.lower().count(search.lower())})
        if len(maps) == 0:
            raise KeyError(f'No color maps found for search term "{search:s}"')
    else:
        maps = sorted(cmlist)

    nmaps = len(maps) + 1

    a = np.linspace(0, 1, 256).reshape(1, -1)
    a = np.vstack((a, a))

    fig = plt.figure(figsize=(7, 10), dpi=128)
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    for i, name in enumerate(maps):
        ax = plt.subplot(nmaps, 1, i + 1)
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=name, origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], name, fontsize=10, horizontalalignment='right')
    plt.show()