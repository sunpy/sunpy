import matplotlib
import pytest

import astropy.units as u

import sunpy.visualization.colormaps as cm
import sunpy.visualization.colormaps.color_tables as ct
from sunpy.tests.helpers import figure_test


# Checks that colormaps are imported by MPL
def test_get_cmap():
    for cmap in cm.cmlist.keys():
        assert cm.cmlist[cmap] == matplotlib.colormaps[cmap]


def test_invalid_show_cmaps():
    with pytest.raises(
            KeyError, match='No color maps found for search term'):
        cm.show_colormaps(search='asdfghjkl')


@pytest.mark.parametrize(
    ('f', 'match'),
    [(ct.aia_color_table, 'Invalid AIA wavelength.'),
     (ct.eit_color_table, 'Invalid EIT wavelength.'),
     (ct.sswidl_lasco_color_table, 'Invalid LASCO number.'),
     (ct.sxt_color_table, 'Invalid SXT filter type number.'),
     (ct.cor_color_table, 'Invalid COR number.'),
     (ct.trace_color_table, 'Invalid TRACE filter waveband passed.'),
     (ct.sot_color_table, r'Invalid \(or not supported\) SOT type.'),
     (ct.iris_sji_color_table, 'Invalid IRIS SJI waveband.'),
     (ct.stereo_hi_color_table, 'Valid HI cameras are 1 and 2'),
     (ct.suvi_color_table, 'Invalid SUVI wavelength.')]
)
def test_invalid_wavelengths(f, match):
    # Check that accessing non-existent color table values raises a ValueError
    with pytest.raises(ValueError, match=match):
        f(-100*u.m)


@figure_test
def test_cmap_visual():
    cm.show_colormaps()
