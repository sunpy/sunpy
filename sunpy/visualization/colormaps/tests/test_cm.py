import matplotlib
import pytest
from numpy import linspace

import astropy.units as u

import sunpy.visualization.colormaps as cm
import sunpy.visualization.colormaps.color_tables as ct
from sunpy.tests.helpers import figure_test

SUIT_FILTERS = [
    "suit_nb01","suit_nb02","suit_nb03","suit_nb04",
    "suit_nb05","suit_nb06","suit_nb07","suit_nb08",
    "suit_bb01","suit_bb02","suit_bb03"
]

@pytest.mark.parametrize("filter_name", SUIT_FILTERS)
def test_suit_colormap_callable(filter_name):
    """Check SUIT colormap exists in cmlist, callable, and returns 256 colors."""
    cmap = cm.cmlist[filter_name]
    assert callable(cmap)
    sample = cmap(linspace(0,1,256))[:,:3]
    assert sample.shape == (256,3)

# Separate check for SUIT invalid filters
@pytest.mark.parametrize("invalid_filter", ["not_a_filter", "", "NB99"])
def test_invalid_suit_filter(invalid_filter):
    """Check that invalid SUIT filter names raise ValueError."""
    import sunpy.visualization.colormaps.color_tables as ct
    with pytest.raises(ValueError, match=r"Invalid Band"):
        ct.suit_color_table(invalid_filter)


@pytest.mark.thread_unsafe(reason="bug fixed in matplotlib dev")
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
