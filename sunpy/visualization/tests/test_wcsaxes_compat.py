import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import WCS

from sunpy.visualization import wcsaxes_compat


@pytest.fixture
def dummy_wcs():
    return WCS({
        'ctype1': 'HPLN-TAN',
        'ctype2': 'HPLT-TAN',
        'cunit1': 'arcsec',
        'cunit2': 'arcsec',
        'cdelt1': 1,
        'cdelt2': 1,
        'crpix1': 50,
        'crpix2': 50,
        'crval1': 0,
        'crval2': 0,
        'mjd-obs': 58849.0,
        'date-obs': '2020-01-01T00:00:00',
    })


def test_is_wcsaxes(dummy_wcs):
    fig = Figure()
    normal_ax = fig.add_subplot(121)
    assert not wcsaxes_compat.is_wcsaxes(normal_ax)

    wcs_ax = fig.add_subplot(122, projection=dummy_wcs)
    assert wcsaxes_compat.is_wcsaxes(wcs_ax)


def test_gca_wcs(dummy_wcs):
    plt.close('all')
    try:
        # Test when figure has no axes: creates a new WCSAxes
        fig = plt.figure()
        assert len(fig.get_axes()) == 0
        ax = wcsaxes_compat.gca_wcs(dummy_wcs, fig=fig)
        assert wcsaxes_compat.is_wcsaxes(ax)
        assert len(fig.get_axes()) == 1

        # Test when figure already has axes: returns current axes
        ax2 = wcsaxes_compat.gca_wcs(dummy_wcs, fig=fig)
        assert ax2 is ax

        # Test with default fig (None)
        plt.close('all')
        fig = plt.figure()
        ax_default = wcsaxes_compat.gca_wcs(dummy_wcs)
        assert ax_default.figure is fig
        assert wcsaxes_compat.is_wcsaxes(ax_default)

        # Test with slices passed
        plt.close('all')
        fig = plt.figure()
        ax_sliced = wcsaxes_compat.gca_wcs(dummy_wcs, fig=fig, slices=('x', 'y'))
        assert wcsaxes_compat.is_wcsaxes(ax_sliced)
    finally:
        plt.close('all')


def test_get_world_transform(dummy_wcs):
    fig = Figure()
    normal_ax = fig.add_subplot(121)
    trans = wcsaxes_compat.get_world_transform(normal_ax)
    assert trans is normal_ax.transData

    wcs_ax = fig.add_subplot(122, projection=dummy_wcs)
    world_trans = wcsaxes_compat.get_world_transform(wcs_ax)
    assert world_trans == wcs_ax.get_transform('world')


def test_default_wcs_grid(dummy_wcs):
    fig = Figure()
    wcs_ax = fig.add_subplot(111, projection=dummy_wcs)
    wcsaxes_compat.default_wcs_grid(wcs_ax)
    assert hasattr(wcs_ax.coords, 'grid')


def test_wcsaxes_heliographic_overlay_stonyhurst(dummy_wcs):
    fig = Figure()
    wcs_ax = fig.add_subplot(111, projection=dummy_wcs)
    overlay = wcsaxes_compat.wcsaxes_heliographic_overlay(
        wcs_ax, grid_spacing=15 * u.deg, system='stonyhurst'
    )
    assert overlay is not None

    # Check axis labels
    assert overlay[0].get_axislabel() == 'Stonyhurst Longitude'
    assert overlay[1].get_axislabel() == 'Stonyhurst Latitude'

    # Check native coords tick positions set to 'bl'
    for c in wcs_ax.coords:
        assert set(c.get_ticks_position()) == {'b', 'l'}


def test_wcsaxes_heliographic_overlay_spacing_tuple(dummy_wcs):
    fig = Figure()
    wcs_ax = fig.add_subplot(111, projection=dummy_wcs)
    overlay = wcsaxes_compat.wcsaxes_heliographic_overlay(
        wcs_ax, grid_spacing=[10, 20] * u.deg, system='stonyhurst'
    )
    assert overlay is not None


def test_wcsaxes_heliographic_overlay_carrington(dummy_wcs):
    fig = Figure()
    wcs_ax = fig.add_subplot(111, projection=dummy_wcs)
    observer = SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.AU, frame='heliographic_stonyhurst', obstime='2020-01-01')
    overlay = wcsaxes_compat.wcsaxes_heliographic_overlay(
        wcs_ax, system='carrington', observer=observer, obstime=Time('2020-01-01')
    )
    assert overlay[0].get_axislabel() == 'Carrington Longitude'
    assert overlay[1].get_axislabel() == 'Carrington Latitude'


def test_wcsaxes_heliographic_overlay_annotate_false(dummy_wcs):
    fig = Figure()
    wcs_ax = fig.add_subplot(111, projection=dummy_wcs)
    overlay = wcsaxes_compat.wcsaxes_heliographic_overlay(wcs_ax, annotate=False)
    assert not overlay[0].get_ticks_visible()
    assert not overlay[1].get_ticks_visible()
    assert not overlay[0].get_ticklabel_visible()
    assert not overlay[1].get_ticklabel_visible()


def test_wcsaxes_heliographic_overlay_title_offset(dummy_wcs):
    fig = Figure()
    wcs_ax = fig.add_subplot(111, projection=dummy_wcs)
    wcs_ax.set_title('Test Title')
    orig_x, orig_y = wcs_ax.title.get_position()
    wcsaxes_compat.wcsaxes_heliographic_overlay(wcs_ax)
    new_x, new_y = wcs_ax.title.get_position()
    assert new_x == orig_x
    assert pytest.approx(new_y) == orig_y + 0.08


def test_wcsaxes_heliographic_overlay_errors(dummy_wcs):
    fig = Figure()
    wcs_ax = fig.add_subplot(111, projection=dummy_wcs)

    # Invalid grid spacing length (3 items)
    with pytest.raises(ValueError, match="grid_spacing must be a Quantity of length one or two"):
        wcsaxes_compat.wcsaxes_heliographic_overlay(wcs_ax, grid_spacing=[10, 20, 30] * u.deg)

    # Invalid system
    with pytest.raises(ValueError, match="system must be 'stonyhurst' or 'carrington'"):
        wcsaxes_compat.wcsaxes_heliographic_overlay(wcs_ax, system='heliocentric')
