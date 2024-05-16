"""
Test Generic Map
"""
import copy

import matplotlib.colors as mcolor
import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.figure import Figure

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy
import sunpy.coordinates
import sunpy.map
from sunpy.coordinates import HeliographicStonyhurst
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test, fix_map_wcs
from sunpy.util.exceptions import SunpyUserWarning

pytestmark = pytest.mark.filterwarnings('ignore:Missing metadata')


@pytest.fixture
def heliographic_test_map():
    m = sunpy.map.Map(get_test_filepath('heliographic_phase_map.fits.gz'))
    return fix_map_wcs(m)


@pytest.fixture
def aia171_test_map_with_mask(aia171_test_map):
    shape = aia171_test_map.data.shape
    mask = np.zeros_like(aia171_test_map.data, dtype=bool)
    mask[0:shape[0] // 2, 0:shape[1] // 2] = True
    return sunpy.map.Map(np.ma.array(
        aia171_test_map.data, mask=mask),
        aia171_test_map.meta)


@figure_test
def test_plot_aia171(aia171_test_map):
    aia171_test_map.plot()


@figure_test
def test_plot_rotated_aia171(aia171_test_map):
    # Check that plotting a rotated map and a rectangle works as expected

    # Set rotation metadata
    aia171_test_map.meta['CROTA2'] = 45
    # Plot map
    aia171_test_map.plot()
    # Plot rectangle
    bottom_left = SkyCoord(
        0 * u.arcsec, 0 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    w = 100 * u.arcsec
    h = 200 * u.arcsec
    aia171_test_map.draw_quadrangle(bottom_left, width=w, height=h)


@figure_test
def test_plot_aia171_clip(aia171_test_map):
    aia171_test_map.plot(clip_interval=(5., 99.)*u.percent)


@figure_test
def test_peek_aia171(aia171_test_map):
    aia171_test_map.peek()


@figure_test
def test_peek_grid_aia171(aia171_test_map):
    aia171_test_map.peek(draw_grid=True)


@figure_test
def test_peek_grid_spacing_aia171(aia171_test_map):
    aia171_test_map.peek(draw_grid=(5, 5) * u.deg)


@figure_test
def test_peek_limb_aia171(aia171_test_map):
    aia171_test_map.peek(draw_limb=True)


@figure_test
def test_draw_grid_aia171(aia171_test_map):
    aia171_test_map.plot()
    aia171_test_map.draw_grid(grid_spacing=(30, 40) * u.deg)
    aia171_test_map.draw_grid(system='carrington', color='blue')


@figure_test
def test_peek_grid_limb_aia171(aia171_test_map):
    aia171_test_map.peek(draw_grid=True, draw_limb=True)


@figure_test
def test_rectangle_aia171_width_height(aia171_test_map):
    aia171_test_map.plot()
    bottom_left = SkyCoord(
        0 * u.arcsec, 0 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    w = 100 * u.arcsec
    h = 100 * u.arcsec
    aia171_test_map.draw_quadrangle(bottom_left, width=w, height=h)


@figure_test
def test_rectangle_aia171_top_right(aia171_test_map):
    aia171_test_map.plot()
    bottom_left = SkyCoord(
        0 * u.arcsec, 0 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    top_right = SkyCoord(
        100 * u.arcsec, 100 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    aia171_test_map.draw_quadrangle(bottom_left, top_right=top_right, label='Rectangle')
    plt.legend()  # Check that the 'Rectangle' label shows up in the legend


@figure_test
def test_quadrangle_aia17_width_height(aia171_test_map):
    aia171_test_map.plot()
    bottom_left = SkyCoord(
        50 * u.deg, -10 * u.deg, frame=HeliographicStonyhurst, obstime=aia171_test_map.date)
    w = 30 * u.deg
    h = 90 * u.deg
    aia171_test_map.draw_quadrangle(bottom_left=bottom_left.frame, width=w, height=h)


@figure_test
def test_quadrangle_aia17_pix_width_height(aia171_test_map):
    aia171_test_map.plot()
    aia171_test_map.draw_quadrangle(bottom_left=(50, 50)*u.pix, width=30*u.pix,
                                    height=50*u.pix, edgecolor="cyan")


@figure_test
def test_quadrangle_aia17_top_right(aia171_test_map):
    aia171_test_map.plot()
    bottom_left = SkyCoord(
        50 * u.deg, -10 * u.deg, frame=HeliographicStonyhurst, obstime=aia171_test_map.date)
    top_right = SkyCoord(
        65 * u.deg, 50 * u.deg, frame=HeliographicStonyhurst, obstime=aia171_test_map.date)
    aia171_test_map.draw_quadrangle(bottom_left, top_right=top_right, label='Quadrangle')
    plt.legend()  # Check that the 'Quadrangle' label shows up in the legend


@figure_test
def test_quadrangle_aia17_pix_top_right(aia171_test_map):
    aia171_test_map.plot()
    aia171_test_map.draw_quadrangle(bottom_left=(50, 50)*u.pix,
                                    top_right=(80, 90)*u.pix, edgecolor='cyan')


@figure_test
def test_quadrangle_aia17_pix_top_right_different_axes(aia171_test_map):
    # Plot the map rotated by 30 degrees
    aia171_test_map.rotate(30*u.deg).plot()
    # Plot a rectangle in the pixel space of the original map
    aia171_test_map.draw_quadrangle(bottom_left=(50, 50)*u.pix,
                                    top_right=(80, 90)*u.pix, edgecolor='cyan')


@figure_test
def test_plot_masked_aia171(aia171_test_map_with_mask):
    aia171_test_map_with_mask.plot()


@figure_test
def test_plot_aia171_superpixel(aia171_test_map):
    aia171_test_map.superpixel((3, 2) * u.pix, offset=(4, 4) * u.pix).plot()


@figure_test
def test_plot_resample(carrington_map):
    # Test that super-pixelling a map preserves the coordinate system correctly.
    # The two plots should have identical coordinate grids
    resamp = carrington_map.resample([10, 5] * u.pix)

    plt.figure()
    ax1 = plt.subplot(121, projection=carrington_map)
    ax2 = plt.subplot(122, projection=resamp)

    carrington_map.plot(axes=ax1)
    resamp.plot(axes=ax2)

    for ax in [ax1, ax2]:
        ax.coords.grid(True, color='tab:red', ls='solid', lw=2, alpha=1)


@figure_test
def test_plot_superpixel(carrington_map):
    # Test that super-pixelling a map preserves the coordinate system correctly.
    # The two plots should have identical coordinate grids
    superpix = carrington_map.superpixel([2, 2] * u.pix)

    plt.figure()
    ax1 = plt.subplot(121, projection=carrington_map)
    ax2 = plt.subplot(122, projection=superpix)

    carrington_map.plot(axes=ax1)
    superpix.plot(axes=ax2)

    for ax in [ax1, ax2]:
        ax.coords.grid(True, color='tab:red', ls='solid', lw=2, alpha=1)


@figure_test
def test_plot_masked_aia171_superpixel(aia171_test_map_with_mask):
    aia171_test_map_with_mask.superpixel(
        (9, 7) * u.pix, offset=(4, 4) * u.pix).plot()


@figure_test
def test_draw_contours_aia(aia171_test_map):
    aia171_test_map.plot()
    aia171_test_map.draw_contours(u.Quantity(np.arange(1, 100, 10), 'percent'))


@figure_test
def test_draw_contours_different_wcs(aia171_test_map):
    aia171_test_map._data = aia171_test_map.data.astype('float32')
    rotated_map = aia171_test_map.rotate(30*u.deg, order=3)
    rotated_map.plot()
    aia171_test_map.draw_contours(u.Quantity(np.arange(1, 100, 10), 'percent'))


@figure_test
def test_draw_contours_aia_fill(aia171_test_map):
    aia171_test_map.plot()
    aia171_test_map.draw_contours(u.Quantity(np.arange(1, 100, 10), 'percent'), fill=True)


@figure_test
def test_heliographic_peek(heliographic_test_map):
    heliographic_test_map.peek()


@figure_test
def test_heliographic_quadrangle_width_height(heliographic_test_map):
    heliographic_test_map.plot()
    bottom_left = SkyCoord(
        60 * u.deg, 50 * u.deg, frame=heliographic_test_map.coordinate_frame)
    w = 13 * u.deg
    h = 13 * u.deg
    heliographic_test_map.draw_quadrangle(bottom_left, width=w, height=h, edgecolor='cyan')


@figure_test
def test_heliographic_quadrangle_top_right(heliographic_test_map):
    heliographic_test_map.plot()
    bottom_left = SkyCoord(
        60 * u.deg, 50 * u.deg, frame=heliographic_test_map.coordinate_frame)
    top_right = SkyCoord(
        80 * u.deg, 90 * u.deg, frame=heliographic_test_map.coordinate_frame)
    heliographic_test_map.draw_quadrangle(bottom_left, top_right=top_right, edgecolor='cyan')


@figure_test
def test_heliographic_grid_annotations(heliographic_test_map):
    heliographic_test_map.plot()
    heliographic_test_map.draw_grid(annotate=False)


def test_plot_norm_error(aia171_test_map):
    # Check that duplicating vmin, vmax raises an error
    norm = mcolor.Normalize(vmin=0, vmax=1)
    with pytest.raises(ValueError, match='Cannot manually specify vmin'):
        aia171_test_map.plot(norm=norm, vmin=0)
    with pytest.raises(ValueError, match='Cannot manually specify vmax'):
        aia171_test_map.plot(norm=norm, vmax=0)


def test_quadrangle_no_wcsaxes(aia171_test_map):
    ax = Figure().add_subplot(projection=None)  # create a non-WCSAxes plot

    bottom_left = SkyCoord(
        [0, 1] * u.arcsec, [0, 1] * u.arcsec, frame=aia171_test_map.coordinate_frame)
    with pytest.raises(TypeError, match='WCSAxes'):
        aia171_test_map.draw_quadrangle(bottom_left, axes=ax)


def test_different_wcs_plot_warning(aia171_test_map, hmi_test_map):
    aia171_test_map.plot()
    with pytest.warns(SunpyUserWarning,
                      match=(r'The map world coordinate system \(WCS\) is different '
                             'from the axes WCS')):
        hmi_test_map.plot(axes=plt.gca())


@figure_test
def test_draw_limb_different_observer(aia171_test_map):
    # Create a new map from the test map with a different observer location
    new_map = copy.deepcopy(aia171_test_map)
    del new_map.meta['haex_obs']
    del new_map.meta['haey_obs']
    del new_map.meta['haez_obs']
    new_map.meta['hgln_obs'] = 45
    new_map.meta['hglt_obs'] = 10

    aia171_test_map.plot()
    new_map.draw_limb(color='red')


@figure_test
def test_draw_limb_heliographic_stonyhurst(aia171_test_map):
    # Create the WCS header for HGS axes
    header = {
        'date-obs': aia171_test_map.date.utc.isot,
        'mjd-obs': 55607.000004,
        'naxis': 2,
        'naxis1': 360,
        'naxis2': 180,
        'ctype1': 'HGLN-CAR',
        'ctype2': 'HGLT-CAR',
        'cdelt1': 1,
        'cdelt2': 1,
        'cunit1': 'deg',
        'cunit2': 'deg',
        'crpix1': 180.5,
        'crpix2': 90.5,
        'crval1': 0,
        'crval2': 0,
        'lonpole': 0,
    }
    fig = Figure()
    ax = fig.add_subplot(projection=WCS(header))
    aia171_test_map.draw_limb(axes=ax, color='red')
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    return fig


@figure_test
def test_plot_autoalign(aia171_test_map):
    aia171_test_map._data = aia171_test_map.data.astype('float32')
    rotated_map = aia171_test_map.rotate(30*u.deg, order=3)

    # Plotting the rotated map on the original projection should appear de-rotated
    fig = Figure()
    ax = fig.add_subplot(projection=aia171_test_map)
    rotated_map.plot(axes=ax, autoalign=True)
    return fig


def test_plot_autoalign_bad_inputs(aia171_test_map):
    with pytest.raises(ValueError):
        aia171_test_map.plot(autoalign='bad')


@figure_test
def test_plot_autoalign_pixel_alignment(aia171_test_map):
    # Verify that autoalign=True does not affect pixel alignment
    x, y = (z.value for z in aia171_test_map.reference_pixel)

    fig = Figure(figsize=(10, 4))

    ax1 = fig.add_subplot(121, projection=aia171_test_map)
    aia171_test_map.plot(axes=ax1, autoalign=False, title='autoalign=False')
    ax1.grid(False)
    ax1.set_xlim(x - 2, x + 2)
    ax1.set_ylim(y - 2, y + 2)

    ax2 = fig.add_subplot(122, projection=aia171_test_map)
    aia171_test_map.plot(axes=ax2, autoalign=True, title='autoalign=True')
    ax2.grid(False)
    ax2.set_xlim(x - 2, x + 2)
    ax2.set_ylim(y - 2, y + 2)

    return fig


def test_plot_unit8(aia171_test_map):
    # Check that plotting a map with uint8 data does not raise an error
    aia171_unit8 = sunpy.map.Map(aia171_test_map.data.astype('uint8'), aia171_test_map.meta)
    aia171_unit8.plot()
