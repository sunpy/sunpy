from textwrap import dedent

import pytest
import numpy as np

import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits

from sunpy.tests.helpers import figure_test
from sunpy.visualization.animator.wcs import ArrayAnimatorWCS


@pytest.fixture
def wcs_4d():
    # This is weirdly indented and with line continuations to make it a valid
    # FITS header string so that no warnings are logged during WCS
    # construction.
    header = dedent("""\
        WCSAXES =                    4 / Number of coordinate axes
        CRPIX1  =                  0.0 / Pixel coordinate of reference point
        CRPIX2  =                  0.0 / Pixel coordinate of reference point
        CRPIX3  =                  0.0 / Pixel coordinate of reference point
        CRPIX4  =                  5.0 / Pixel coordinate of reference point
        CDELT1  =                  0.4 / [min] Coordinate increment at reference point
        CDELT2  =                2E-11 / [m] Coordinate increment at reference point
        CDELT3  =   0.0027777777777778 / [deg] Coordinate increment at reference point
        CDELT4  =   0.0013888888888889 / [deg] Coordinate increment at reference point
        CUNIT1  = 'min'                / Units of coordinate increment and value
        CUNIT2  = 'm'                  / Units of coordinate increment and value
        CUNIT3  = 'deg'                / Units of coordinate increment and value
        CUNIT4  = 'deg'                / Units of coordinate increment and value
        CTYPE1  = 'TIME'               / Coordinate type code
        CTYPE2  = 'WAVE'               / Vacuum wavelength (linear)
        CTYPE3  = 'HPLT-TAN'           / Coordinate type codegnomonic projection
        CTYPE4  = 'HPLN-TAN'           / Coordinate type codegnomonic projection
        CRVAL1  =                  0.0 / [min] Coordinate value at reference point
        CRVAL2  =                  0.0 / [m] Coordinate value at reference point
        CRVAL3  =                  0.0 / [deg] Coordinate value at reference point
        CRVAL4  =                  0.0 / [deg] Coordinate value at reference point
        LONPOLE =                180.0 / [deg] Native longitude of celestial pole
        LATPOLE =                  0.0 / [deg] Native latitude of celestial pole
        """)
    return WCS(header=fits.Header.fromstring(header, sep='\n'))


@pytest.mark.parametrize("data, slices, dim", (
    (np.arange(120).reshape((5, 4, 3, 2)), [0, 0, 'x', 'y'], 2),
    (np.arange(120).reshape((5, 4, 3, 2)), [0, 'x', 0, 'y'], 2),
    (np.arange(120).reshape((5, 4, 3, 2)), ['x', 0, 0, 'y'], 2),
    (np.arange(120).reshape((5, 4, 3, 2)), ['y', 0, 'x', 0], 2),
    (np.arange(120).reshape((5, 4, 3, 2)), [0, 0, 0, 'x'], 1),
))
def test_construct_array_animator(wcs_4d, data, slices, dim):
    if dim == 1:
        pytest.importorskip("astropy", minversion="4.0.0")

    array_animator = ArrayAnimatorWCS(data, wcs_4d, slices)

    assert isinstance(array_animator, ArrayAnimatorWCS)
    assert array_animator.plot_dimensionality == dim
    assert array_animator.num_sliders == data.ndim - dim
    for i, (wslice, arange) in enumerate(zip(slices, array_animator.axis_ranges[::-1])):
        if wslice not in ['x', 'y']:
            assert callable(arange)
            a = arange(0)
            if "pos" in wcs_4d.world_axis_physical_types[i]:
                assert isinstance(a, u.Quantity)
                assert u.allclose(a, 0 * u.pix)
            else:
                assert isinstance(a, u.Quantity)
                assert a.value == wcs_4d.pixel_to_world_values(*[0]*wcs_4d.world_n_dim)[i]
                assert a.unit == wcs_4d.world_axis_units[i]
        else:
            assert arange is None


def test_constructor_errors(wcs_4d):
    # WCS is not BaseLowLevelWCS
    with pytest.raises(ValueError, match="provided that implements the astropy WCS API."):
        ArrayAnimatorWCS(np.arange(25).reshape((5, 5)), {}, ['x', 'y'])

    # Data has wrong number of dimensions
    with pytest.raises(ValueError, match="Dimensionality of the data and WCS object do not match."):
        ArrayAnimatorWCS(np.arange(25).reshape((5, 5)), wcs_4d, ['x', 'y'])

    # Slices is wrong length
    with pytest.raises(ValueError, match="slices should be the same length"):
        ArrayAnimatorWCS(np.arange(16).reshape((2, 2, 2, 2)), wcs_4d, ['x', 'y'])

    # x not in slices
    with pytest.raises(ValueError, match="slices should contain at least"):
        ArrayAnimatorWCS(np.arange(16).reshape((2, 2, 2, 2)), wcs_4d, [0, 0, 0, 'y'])


@figure_test
def test_array_animator_wcs_2d_simple_plot(wcs_4d):
    data = np.arange(120).reshape((5, 4, 3, 2))
    a = ArrayAnimatorWCS(data, wcs_4d, [0, 0, 'x', 'y'])
    return a.fig


@figure_test
def test_array_animator_wcs_2d_update_plot(wcs_4d):
    data = np.arange(120).reshape((5, 4, 3, 2))
    a = ArrayAnimatorWCS(data, wcs_4d, [0, 0, 'x', 'y'])
    a.update_plot(1, a.im, a.sliders[0]._slider)
    return a.fig


@figure_test
def test_array_animator_wcs_2d_transpose_update_plot(wcs_4d):
    data = np.arange(120).reshape((5, 4, 3, 2))
    a = ArrayAnimatorWCS(data, wcs_4d, [0, 0, 'y', 'x'], colorbar=True)
    a.update_plot(1, a.im, a.sliders[0]._slider)
    return a.fig


@figure_test
def test_array_animator_wcs_2d_colorbar_buttons(wcs_4d):
    data = np.arange(120).reshape((5, 4, 3, 2))
    bf = [lambda x: x]*10
    bl = ['h']*10
    a = ArrayAnimatorWCS(data, wcs_4d, [0, 0, 'y', 'x'],
                         colorbar=True, button_func=bf, button_labels=bl)
    a.update_plot(1, a.im, a.sliders[0]._slider)
    return a.fig


@figure_test
def test_array_animator_wcs_2d_colorbar_buttons_default_labels(wcs_4d):
    data = np.arange(120).reshape((5, 4, 3, 2))
    bf = [lambda x: x]*10
    bl = ['h']*10
    a = ArrayAnimatorWCS(data, wcs_4d, [0, 0, 'y', 'x'], colorbar=True, button_func=bf)
    a.update_plot(1, a.im, a.sliders[0]._slider)
    return a.fig


@figure_test
def test_array_animator_wcs_2d_extra_sliders(wcs_4d):
    def vmin_slider(val, im, slider):
        im.set_clim(vmin=val)

    def vmax_slider(val, im, slider):
        im.set_clim(vmax=val)

    data = np.arange(120).reshape((5, 4, 3, 2))
    a = ArrayAnimatorWCS(data, wcs_4d, [0, 0, 'y', 'x'], colorbar=True,
                         slider_functions=[vmin_slider, vmax_slider],
                         slider_ranges=[[0, 100], [0, 100]])
    a.update_plot(1, a.im, a.sliders[0]._slider)
    return a.fig


@figure_test
def test_array_animator_wcs_1d_update_plot(wcs_4d):
    pytest.importorskip("astropy", minversion="4.0.0")
    data = np.arange(120).reshape((5, 4, 3, 2))
    a = ArrayAnimatorWCS(data, wcs_4d, [0, 0, 'x', 0], ylabel="Y axis!")
    a.update_plot(1, a.im, a.sliders[0]._slider)
    return a.fig


@figure_test
def test_array_animator_wcs_coord_params(wcs_4d):

    coord_params = {
        'hpln': {
            'format_unit': u.deg,
            'major_formatter': 'hh:mm:ss',
            'axislabel': 'Longitude',
            'ticks': {'spacing': 10*u.arcsec}
        }
    }

    data = np.arange(120).reshape((5, 4, 3, 2))
    a = ArrayAnimatorWCS(data, wcs_4d, [0, 0, 'x', 'y'], coord_params=coord_params)
    return a.fig


@figure_test
def test_array_animator_wcs_coord_params_grid(wcs_4d):
    pytest.importorskip("astropy", minversion="4.0.0")

    coord_params = {
        'hpln': {
            'format_unit': u.deg,
            'major_formatter': 'hh:mm:ss',
            'axislabel': 'Longitude',
            'grid': True
        }
    }

    data = np.arange(120).reshape((5, 4, 3, 2))
    a = ArrayAnimatorWCS(data, wcs_4d, [0, 0, 'x', 'y'], coord_params=coord_params)
    return a.fig
