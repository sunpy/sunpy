import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.io
from sunpy.data.test import get_test_filepath
from sunpy.map import Map
from sunpy.tests.helpers import figure_test
from sunpy.visualization import drawing


@pytest.fixture
def aia171_test_map():
    return Map(get_test_filepath('aia_171_level1.fits'))


@pytest.fixture
def heliographic_test_map():
    (data, header), = sunpy.io._file_tools.read_file(get_test_filepath('heliographic_phase_map.fits.gz'))
    # Fix unit strings to prevent some astropy fits fixing warnings
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    # Set observer location to avoid warnings later
    header['HGLN_OBS'] = 0.0
    return Map((data, header))


@figure_test
def test_draw_equator_aia171(aia171_test_map):
    fig = plt.figure()
    axes = fig.add_subplot(projection=aia171_test_map)
    aia171_test_map.plot()
    drawing.equator(axes)


@figure_test
def test_draw_prime_meridian_aia171(aia171_test_map):
    fig = plt.figure()
    axes = fig.add_subplot(projection=aia171_test_map)
    aia171_test_map.plot()
    drawing.prime_meridian(axes)


@figure_test
def test_heliographic_equator_prime_meridian(heliographic_test_map):
    fig = plt.figure()
    axes = fig.add_subplot(projection=heliographic_test_map)
    heliographic_test_map.plot()
    drawing.equator(axes, color="blue")
    drawing.prime_meridian(axes, color="red")


def test_prime_meridian_error():
    axes = plt.subplot(projection=WCS())
    with pytest.raises(ValueError, match="does not have an observer"):
        drawing.prime_meridian(axes)


def test_limb_invisible(aia171_test_map):
    aia_obs = aia171_test_map.observer_coordinate
    # Create a new observer on the opposite side of the Sun
    new_obs = SkyCoord(lon=aia_obs.lon + 180*u.deg,
                       lat=-aia_obs.lat,
                       radius=aia_obs.radius / 10,
                       frame=aia_obs.replicate_without_data())

    ax = Figure().add_subplot(111, projection=aia171_test_map)
    # Original observer
    visible, hidden = drawing.limb(ax, aia_obs)
    assert visible is not None
    assert hidden is None
    # Far side observer
    visible, hidden = drawing.limb(ax, new_obs)
    assert visible is None
    assert hidden is not None
