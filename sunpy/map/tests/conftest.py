import os

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

import sunpy.coordinates
import sunpy.data.test
import sunpy.io
import sunpy.map
import sunpy.sun
from sunpy.coordinates import sun

testpath = sunpy.data.test.rootdir


@pytest.fixture
def test_map(request):
    return request.getfixturevalue(request.param)


@pytest.fixture
def hmi_test_map():
    (data, header), = sunpy.io.read_file(os.path.join(testpath, 'resampled_hmi.fits'))

    # Get rid of the blank keyword to prevent some astropy fits fixing warnings
    header.pop('CRDER2')
    header.pop('CRDER1')
    return sunpy.map.Map((data, header))


@pytest.fixture
def aia171_test_map():
    (data, header), = sunpy.io.read_file(os.path.join(testpath, 'aia_171_level1.fits'))

    # Get rid of the blank keyword to prevent some astropy fits fixing warnings
    header.pop('BLANK')
    return sunpy.map.Map((data, header))


@pytest.fixture
def aia171_roll_map(aia171_test_map):
    return aia171_test_map.rotate(-45*u.deg)


@pytest.fixture
def heliographic_test_map():
    (data, header), = sunpy.io.read_file(os.path.join(testpath, 'heliographic_phase_map.fits.gz'))

    # Fix unit strings to prevent some astropy fits fixing warnings
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    # Set observer location to avoid warnings later
    header['HGLN_OBS'] = 0.0
    return sunpy.map.Map((data, header))


@pytest.fixture
def aia171_test_map_with_mask(aia171_test_map):
    shape = aia171_test_map.data.shape
    mask = np.zeros_like(aia171_test_map.data, dtype=bool)
    mask[0:shape[0] // 2, 0:shape[1] // 2] = True
    return sunpy.map.Map(np.ma.array(aia171_test_map.data, mask=mask), aia171_test_map.meta)


@pytest.fixture
def generic_map():
    data = np.ones([6, 6], dtype=np.float64)
    dobs = Time('1970-01-01T00:00:00')
    l0 = sun.L0(dobs).to_value(u.deg)
    b0 = sun.B0(dobs).to_value(u.deg)
    dsun = sun.earth_distance(dobs).to_value(u.m)
    header = {
        'CRVAL1': 0,
        'CRVAL2': 0,
        'CRPIX1': 5,
        'CRPIX2': 5,
        'CDELT1': 10,
        'CDELT2': 10,
        'CUNIT1': 'arcsec',
        'CUNIT2': 'arcsec',
        'CTYPE1': 'HPLN-TAN',
        'CTYPE2': 'HPLT-TAN',
        'PC1_1': 0,
        'PC1_2': -1,
        'PC2_1': 1,
        'PC2_2': 0,
        'NAXIS1': 6,
        'NAXIS2': 6,
        'date-obs': dobs.isot,
        'crln_obs': l0,
        'crlt_obs': b0,
        "dsun_obs": dsun,
        'mjd-obs': 40587.0,
        'obsrvtry': 'Foo',
        'detector': 'bar',
        'wavelnth': 10,
        'waveunit': 'm',
        'bunit': 'ct/s',
    }
    return sunpy.map.Map((data, header))


@pytest.fixture
def simple_map():
    # A 3x3 map, with it's center at (0, 0), and scaled differently in
    # each direction
    data = np.arange(9).reshape((3, 3))
    ref_coord = SkyCoord(0.0, 0.0, frame='helioprojective', obstime='now', unit='deg',
                         observer=SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.AU,
                                           frame='heliographic_stonyhurst'))
    ref_pix = [1, 1] * u.pix
    scale = [2, 1] * u.arcsec / u.pix
    header = sunpy.map.make_fitswcs_header(data, ref_coord, reference_pixel=ref_pix, scale=scale)
    return sunpy.map.Map(data, header)


@pytest.fixture
def eit_test_map():
    """
    Load SunPy's test EIT image.
    """
    testpath = sunpy.data.test.rootdir
    eit_file = os.path.join(testpath, "EIT", "efz20040301.020010_s.fits")
    return sunpy.map.Map(eit_file)
