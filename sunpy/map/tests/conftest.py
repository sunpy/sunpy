import pathlib
import warnings

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.time import Time

import sunpy.coordinates
import sunpy.map
import sunpy.sun
from sunpy.coordinates import sun
from sunpy.data.test import get_dummy_map_from_header, get_test_filepath
from sunpy.io._file_tools import read_file


@pytest.fixture
def test_map(request):
    return request.getfixturevalue(request.param)


@pytest.fixture
def hmi_test_map():
    (data, header), = read_file(get_test_filepath('resampled_hmi.fits'))

    # Get rid of the blank keyword to prevent some astropy fits fixing warnings
    header.pop('BLANK')
    header.pop('CRDER2')
    header.pop('CRDER1')
    return sunpy.map.Map((data, header))


@pytest.fixture
def aia171_test_map():
    (data, header), = read_file(get_test_filepath('aia_171_level1.fits'))

    # Get rid of the blank keyword to prevent some astropy fits fixing warnings
    header.pop('BLANK')
    return sunpy.map.Map((data, header))


@pytest.fixture
def adjusted_test_maps(aia171_test_map, hmi_test_map):
    # The test maps have wildly different observation times, which throws off compositing
    hmi_test_map.meta['date-obs'] = aia171_test_map.meta['date-obs']
    hmi_test_map.meta['t_obs'] = aia171_test_map.meta['t_obs']
    # Also set the HMI observer location to be the same as the AIA observer location
    del hmi_test_map.meta['crln_obs']
    del hmi_test_map.meta['crlt_obs']
    hmi_test_map.meta['hgln_obs'] = aia171_test_map.observer_coordinate.lon.to_value('deg')
    hmi_test_map.meta['hglt_obs'] = aia171_test_map.observer_coordinate.lat.to_value('deg')
    return aia171_test_map, hmi_test_map


@pytest.fixture
def aia171_roll_map(aia171_test_map):
    return aia171_test_map.rotate(-45*u.deg)


@pytest.fixture
def heliographic_test_map():
    (data, header), = read_file(get_test_filepath('heliographic_phase_map.fits.gz'))

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
    data = np.arange(36, dtype=np.float64).reshape((6, 6))
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
        'bunit': 'DN/s',
    }
    return sunpy.map.Map((data, header))


@pytest.fixture
def simple_map():
    """
    A simple 9x9 map, with its center at (0, 0),
    and scaled differently in each direction.
    """
    data = np.arange(81).reshape((9, 9))
    ref_coord = SkyCoord(0.0, 0.0, frame='helioprojective', obstime='2020-01-01 00:00:00', unit='deg',
                         observer=SkyCoord(0 * u.deg, 0 * u.deg, 1 * u.AU, frame='heliographic_stonyhurst'))
    ref_pix = [4, 4] * u.pix
    scale = [2, 1] * u.arcsec / u.pix
    header = sunpy.map.make_fitswcs_header(data, ref_coord, reference_pixel=ref_pix, scale=scale)
    return sunpy.map.Map(data, header)


@pytest.fixture
def carrington_map():
    # This is a 20 x 20 map in a Carrington frame, with the reference pixel *not* at the
    # equator. This results in a non-linear transformation between pixel and world
    # coordinates, so is ideal for testing situations where the non-linearity matters
    data = np.arange(20**2).reshape((20, 20))
    obstime = '2020-01-01'
    observer = SkyCoord(0*u.deg, 0*u.deg, frame='heliographic_stonyhurst', obstime=obstime)
    ref_coord = SkyCoord(120*u.deg, -70*u.deg, frame='heliographic_carrington', obstime=obstime, observer=observer)
    ref_pix = [0, 0] * u.pix
    scale = [2, 1] * u.deg / u.pix
    header = sunpy.map.make_fitswcs_header(data, ref_coord, reference_pixel=ref_pix, scale=scale)
    return sunpy.map.Map(data, header)


@pytest.fixture
def eit_test_map():
    """
    Load SunPy's test EIT image.
    """
    return get_dummy_map_from_header(get_test_filepath("EIT_header/efz20040301.020010_s.header"))


@pytest.fixture
def sample_171():
    from sunpy.data.sample import AIA_171_IMAGE

    # Need to bypass
    # VerifyWarning: Invalid 'BLANK' keyword in header.
    # The 'BLANK' keyword is only applicable to integer data, and will be ignored in this HDU.
    with fits.open(AIA_171_IMAGE) as hdu:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=VerifyWarning)
            return sunpy.map.Map(hdu[1].data, hdu[1].header)


@pytest.fixture
def sample_hmi():
    from sunpy.data.sample import HMI_LOS_IMAGE

    return sunpy.map.Map(HMI_LOS_IMAGE)


@pytest.fixture
def swap_test_map():
    return get_dummy_map_from_header(get_test_filepath("swap_lv1_20140606_000113.header"))


@pytest.fixture
def rhessi_test_map():
    return sunpy.map.Map(get_test_filepath("hsi_image_20101016_191218.fits"))


@pytest.fixture
def lasco_test_map():
    return get_dummy_map_from_header(get_test_filepath("lasco_c2_25299383_s.header"))


ALL_MAP_SOURCE_FILES = [
    "aia_171_level1.fits",
    "resampled_hmi.fits",
    "hsi_image_20101016_191218.fits",
    "adapt.header",
    "FGMG4_20110214_030443.7.header",
    "HinodeXRT.header",
    "iris_l2_20130801_074720_4040000014_SJI_1400_t000.header",
    "20181209_180305_kcor_l2.header",
    "swap_lv1_20140606_000113.header",
    "tsi20010130_025823_a2.header",
    "punch.header",
    "lasco_c2_25299383_s.header",
    "mdi.fd_Ic.20101015_230100_TAI.data.header",
    "solo_L1_eui-fsi304-image_20201021T145510206_V03.header",
    "cor1_20090615_000500_s4c1A.header",
    "euvi_20090615_000900_n4euA_s.header",
    "hi_20110910_114721_s7h2A.header",
    "SUT_T24_0847_000444_Lev1.0_2024-06-28T18.21.33.178_0971NB03.header",
    "dr_suvi-l2-ci195_g16_s20190403T093200Z_e20190403T093600Z_v1-0-0_rebinned.header",
    "YohkohSXT.header",
    "gong_halpha.header",
    "gong_synoptic.header",
    "EIT_header/efz20040301.020010_s.header",
    "hxi_imgcube_01e02t_20230501_130758_HXI_CLEAN_0.header",
    "solo_L2_metis-vl-tb_20220322T211301_V01.header",
    "solo_L2_phi-hrt-icnt_20220307T000009_V202208311927_0243070101.header",
]


@pytest.fixture(params=ALL_MAP_SOURCE_FILES)
def all_test_maps(request):
    """
    Fixture parameterized across all supported map sources.
    """
    p = pathlib.Path(get_test_filepath(request.param))
    return get_dummy_map_from_header(p) if p.suffix == ".header" else sunpy.map.Map(p)

