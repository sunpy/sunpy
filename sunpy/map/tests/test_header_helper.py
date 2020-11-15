import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates import frames
from sunpy.util.metadata import MetaDict


# test setups
@pytest.fixture
def map_data():
    return np.random.rand(10, 10)


@pytest.fixture
def hpc_test_header():
    return SkyCoord(0*u.arcsec, 0*u.arcsec, observer='earth',
                    obstime='2013-10-28 00:00', frame=frames.Helioprojective)


@pytest.fixture
def hgc_test_header():
    return SkyCoord(70*u.deg, -30*u.deg, observer='earth',
                    obstime='2013-10-28 00:00', frame=frames.HeliographicCarrington)


@pytest.fixture
def hgs_test_header():
    return SkyCoord(-50*u.deg, 50*u.deg, observer='earth',
                    obstime='2013-10-28 00:00', frame=frames.HeliographicStonyhurst)


@pytest.fixture
def hcc_test_header():
    return SkyCoord(-72241*u.km, 361206.1*u.km, 589951.4*u.km,
                    obstime='2013-10-28 00:00', frame=frames.Heliocentric)


@pytest.fixture
def hpc_test_header_notime():
    return SkyCoord(0*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective)


# tests
def test_metakeywords():
    meta = sunpy.map.meta_keywords()
    assert isinstance(meta, dict)


def test_rotation_angle(map_data, hpc_test_header):
    header = sunpy.map.make_fitswcs_header(map_data, hpc_test_header,
                                           rotation_angle=90*u.deg)
    wcs = WCS(header)
    np.testing.assert_allclose(wcs.wcs.pc, [[0, -1], [1, 0]], atol=1e-5)


def test_rotation_matrix(map_data, hpc_test_header):
    header = sunpy.map.make_fitswcs_header(map_data, hpc_test_header,
                                           rotation_matrix=np.array([[1, 0], [0, 1]]))
    wcs = WCS(header)
    np.testing.assert_allclose(wcs.wcs.pc, [[1, 0], [0, 1]], atol=1e-5)


def test_make_fits_header(map_data, hpc_test_header, hgc_test_header,
                          hgs_test_header, hcc_test_header, hpc_test_header_notime):

    # Check that different coordinate frames return header MetaDict or not in the case of HCC
    assert isinstance(sunpy.map.make_fitswcs_header(map_data, hpc_test_header), MetaDict)
    assert isinstance(sunpy.map.make_fitswcs_header(map_data, hgc_test_header), MetaDict)
    assert isinstance(sunpy.map.make_fitswcs_header(map_data, hgs_test_header), MetaDict)
    # Raise the HCC error
    with pytest.raises(ValueError):
        sunpy.map.make_fitswcs_header(map_data, hcc_test_header)

    # Check for when coordinate argument isn't given as an `astropy.coordinate.SkyCoord`
    with pytest.raises(ValueError):
        sunpy.map.make_fitswcs_header(map_data, map_data)

    # Check for when an observation time isn't given
    with pytest.raises(ValueError):
        sunpy.map.make_fitswcs_header(map_data, hpc_test_header_notime)

    # Check that correct information is in header MetaDict including observer for HPC
    header = sunpy.map.make_fitswcs_header(map_data, hpc_test_header)
    assert header['crval1'] == 0
    assert header['crpix1'] == 5.5
    assert header['ctype1'] == 'HPLN-TAN'
    assert u.allclose(header['dsun_obs'], hpc_test_header.frame.observer.radius.to_value(u.m))
    assert isinstance(WCS(header), WCS)

    # Check no observer info for HGS
    header = sunpy.map.make_fitswcs_header(map_data, hgs_test_header)
    assert header.get('dsun_obs') is None
    assert isinstance(WCS(header), WCS)

    # Check for observer info for HGC
    header = sunpy.map.make_fitswcs_header(map_data, hgc_test_header)
    assert u.allclose(header['dsun_obs'], hgc_test_header.frame.observer.radius.to_value(u.m))
    assert isinstance(WCS(header), WCS)

    # Check arguments not given as astropy Quantities
    with pytest.raises(TypeError):
        header = sunpy.map.make_fitswcs_header(map_data, hpc_test_header, reference_pixel=[0, 0])
        header = sunpy.map.make_fitswcs_header(map_data, hpc_test_header, scale=[0, 0])

    # Check arguments of reference_pixel and scale have to be given in astropy units of pix, and arcsec/pix
    with pytest.raises(u.UnitsError):
        header = sunpy.map.make_fitswcs_header(
            map_data, hpc_test_header, reference_pixel=u.Quantity([0, 0]))
        header = sunpy.map.make_fitswcs_header(map_data, hpc_test_header, scale=u.Quantity([0, 0]))
        header = sunpy.map.make_fitswcs_header(
            map_data, hpc_test_header, scale=u.Quantity([0, 0]*u.arcsec))

    # Check keyword helper arguments
    header = sunpy.map.make_fitswcs_header(map_data, hpc_test_header, instrument='test name')
    assert header['instrume'] == 'test name'

    # Check returned MetaDict will make a `sunpy.map.Map`
    map_test = sunpy.map.Map(map_data, header)
    assert isinstance(map_test, sunpy.map.mapbase.GenericMap)


def test_HGS_CAR_header():
    # This tests both non-HPC and non-TAN header generation.
    new_data = np.empty((72, 144))
    new_frame = SkyCoord(0*u.deg, 0*u.deg, obstime="2019-06-16", frame="heliographic_stonyhurst")
    new_header = sunpy.map.make_fitswcs_header(new_data, new_frame,
                                               scale=[2.5, 2.5]*u.deg/u.pix,
                                               projection_code="CAR")

    assert new_header['ctype1'] == "HGLN-CAR"
    assert new_header['ctype2'] == "HGLT-CAR"
    assert new_header['cunit1'] == "deg"
    assert new_header['cunit2'] == "deg"

    assert "hgln_obs" not in new_header
    assert "hglt_obs" not in new_header
    assert "dsun_obs" not in new_header

    wcs = WCS(new_header)
    pix_coords = np.mgrid[0:new_data.shape[1], 0:new_data.shape[0]]
    world_coords = wcs.pixel_to_world(*pix_coords)

    # In this projection the coordinate system is symmetric around the
    # reference pixel. because the reference pixel should be in the center of
    # the array it means the array should be symmetric apart from the sign
    # flip. This is therefore testing that the reference pixel is correctly in
    # the center of the array.
    assert u.allclose(world_coords.lon, world_coords.lon[::-1, ::-1]*-1)
    assert u.allclose(world_coords.lat, world_coords.lat[::-1, ::-1]*-1)


def test_latlon_order():
    # This tests that CRVAL1, CRVAL2 (i.e. lon,lat) correspond to CTYPE1, CTYPE2 (i.e. HGLN, HGLT)
    data = np.zeros((100, 100))
    coord = SkyCoord(20*u.arcsec, -10*u.arcsec, obstime='2013-10-28', frame=frames.Helioprojective)
    header = sunpy.map.make_fitswcs_header(data, coord)
    # check LON
    assert ('LN' in header['ctype1']) and header['crval1'] == coord.spherical.lon.to_value()
    # check LAT
    assert 'LT' in header['ctype2'] and header['crval2'] == coord.spherical.lat.to_value()
