import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates import frames, sun
from sunpy.map import make_fitswcs_header
from sunpy.util.metadata import MetaDict


@pytest.fixture
def map_data():
    return np.random.rand(20, 10)


@pytest.fixture
def hpc_coord():
    return SkyCoord(0*u.arcsec, 100*u.arcsec, observer='earth', rsun=7.1e8*u.m,
                    obstime='2013-10-28 00:00', frame=frames.Helioprojective)


@pytest.fixture
def hgc_coord():
    return SkyCoord(70*u.deg, -30*u.deg, 1.5*u.AU, observer='self', rsun=7.2e8*u.m,
                    obstime='2013-10-28 00:00', frame=frames.HeliographicCarrington)


@pytest.fixture
def hgs_coord():
    return SkyCoord(-50*u.deg, 50*u.deg, rsun=7.3e8*u.m,
                    obstime='2013-10-28 00:00', frame=frames.HeliographicStonyhurst)


@pytest.fixture
def hcc_coord():
    return SkyCoord(-72241*u.km, 361206.1*u.km, 589951.4*u.km, observer='earth',
                    obstime='2013-10-28 00:00', frame=frames.Heliocentric)


@pytest.fixture
def hpc_header(map_data, hpc_coord):
    return make_fitswcs_header(map_data, hpc_coord)


@pytest.fixture
def hgc_header(map_data, hgc_coord):
    return make_fitswcs_header(map_data, hgc_coord, projection_code='CAR')


@pytest.fixture
def hgs_header(map_data, hgs_coord):
    return make_fitswcs_header(map_data, hgs_coord, projection_code='CAR')


@pytest.fixture
def hpc_coord_notime():
    return SkyCoord(0*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective)


def test_metakeywords():
    meta = sunpy.map.meta_keywords()
    assert isinstance(meta, dict)


def test_rotation_angle(map_data, hpc_coord):
    header = make_fitswcs_header(map_data, hpc_coord,
                                 rotation_angle=90*u.deg)
    wcs = WCS(header)
    np.testing.assert_allclose(wcs.wcs.pc, [[0, -1], [1, 0]], atol=1e-5)


def test_rotation_matrix(map_data, hpc_coord):
    header = make_fitswcs_header(map_data, hpc_coord,
                                 rotation_matrix=np.array([[1, 0], [0, 1]]))
    wcs = WCS(header)
    np.testing.assert_allclose(wcs.wcs.pc, [[1, 0], [0, 1]], atol=1e-5)


def test_hpc_header(hpc_header, hpc_coord):
    assert isinstance(hpc_header, MetaDict)

    assert hpc_header['naxis1'] == 10
    assert hpc_header['naxis2'] == 20
    assert hpc_header['crval1'] == 0
    assert hpc_header['crpix1'] == 5.5
    assert hpc_header['ctype1'] == 'HPLN-TAN'
    assert hpc_header['crval2'] == 100.
    assert hpc_header['crpix2'] == 10.5
    assert hpc_header['ctype2'] == 'HPLT-TAN'

    assert hpc_header['lonpole'] == 180.
    assert u.allclose(hpc_header['rsun_ref'] * u.m, hpc_coord.rsun)

    # Check for observer info for HPC
    assert u.allclose(hpc_header['hgln_obs'] * u.deg, hpc_coord.observer.lon)
    assert u.allclose(hpc_header['hglt_obs'] * u.deg, hpc_coord.observer.lat)
    assert u.allclose(hpc_header['dsun_obs'] * u.m, hpc_coord.observer.radius)
    assert u.allclose(hpc_header['rsun_obs'] * u.arcsec,
                      sun._angular_radius(hpc_coord.rsun, hpc_coord.observer.radius))

    assert isinstance(WCS(hpc_header), WCS)


def test_hgc_header(hgc_header, hgc_coord):
    assert isinstance(hgc_header, MetaDict)

    assert hgc_header['naxis1'] == 10
    assert hgc_header['naxis2'] == 20
    assert hgc_header['crval1'] == 70
    assert hgc_header['crpix1'] == 5.5
    assert hgc_header['ctype1'] == "CRLN-CAR"
    assert hgc_header['crval2'] == -30
    assert hgc_header['crpix2'] == 10.5
    assert hgc_header['ctype2'] == "CRLT-CAR"
    assert hgc_header['cunit1'] == "deg"
    assert hgc_header['cunit2'] == "deg"

    assert hgc_header['lonpole'] == 180.  # for negative reference latitude in a CAR projection
    assert u.allclose(hgc_header['rsun_ref'] * u.m, hgc_coord.rsun)

    # Check for observer info for HGC (which is set to "self")
    assert u.allclose(hgc_header['crln_obs'] * u.deg, hgc_coord.lon)
    assert u.allclose(hgc_header['crlt_obs'] * u.deg, hgc_coord.lat)
    assert u.allclose(hgc_header['dsun_obs'] * u.m, hgc_coord.radius)
    assert u.allclose(hgc_header['rsun_obs'] * u.arcsec,
                      sun._angular_radius(hgc_coord.rsun, hgc_coord.radius))

    assert isinstance(WCS(hgc_header), WCS)


def test_hgs_header(hgs_header, hgs_coord):
    assert isinstance(hgs_header, MetaDict)

    assert hgs_header['naxis1'] == 10
    assert hgs_header['naxis2'] == 20
    assert hgs_header['crval1'] == -50
    assert hgs_header['crpix1'] == 5.5
    assert hgs_header['ctype1'] == "HGLN-CAR"
    assert hgs_header['crval2'] == 50
    assert hgs_header['crpix2'] == 10.5
    assert hgs_header['ctype2'] == "HGLT-CAR"
    assert hgs_header['cunit1'] == "deg"
    assert hgs_header['cunit2'] == "deg"

    assert hgs_header['lonpole'] == 0.  # for positive reference latitude in a CAR projection
    assert u.allclose(hgs_header['rsun_ref'] * u.m, hgs_coord.rsun)

    # Check no observer info for HGS
    assert "hgln_obs" not in hgs_header
    assert "hglt_obs" not in hgs_header
    assert 'dsun_obs' not in hgs_header
    assert 'rsun_obs' not in hgs_header

    assert isinstance(WCS(hgs_header), WCS)


def test_instrument_keyword(map_data, hpc_coord):
    header = make_fitswcs_header(map_data, hpc_coord, instrument='test name')
    assert header['instrume'] == 'test name'

    # Check returned MetaDict will make a `sunpy.map.Map`
    map_test = sunpy.map.Map(map_data, header)
    assert isinstance(map_test, sunpy.map.mapbase.GenericMap)


def test_invalid_inputs(map_data, hcc_coord, hpc_coord_notime, hpc_coord):
    # Raise the HCC error
    with pytest.raises(ValueError):
        make_fitswcs_header(map_data, hcc_coord)

    # Check for when coordinate argument isn't given as an `astropy.coordinate.SkyCoord`
    with pytest.raises(ValueError):
        make_fitswcs_header(map_data, map_data)

    # Check for when an observation time isn't given
    with pytest.raises(ValueError):
        make_fitswcs_header(map_data, hpc_coord_notime)

    # Check arguments not given as astropy Quantities
    with pytest.raises(TypeError):
        header = make_fitswcs_header(map_data, hpc_coord, reference_pixel=[0, 0])
        header = make_fitswcs_header(map_data, hpc_coord, scale=[0, 0])

    # Check arguments of reference_pixel and scale have to be given in astropy units of pix, and arcsec/pix
    with pytest.raises(u.UnitsError):
        header = make_fitswcs_header(map_data, hpc_coord, reference_pixel=u.Quantity([0, 0]))
        header = make_fitswcs_header(map_data, hpc_coord, scale=u.Quantity([0, 0]))
        header = make_fitswcs_header(map_data, hpc_coord, scale=u.Quantity([0, 0]*u.arcsec))
