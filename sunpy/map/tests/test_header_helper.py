import re

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates import frames, sun
from sunpy.map import make_fitswcs_header
from sunpy.map.header_helper import make_heliographic_header, make_hpr_header
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


def test_scale_conversion(map_data, hpc_coord):
    # The header will have cunit1/2 of arcsec
    header = make_fitswcs_header(map_data, hpc_coord, scale=[1, 2] * u.arcmin / u.pix)
    assert header['cdelt1'] == 60
    assert header['cdelt2'] == 120


def test_default_rotation(map_data, hpc_coord):
    header = make_fitswcs_header(map_data, hpc_coord)
    wcs = WCS(header)
    np.testing.assert_allclose(wcs.wcs.pc, [[1, 0], [0, 1]], atol=1e-5)


def test_rotation_angle(map_data, hpc_coord):
    header = make_fitswcs_header(map_data, hpc_coord,
                                 rotation_angle=90*u.deg)
    wcs = WCS(header)
    np.testing.assert_allclose(wcs.wcs.pc, [[0, -1], [1, 0]], atol=1e-5)


def test_rotation_angle_rectangular_pixels(map_data, hpc_coord):
    header = make_fitswcs_header(map_data, hpc_coord, scale=[2, 5] * u.arcsec / u.pix,
                                 rotation_angle=45*u.deg)
    wcs = WCS(header)
    np.testing.assert_allclose(wcs.wcs.pc, np.sqrt(0.5) * np.array([[1, -2.5], [0.4, 1]]), atol=1e-5)


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
    instrument_kwargs = {
        'instrument': 'test name',
        'observatory': 'test observatory',
        'telescope': 'test telescope',
        'detector': 'test detector',
        'wavelength': 171 * u.Angstrom,
        'exposure': 2 * u.s,
        'unit': u.Unit('ct s-1'),
    }
    header = make_fitswcs_header(map_data, hpc_coord, **instrument_kwargs)
    assert header['instrume'] == instrument_kwargs['instrument']
    assert header['obsrvtry'] == instrument_kwargs['observatory']
    assert header['telescop'] == instrument_kwargs['telescope']
    assert header['detector'] == instrument_kwargs['detector']
    assert header['wavelnth'] == instrument_kwargs['wavelength'].to_value()
    assert header['waveunit'] == instrument_kwargs['wavelength'].unit.to_string("fits")
    assert header['exptime'] == instrument_kwargs['exposure'].to_value('s')
    assert header['bunit'] == instrument_kwargs['unit'].to_string("fits")

    # Check returned MetaDict will make a `sunpy.map.Map`
    map_test = sunpy.map.Map(map_data, header)
    assert isinstance(map_test, sunpy.map.mapbase.GenericMap)


def test_quantity_input(map_data, hpc_coord):
    # Test that unit information is extracted correctly when data array is a quantity
    map_unit = u.Unit('ct / s')
    map_quantity = u.Quantity(map_data, map_unit)
    header = make_fitswcs_header(map_quantity, hpc_coord)
    assert header['bunit'] == map_unit.to_string('fits')
    # In cases where unit is specified and data is a quantity, specified unit will take
    # precedence
    override_unit = u.Unit('erg cm-2 s-1')
    header = make_fitswcs_header(map_quantity, hpc_coord, unit=override_unit)
    assert header['bunit'] == override_unit.to_string('fits')


def test_unit_as_string(map_data, hpc_coord):
    # Test that unit can be passed in as a string
    map_unit = u.Unit('ct / (pix s)')
    header = make_fitswcs_header(map_data, hpc_coord, unit=map_unit.to_string())
    assert header['bunit'] == map_unit.to_string(format='fits')


@pytest.mark.parametrize(('input_unit','output_string'),
                         [('DN', 'DN'), (u.DN, 'DN'), ('DN / s', 'DN / s')])
def test_make_fitswcs_header_handles_dn(input_unit, output_string, map_data, hpc_coord):
    header = make_fitswcs_header(map_data, hpc_coord, unit=input_unit)
    assert header['bunit'] == output_string

@pytest.mark.parametrize(
    ("coordinate_input", "expected_error_message"),
    [
        ("hcc_coord", "This function does not currently support heliocentric coordinates."),
        ("map_data", "coordinate needs to be a coordinate frame or an SkyCoord instance."),
        ("hpc_coord_notime", "The coordinate needs an observation time, `obstime`.")
    ]
)
def test_invalid_inputs(coordinate_input, expected_error_message, map_data, hcc_coord, hpc_coord_notime, hpc_coord):
    coordinates = {
        "hcc_coord": hcc_coord,
        "map_data": map_data,
        "hpc_coord_notime": hpc_coord_notime,
    }
    coord = coordinates[coordinate_input]


    with pytest.raises(ValueError, match=re.escape(expected_error_message)):
        make_fitswcs_header(map_data, coord)

    # Check arguments not given as astropy Quantities
    with pytest.raises(TypeError):
        make_fitswcs_header(map_data, hpc_coord, reference_pixel=[0, 0])
    with pytest.raises(TypeError):
        make_fitswcs_header(map_data, hpc_coord, scale=[0, 0])

    # Check arguments of reference_pixel and scale have to be given in astropy units of pix, and arcsec/pix
    with pytest.raises(u.UnitsError):
        make_fitswcs_header(map_data, hpc_coord, reference_pixel=u.Quantity([0, 0]))
    with pytest.raises(u.UnitsError):
        make_fitswcs_header(map_data, hpc_coord, scale=u.Quantity([0, 0]))
    with pytest.raises(u.UnitsError):
        make_fitswcs_header(map_data, hpc_coord, scale=u.Quantity([0, 0]*u.arcsec))


@pytest.mark.parametrize('frame', ['carrington', 'stonyhurst'])
# Second case here chosen to produce non-square pixels
@pytest.mark.parametrize('shape', [[90, 180], [240, 100]])
@pytest.mark.parametrize('projection_code', ['CAR', 'CEA'])
def test_make_heliographic_header(aia171_test_map, shape, projection_code, frame):
    header = make_heliographic_header(aia171_test_map.date, aia171_test_map.observer_coordinate, shape, frame=frame, projection_code=projection_code)
    carr_map = aia171_test_map.reproject_to(header)

    # Check upper right and lower left coordinates are as expected
    ll_coord = carr_map.pixel_to_world(-0.5 * u.pix, -0.5 * u.pix)
    assert ll_coord.lon in [-180 * u.deg, 180*u.deg]
    assert ll_coord.lat == -90 * u.deg
    assert ll_coord.frame.name == f"heliographic_{frame}"

    ur_coord = carr_map.pixel_to_world((shape[1] - 0.5) * u.pix, (shape[0] - 0.5) * u.pix)
    assert ur_coord.lon in [-180 * u.deg, 180*u.deg]
    assert ur_coord.lat == 90 * u.deg
    assert ur_coord.frame.name == f"heliographic_{frame}"


def test_make_heliographic_header_invalid_inputs(aia171_test_map):
    with pytest.raises(ValueError, match='projection_code must be one of'):
        make_heliographic_header(aia171_test_map.date, aia171_test_map.observer_coordinate, [90, 180], frame='carrington', projection_code='blah')

    with pytest.raises(ValueError, match='frame must be one of'):
        make_heliographic_header(aia171_test_map.date, aia171_test_map.observer_coordinate, [90, 180], frame='blah')

    with pytest.raises(TypeError):
        make_heliographic_header(aia171_test_map.date, aia171_test_map.observer_coordinate, [90, 180], frame='carrington', map_center_longitude=0)

    with pytest.raises(u.UnitsError):
        make_heliographic_header(aia171_test_map.date, aia171_test_map.observer_coordinate, [90, 180], frame='carrington', map_center_longitude=0*u.m)

    # Test new keyword propagates correctly to header reference pixel coordinate
    # hgc/hgs = heliographic carrington/stonyhurst
    header_test_hgs = make_heliographic_header(aia171_test_map.date, aia171_test_map.observer_coordinate, [90, 180], frame='carrington')
    header_test_hgc = make_heliographic_header(aia171_test_map.date, aia171_test_map.observer_coordinate, [90, 180], frame='carrington', map_center_longitude=180*u.deg)
    assert header_test_hgs['crval1'] == 0.0
    assert header_test_hgc['crval1'] == 180.0
    # Test keyword wraps correctly to range [0,360]
    header_test_above = make_heliographic_header(aia171_test_map.date, aia171_test_map.observer_coordinate, [90, 180], frame='carrington', map_center_longitude=361*u.deg)
    header_test_below = make_heliographic_header(aia171_test_map.date, aia171_test_map.observer_coordinate, [90, 180], frame='carrington', map_center_longitude=-1*u.deg)
    assert header_test_above['crval1'] == 1.0
    assert header_test_below['crval1'] == 359.0


def test_make_hpr_header(hgs_coord):
    header = make_hpr_header(hgs_coord, (100, 720), theta_binsize=0.1*u.deg)
    assert header['ctype1'] == 'HRLN-CAR'
    assert header['ctype2'] == 'HRLT-CAR'
    assert header['naxis1'] == 720
    assert header['naxis2'] == 100
    assert header['crpix1'] == 360.5
    assert header['crpix2'] == 900.5
    assert header['crval1'] == 180.0
    assert header['crval2'] == 0.0
    assert header['hgln_obs'] == -50.0
    assert header['hglt_obs'] == 50.0
    assert Time(header['date-obs']) == hgs_coord.obstime
