from regions import PointSkyRegion, PolygonSkyRegion

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sunpy.net.hek.utils import get_unit, parse_chaincode


def test_get_unit():
    erg_per_cm3 = u.def_unit("ergs/cm^3", u.erg/u.ml)

    u1 = get_unit('DN/sec/pixel')
    u2 = get_unit('ergs per cubic centimeter')
    u3 = get_unit('m/s/s')

    assert u1 == u.Unit('DN / (pix s)')
    assert u2 == erg_per_cm3
    assert u3 == u.Unit('m / s2')

def test_parse_chaincode_helioprojective():
    value = "POLYGON((10 20, 30 40, 50 60))"
    attribute = {"frame": "helioprojective"}
    unit = "arcsec"
    time = Time('2011-08-09 12:04:00.000')

    result = parse_chaincode(value, attribute, unit, time)

    vertices = SkyCoord([10, 30, 50] * u.arcsec, [20, 40, 60] * u.arcsec, obstime=time, observer='earth', frame="helioprojective")
    region = PolygonSkyRegion(vertices=vertices)

    assert region == result
    assert result.vertices.separation(region.vertices).max() < 1 * u.arcsec

def test_parse_chaincode_heliocentric():
    value = "POLYGON((10 20, 30 40, 50 60))"
    attribute = {"frame": "heliocentric"}
    unit = "deg"
    time = Time('2011-08-09 12:04:00.000')

    result = parse_chaincode(value, attribute, unit, time)

    vertices = SkyCoord([10, 30, 50] * u.R_sun, [20, 40, 60] * u.deg, [1, 1, 1] * u.AU, obstime=time, observer='earth', representation_type="cylindrical", frame="heliocentric")
    region = PolygonSkyRegion(vertices=vertices)

    assert region == result
    assert result.vertices.separation(region.vertices).max() < 1 * u.deg

def test_parse_chaincode_icrs():
    value = "POLYGON((10 20, 30 40, 50 60))"
    attribute = {"frame": "icrs"}
    unit = "deg"
    time = Time('2011-08-09 12:04:00.000')

    result = parse_chaincode(value, attribute, unit, time)

    vertices = SkyCoord([10, 30, 50] * u.deg, [20, 40, 60] * u.deg, obstime=time, observer='earth', frame="icrs")
    region = PolygonSkyRegion(vertices=vertices)

    assert region == result
    assert result.vertices.separation(region.vertices).max() < 1 * u.deg

def test_parse_chaincode_points():
    value = "POINT(-107 27)"
    attribute = {"is_point": True, "frame": "heliographic_carrington"}
    unit = "deg"
    time = Time('2011-08-09 12:04:00.000')

    result = parse_chaincode(value, attribute, unit, time)

    center_sky = SkyCoord(-107 * u.deg, 27 * u.deg, obstime=time, observer='earth', frame="heliographic_carrington")
    region = PointSkyRegion(center=center_sky)

    assert region == result
