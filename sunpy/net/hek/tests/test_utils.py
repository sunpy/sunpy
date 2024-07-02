from regions import PointSkyRegion, PolygonSkyRegion

from astropy import units as u
from astropy.coordinates import SkyCoord

from sunpy.net import Fido
from sunpy.net import attrs as a
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

    result = parse_chaincode(value, attribute, unit)

    vertices = SkyCoord([10, 30, 50] * u.arcsec, [20, 40, 60] * u.arcsec, frame="helioprojective")
    region = PolygonSkyRegion(vertices=vertices)

    assert region == result
    assert result.vertices.separation(region.vertices).max() < 1 * u.arcsec

def test_parse_chaincode_heliocentric():
    value = "POLYGON((10 20, 30 40, 50 60))"
    attribute = {"frame": "heliocentric"}
    unit = "deg"

    result = parse_chaincode(value, attribute, unit)

    vertices = SkyCoord([10, 30, 50] * u.R_sun, [20, 40, 60] * u.deg, [1, 1, 1] * u.AU, representation_type="cylindrical", frame="heliocentric")
    region = PolygonSkyRegion(vertices=vertices)

    assert region == result
    assert result.vertices.separation(region.vertices).max() < 1 * u.deg

def test_parse_chaincode_icrs():
    value = "POLYGON((10 20, 30 40, 50 60))"
    attribute = {"frame": "icrs"}
    unit = "deg"

    result = parse_chaincode(value, attribute, unit)

    vertices = SkyCoord([10, 30, 50] * u.deg, [20, 40, 60] * u.deg, frame="icrs")
    region = PolygonSkyRegion(vertices=vertices)

    assert region == result
    assert result.vertices.separation(region.vertices).max() < 1 * u.deg

def test_parse_chaincode_points():
    value = "POINT(-107 27)"
    attribute = {"is_point": True, "frame": "heliographic_carrington"}
    unit = "deg"

    result = parse_chaincode(value, attribute, unit)

    center_sky = SkyCoord(-107 * u.deg, 27 * u.deg, frame="heliographic_carrington")
    region = PointSkyRegion(center=center_sky)

    assert region == result

def test_merging_event_coords():
    tstart = '2011/08/09 07:23:56'
    tend = '2011/08/09 12:40:29'
    result = Fido.search(a.Time(tstart,tend), a.hek.EventType('CH'))

    coord1 = -2.91584*u.arcsec
    coord2 = 940.667*u.arcsec
    frame='helioprojective'
    event_coord = SkyCoord(coord1, coord2, frame=frame)

    assert result[0]['event_coord'][0] == event_coord
