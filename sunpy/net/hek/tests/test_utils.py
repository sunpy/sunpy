import json

import numpy as np
import pytest
from regions import PolygonSkyRegion

from astropy import units as u
from astropy.coordinates import SkyCoord

from sunpy.net import attrs, hek
from sunpy.net.hek.tests.test_hek import HEKResult
from sunpy.net.hek.utils import COORD_FILE_PATH, UNIT_FILE_PATH, get_unit, parse_chaincode

_hek_result = HEKResult()

@pytest.fixture
def read_unit_attributes():
    with open(UNIT_FILE_PATH) as unit_file:
        unit_properties = json.load(unit_file)

    return unit_properties

@pytest.fixture
def read_coord_attributes():
    with open(COORD_FILE_PATH) as coord_file:
        coord_properties = json.load(coord_file)

    return coord_properties

@pytest.mark.remote_data
def test_missing_times():
    # Check for https://github.com/sunpy/sunpy/pull/7627#issuecomment-2113451964
    client = hek.HEKClient()
    results = client.search(attrs.Time('2024-05-10', '2024-05-12'), attrs.hek.AR.NOAANum == 13664)
    assert isinstance(results["event_peaktime"][0], np.ma.core.MaskedConstant)
    assert results["event_peaktime"][6].isot == "2024-05-10T16:08:00.000"

@pytest.mark.remote_data
def test_astropy_unit_parsing(read_unit_attributes, read_coord_attributes):
    client = hek.HEKClient()
    tstart = '2014/10/24 20:50'
    tend = '2014/10/25 00:14'
    event_type = 'FL'
    result = client.search(attrs.Time(tstart, tend), attrs.hek.EventType(event_type))
    unit_properties = read_unit_attributes
    coord_properties = read_coord_attributes
    unit_attributes_with_unit = [ prop for prop in unit_properties["attributes"] if prop.get("unit_prop",None) is not None]
    coord_attributes_with_unit = [prop for prop in coord_properties["attributes"] if not prop.get("is_chaincode", False) and not prop.get("is_unit_prop",False)]

    for attribute in unit_attributes_with_unit + coord_attributes_with_unit:
        if attribute["name"] in result.colnames:
            assert all([value in ['', None] or isinstance(value, u.Quantity) for value in result[attribute['name']]])


@pytest.mark.remote_data
def test_chaincode_parsing(read_coord_attributes):
    client = hek.HEKClient()
    tstart = '2014/10/24 20:50'
    tend = '2014/10/25 00:14'
    event_type = 'FL'
    result = client.search(attrs.Time(tstart, tend), attrs.hek.EventType(event_type))
    coord_properties = read_coord_attributes
    chaincode_properties = [prop for prop in coord_properties["attributes"] if prop.get("is_chaincode", False)]

    for attribute in chaincode_properties:
        if attribute["name"] in result.colnames:
            assert all([value in ['', None] or isinstance(value, PolygonSkyRegion) for value in result[attribute['name']]])

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
