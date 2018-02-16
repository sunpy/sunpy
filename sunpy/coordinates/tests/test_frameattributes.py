# -*- coding: utf-8 -*-

import datetime

import pytest

import astropy.units as u
from astropy.time import Time
from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import ICRS, get_body_barycentric

from sunpy.time import parse_time
from ..frames import Helioprojective, HeliographicStonyhurst
from ..frameattributes import TimeFrameAttributeSunPy, ObserverCoordinateAttribute
from sunpy.coordinates import get_earth, frames


@pytest.fixture
def attr():
    return TimeFrameAttributeSunPy()


def test_now(attr):
    """ We can't actually test the value independantly """
    result, converted = attr.convert_input('now')

    assert isinstance(result, Time)
    assert converted


def test_none(attr):
    """ We can't actually test the value independantly """
    result, converted = attr.convert_input(None)

    assert result is None
    assert not converted


@pytest.mark.parametrize('input', [
    Time('2012-01-01 00:00:00'), '2012/01/01T00:00:00', '20120101000000', '2012/01/01 00:00:00'
])
def test_convert(attr, input):
    result, converted = attr.convert_input(input)

    output = Time('2012-01-01 00:00:00')

    assert isinstance(result, Time)
    assert result == output


@pytest.mark.parametrize('input', [
    Time('2012-01-01 00:00:00'), '2012/01/01T00:00:00', '20120101000000', '2012/01/01 00:00:00'
])
def test_on_frame(input):
    hpc1 = Helioprojective(obstime=input)

    output = Time('2012-01-01 00:00:00')

    assert isinstance(hpc1.obstime, Time)
    assert hpc1.obstime == output


def test_non_string():
    output = datetime.datetime.now()

    hpc1 = Helioprojective(obstime=output)

    assert isinstance(hpc1.obstime, Time)
    assert hpc1.obstime == output


def test_on_frame_error():
    with pytest.raises(ValueError):
        Helioprojective(obstime='ajshdasjdhk')


def test_on_frame_error2():
    with pytest.raises(ValueError):
        Helioprojective(obstime=17263871263)


# ObserverCoordinateAttribute


def test_string_coord():

    oca = ObserverCoordinateAttribute(HeliographicStonyhurst)

    obstime = "2011-01-01"
    coord = oca._convert_string_to_coord("earth", obstime)

    assert isinstance(coord, HeliographicStonyhurst)

    assert coord.obstime == parse_time(obstime)


def test_coord_get():

    # Test default (instance=None)
    obs = Helioprojective.observer
    assert obs is "earth"

    # Test get
    obstime = "2013-04-01"
    obs = Helioprojective(observer="earth", obstime=obstime).observer
    earth = get_earth(obstime)
    assert isinstance(obs, HeliographicStonyhurst)
    assert_quantity_allclose(obs.lon, earth.lon)
    assert_quantity_allclose(obs.lat, earth.lat)
    assert_quantity_allclose(obs.radius, earth.radius)

    # Test get
    obstime = "2013-04-01"
    obs = Helioprojective(obstime=obstime).observer
    earth = get_earth(obstime)
    assert isinstance(obs, HeliographicStonyhurst)
    assert_quantity_allclose(obs.lon, earth.lon)
    assert_quantity_allclose(obs.lat, earth.lat)
    assert_quantity_allclose(obs.radius, earth.radius)

    # Test get mars
    obstime = Time(parse_time("2013-04-01"))
    obs = Helioprojective(observer="mars", obstime=obstime).observer
    out_icrs = ICRS(get_body_barycentric("mars", obstime))
    mars = out_icrs.transform_to(HeliographicStonyhurst(obstime=obstime))

    assert isinstance(obs, HeliographicStonyhurst)
    assert_quantity_allclose(obs.lon, mars.lon)
    assert_quantity_allclose(obs.lat, mars.lat)
    assert_quantity_allclose(obs.radius, mars.radius)


def test_default_hcc_observer():
    h = frames.Heliocentric()
    assert h.observer is "earth"

    h = frames.Heliocentric(observer="mars")
    assert h.observer is "mars"


def test_obstime_hack():
    """
    Test that the obstime can be updated in place, this is used in the transform pipeline.
    """
    h = frames.Heliocentric()
    assert h.observer is "earth"

    obstime = "2011-01-01"
    h._obstime = obstime

    assert isinstance(h.observer, frames.HeliographicStonyhurst)

    earth = get_earth(obstime)
    obs = h._observer
    assert isinstance(obs, HeliographicStonyhurst)
    assert_quantity_allclose(obs.lon, earth.lon)
    assert_quantity_allclose(obs.lat, earth.lat)
    assert_quantity_allclose(obs.radius, earth.radius)


"""
These two tests are to make sure that during the transformation stack the value
of observer is correctly calculated.
"""


def test_default_observer_transform_hcc():
    center = frames.HeliographicStonyhurst(0 * u.deg, 0 * u.deg, obstime="2017-07-11 15:00")
    hpc = center.transform_to(frames.Heliocentric(obstime="2017-07-11 15:00"))

    assert_quantity_allclose(hpc.y, -48471.1283979 * u.km)


def test_default_observer_transform_hpc():
    center = frames.HeliographicStonyhurst(0 * u.deg, 0 * u.deg, obstime="2017-07-11 15:00")
    hpc = center.transform_to(frames.Helioprojective(obstime="2017-07-11 15:00"))

    assert_quantity_allclose(hpc.Ty, -66.04425197 * u.arcsec)
