
import pytest

import astropy.units as u
from astropy.coordinates import ICRS, HeliocentricMeanEcliptic, get_body_barycentric
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from sunpy.coordinates import frames, get_earth
from sunpy.time import parse_time
from ..frameattributes import ObserverCoordinateAttribute, TimeFrameAttributeSunPy
from ..frames import HeliocentricInertial, HeliographicCarrington, HeliographicStonyhurst, Helioprojective


@pytest.fixture
def attr():
    return TimeFrameAttributeSunPy()


@pytest.fixture
def oca():
    return ObserverCoordinateAttribute(HeliographicStonyhurst)


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
    output = parse_time('now')

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


def test_string_coord(oca):
    obstime = "2011-01-01"
    coord = oca._convert_string_to_coord("earth", obstime)

    assert isinstance(coord, HeliographicStonyhurst)

    assert coord.obstime == parse_time(obstime)


def test_observer_not_hgs_sunpy(oca):
    observer = HeliocentricInertial(0*u.deg, 0*u.deg, 1*u.AU, obstime='2001-01-01')
    result, converted = oca.convert_input(observer)

    assert isinstance(result, HeliographicStonyhurst)
    assert result.obstime == observer.obstime
    assert converted


def test_observer_not_hgs_astropy(oca):
    observer = HeliocentricMeanEcliptic(0*u.deg, 0*u.deg, 1*u.AU, obstime='2001-01-01')
    result, converted = oca.convert_input(observer)

    assert isinstance(result, HeliographicStonyhurst)
    assert result.obstime == observer.obstime
    assert converted


def test_coord_get():

    # Test default (instance=None)
    obs = Helioprojective.observer
    assert obs is None

    # Test get
    obstime = "2013-04-01"
    obs = Helioprojective(observer="earth", obstime=obstime).observer
    earth = get_earth(obstime)
    assert isinstance(obs, HeliographicStonyhurst)
    assert_quantity_allclose(obs.lon, earth.lon)
    assert_quantity_allclose(obs.lat, earth.lat)
    assert_quantity_allclose(obs.radius, earth.radius)

    assert hasattr(obs, "object_name")
    assert obs.object_name == "earth"
    assert str(obs) == "<HeliographicStonyhurst Coordinate for 'earth'>"

    # Test get
    obstime = "2013-04-01"
    obs = Helioprojective(observer="earth", obstime=obstime).observer
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

    assert hasattr(obs, "object_name")
    assert obs.object_name == "mars"
    assert str(obs) == "<HeliographicStonyhurst Coordinate for 'mars'>"


def test_observer_self_get():
    hgc_noobstime = HeliographicCarrington(observer="self")
    assert hgc_noobstime.observer == "self"

    hgc_obstime = HeliographicCarrington(observer="self", obstime="2001-01-01")
    assert hgc_obstime.observer == "self"


def test_default_hcc_observer():
    h = frames.Heliocentric()
    assert h.observer is None

    h = frames.Heliocentric(observer="mars")
    assert h.observer == "mars"


def test_obstime_hack():
    """
    Test that the obstime can be updated in place, this is used in the transform pipeline.
    """
    h = frames.Heliocentric(observer="earth")

    obstime = "2011-01-01"
    h._obstime = obstime

    assert isinstance(h.observer, frames.HeliographicStonyhurst)

    earth = get_earth(obstime)
    obs = h._observer
    assert isinstance(obs, HeliographicStonyhurst)
    assert_quantity_allclose(obs.lon, earth.lon)
    assert_quantity_allclose(obs.lat, earth.lat)
    assert_quantity_allclose(obs.radius, earth.radius)
