import numpy as np
import pytest
from hypothesis import given, settings

import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.tests.helper import assert_quantity_allclose

from sunpy.coordinates import NorthOffsetFrame
from sunpy.coordinates.tests.helpers import assert_longitude_allclose
from sunpy.coordinates.tests.strategies import latitudes, longitudes


def test_null():
    """
    test init of a frame where the origins are the same.
    """
    off = NorthOffsetFrame(north=SkyCoord(0*u.deg, 90*u.deg,
                                          frame='heliographic_stonyhurst'))
    assert isinstance(off, SkyOffsetFrame)
    assert off.origin.lat == 0*u.deg
    assert off.origin.lon == 0*u.deg


@given(lon=longitudes(), lat=latitudes())
@settings(deadline=1000)
def test_transform(lon, lat):
    """
    Test that the north pole in the new frame transforms back to the given
    north argument.
    """
    north = SkyCoord(lon=lon, lat=lat, frame='heliographic_stonyhurst')
    off = NorthOffsetFrame(north=north)
    t_north = SkyCoord(lon=0*u.deg, lat=90*u.deg, frame=off)
    t_north = t_north.transform_to('heliographic_stonyhurst')
    assert_quantity_allclose(north.separation(t_north), 0*u.deg, atol=1e-6*u.deg)


def test_south_pole():
    s = SkyCoord(-10*u.deg, 0*u.deg, frame='heliographic_stonyhurst')
    off = NorthOffsetFrame(north=s)
    assert_longitude_allclose(off.origin.lon, 170*u.deg)
    assert_quantity_allclose(off.origin.lat, -90*u.deg)


def test_cartesian():
    # Test a north pole specified in Cartesian coordinates
    s = SkyCoord(1, 2, 3, unit=u.m,
                 representation_type='cartesian', frame='heliographic_stonyhurst')
    off = NorthOffsetFrame(north=s)
    assert_longitude_allclose(off.origin.lon, np.arctan2(s.y, s.x))
    assert_quantity_allclose(off.origin.lat, np.arctan2(s.z, np.sqrt(s.x**2 + s.y**2)) - 90*u.deg)


def test_error():
    with pytest.raises(TypeError):
        NorthOffsetFrame()
