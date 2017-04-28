import pytest
import numpy as np
import hypothesis.strategies as st
from hypothesis import given

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import SkyCoord, SkyOffsetFrame

from sunpy.coordinates import NorthOffsetFrame


@st.composite
def latitude(draw, lat=st.floats(min_value=-90, max_value=90,
                                 allow_nan=False, allow_infinity=False)):
    return draw(lat) * u.deg


@st.composite
def lonitude(draw, lon=st.floats(min_value=-180, max_value=180,
                                 allow_nan=False, allow_infinity=False)):
    return draw(lon) * u.deg


def test_null():
    """
    test init of a frame where the origins are the same.
    """
    off = NorthOffsetFrame(north=SkyCoord(0*u.deg, 90*u.deg,
                                          frame='heliographic_stonyhurst'))
    assert isinstance(off, SkyOffsetFrame)
    assert off.origin.lat == 0*u.deg
    assert off.origin.lon == 0*u.deg


@given(lon=lonitude(), lat=latitude())
def test_transform(lon, lat):
    """
    Test that the north pole in the new frame transforms back to the given
    north argument.
    """
    north = SkyCoord(lon=lon, lat=lat, frame='heliographic_stonyhurst')
    off = NorthOffsetFrame(north=north)
    t_north = SkyCoord(lon=0*u.deg, lat=90*u.deg, frame=off)
    t_north = t_north.transform_to('heliographic_stonyhurst')
    assert_quantity_allclose(north.lon, t_north.lon, atol=1e6*u.deg)
    assert_quantity_allclose(north.lat, t_north.lat, atol=1e6*u.deg)


def test_south_pole():
    s = SkyCoord(-10*u.deg, 0*u.deg, frame='heliographic_stonyhurst')
    off = NorthOffsetFrame(north=s)
    assert_quantity_allclose(off.origin.lon, 170*u.deg)
    assert_quantity_allclose(off.origin.lat, -90*u.deg)


def test_error():
    with pytest.raises(TypeError):
        NorthOffsetFrame()
