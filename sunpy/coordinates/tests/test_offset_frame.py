import pytest
import numpy as np

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.coordinates import SkyCoord, SkyOffsetFrame

from sunpy.coordinates import NorthOffsetFrame


def test_null():
    """
    test init of a frame where the origins are the same.
    """
    off = NorthOffsetFrame(north=SkyCoord(0*u.deg, 90*u.deg,
                                          frame='heliographic_stonyhurst'))
    assert isinstance(off, SkyOffsetFrame)
    assert off.origin.lat == 0*u.deg
    assert off.origin.lon == 0*u.deg


def test_transform(lon=0*u.deg, lat=0*u.deg):
    """
    Test that the north pole in the new frame transforms back to the given
    north argument.
    """
    north = SkyCoord(lon=lon, lat=lat, frame='heliographic_stonyhurst')
    off = NorthOffsetFrame(north=north)
    t_north = SkyCoord(lon=0*u.deg, lat=90*u.deg, frame=off)
    t_north = t_north.transform_to('heliographic_stonyhurst')
    assert_quantity_allclose(north.lon, t_north.lon, rtol=1e6)
    assert_quantity_allclose(north.lat, t_north.lat, rtol=1e6)

