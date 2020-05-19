"""
Provide a set of hypothesis strategies for various coordinates-related tests.
"""
import hypothesis.strategies as st

import astropy.units as u
from astropy.coordinates import Latitude, Longitude

from sunpy.time import parse_time


@u.quantity_input
@st.composite
def latitudes(draw, min_lat: u.deg = -90*u.deg, max_lat: u.deg = 90*u.deg):
    lat = st.floats(min_value=min_lat.to_value(u.deg),
                    max_value=max_lat.to_value(u.deg),
                    allow_nan=False, allow_infinity=False)
    return Latitude(draw(lat) * u.deg)


@u.quantity_input
@st.composite
def longitudes(draw, min_lon: u.deg = -180*u.deg, max_lon: u.deg = 180*u.deg,
               wrap_angle: u.deg = 180*u.deg):
    lon = st.floats(min_value=min_lon.to_value(u.deg),
                    max_value=max_lon.to_value(u.deg),
                    allow_nan=False, allow_infinity=False)
    return Longitude(draw(lon) * u.deg, wrap_angle=wrap_angle)


@st.composite
def times(draw, min_time='1960-01-01', max_time='2024-01-01'):
    days = st.floats(min_value=0,
                     max_value=(parse_time(max_time) - parse_time(min_time)).to(u.day).value,
                     allow_nan=False, allow_infinity=False)
    return parse_time(min_time) + draw(days) * u.day
