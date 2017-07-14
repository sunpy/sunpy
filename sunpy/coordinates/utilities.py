# -*- coding: utf-8 -*-
from __future__ import absolute_import, division
import datetime

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body_barycentric, PrecessedGeocentric, Angle
from astropy.coordinates.representation import CartesianRepresentation

from sunpy.time import parse_time

from .frames import HeliographicStonyhurst as HGS

__all__ = ['get_earth', 'get_sun_B0', 'get_sun_P']


def get_earth(time='now'):
    """
    Return a SkyCoord for the location of the Earth at a specified time in the
    HeliographicStonyhurst frame.  The longitude will be 0 by definition.

    Parameters
    ----------
    time : various
        Time to use in a parse_time-compatible format

    Returns
    -------
    out : `~astropy.coordinates.SkyCoord`
        SkyCoord for the location of the Earth in the HeliographicStonyhurst frame
    """
    obstime = Time(parse_time(time))

    earth_icrs = get_body_barycentric('earth', obstime)
    earth = SkyCoord(earth_icrs, frame='icrs', obstime=obstime).transform_to(HGS)

    # Explicitly set the longitude to 0
    earth = SkyCoord(0*u.deg, earth.lat, earth.radius, frame=earth)

    return earth


def get_sun_B0(time='now'):
    """
    Return the B0 angle for the Sun at a specified time, which is the heliographic latitude of the
    Sun-disk center as seen from Earth.  The range of B0 is +/-7.23 degrees.

    Parameters
    ----------
    time : various
        Time to use in a parse_time-compatible format

    Returns
    -------
    out : `~astropy.coordinates.Angle`
        The position angle
    """
    return Angle(get_earth(time).lat)


def get_sun_P(time='now'):
    """
    Return the position (P) angle for the Sun at a specified time, which is the angle between
    geocentric north and solar north as seen from Earth, measured eastward from geocentric north.
    The range of P is +/-26.3 degrees.

    Parameters
    ----------
    time : various
        Time to use in a parse_time-compatible format

    Returns
    -------
    out : `~astropy.coordinates.Angle`
        The position angle
    """
    obstime = Time(parse_time(time))

    sun_center = SkyCoord(0*u.deg, 0*u.deg, 0*u.km, frame=HGS, obstime=obstime)
    sun_north_pole = SkyCoord(0*u.deg, 90*u.deg, 690000*u.km, frame=HGS, obstime=obstime)

    # Represent the two north vectors in precessed geocentric coordinates
    geocentric = PrecessedGeocentric(equinox=obstime, obstime=obstime)
    sky_normal_vector = sun_center.transform_to(geocentric).data.to_cartesian()
    sun_north_pole_vector = sun_north_pole.transform_to(geocentric).data.to_cartesian()
    earth_north_pole_vector = CartesianRepresentation(0, 0, 1)

    # Use cross products to obtain the sky projections of the two north vectors
    sun_in_sky = sun_north_pole_vector.cross(sky_normal_vector)
    earth_in_sky = earth_north_pole_vector.cross(sky_normal_vector)

    # Use a cross product to calculate the signed angle between the two projected north vectors
    cross_vector = sun_in_sky.cross(earth_in_sky) / (sun_in_sky.norm() * earth_in_sky.norm())
    angle = np.arcsin(cross_vector.dot(sky_normal_vector) / sky_normal_vector.norm()).to('deg')

    return Angle(angle)
