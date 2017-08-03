from __future__ import division

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates.baseframe import BaseCoordinateFrame

from sunpy.extern import six
from sunpy.time import parse_time
from sunpy.coordinates.ephemeris import get_earth, get_body_heliographic_stonyhurst
from sunpy.coordinates import frames

__author__ = ["Jose Ivan Campos Rozo", "Stuart Mumford", "Jack Ireland"]
__all__ = ['diff_rot', 'solar_rotate_coordinate']


@u.quantity_input(duration=u.s, latitude=u.degree)
def diff_rot(duration, latitude, rot_type='howard', frame_time='sidereal'):
    """
    This function computes the change in longitude over days in degrees.

    Parameters
    -----------
    duration : `~astropy.units.Quantity`
        Number of seconds to rotate over.
    latitude : `~astropy.units.Quantity`
        heliographic coordinate latitude in Degrees.
    rot_type : {'howard' | 'snodgrass' | 'allen'}
        howard : Use values for small magnetic features from Howard et al.
        snodgrass : Use Values from Snodgrass et. al
        allen : Use values from Allen, Astrophysical Quantities, and simpler equation.
    frame_time : {'sidereal' | 'synodic'}
        Choose 'type of day' time reference frame.

    Returns
    -------
    longitude_delta : `~astropy.units.Quantity`
        The change in longitude over days (units=degrees)

    Notes
    -----
    * IDL code equivalent: http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/diff_rot.pro
    * Howard rotation: http://adsabs.harvard.edu/abs/1990SoPh..130..295H
    * A review of rotation parameters (including Snodgrass values): http://link.springer.com/article/10.1023%2FA%3A1005226402796

    Examples
    --------
    Default rotation calculation over two days at 30 degrees latitude:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from sunpy.physics.differential_rotation import diff_rot
    >>> rotation = diff_rot(2 * u.day, 30 * u.deg)

    Default rotation over two days for a number of latitudes:

    >>> rotation = diff_rot(2 * u.day, np.linspace(-70, 70, 20) * u.deg)

    With rotation type 'allen':

    >>> rotation = diff_rot(2 * u.day, np.linspace(-70, 70, 20) * u.deg, 'allen')
    """

    latitude = latitude.to(u.deg)
    delta_seconds = duration.to(u.s).value
    delta_days = delta_seconds / 24.0 / 3600.0

    sin2l = (np.sin(latitude))**2
    sin4l = sin2l**2

    rot_params = {'howard': [2.894, -0.428, -0.370],
                  'snodgrass': [2.851, -0.343, -0.474]
                  }

    if rot_type not in ['howard', 'allen', 'snodgrass']:
        raise ValueError("""rot_type must equal one of
                        { 'howard' | 'allen' | 'snodgrass' }""")

    elif rot_type == 'allen':
        rotation_deg = delta_days * (14.44 - (3.0 * sin2l))

    else:
        A, B, C = rot_params[rot_type]

        # This is in micro-radians / sec
        rotation_rate = A + B * sin2l + C * sin4l
        rotation_deg = rotation_rate * 1e-6 * delta_seconds / np.deg2rad(1)

    if frame_time == 'synodic':
        rotation_deg -= 0.9856 * delta_days

    return rotation_deg * u.deg


def solar_rotate_coordinate(coordinate,
                            new_observer_time,
                            new_observer_location="earth", **diff_rot_kwargs):
    """Given a coordinate on the Sun, calculate where that coordinate maps to
    at some later or earlier time, given the solar rotation profile.  Note that
    if the new observer location is defined using a BaseCoordinateFrame or
    SkyCoord, then it is assumed that the new observer location is correct for
    the new observer time that was also passed in.

    Parameters
    ----------
    coordinate : `~sunpy.coordinates`
        a sunpy coordinate

    new_observer_time : `sunpy.time.time`
        date/time at which the input co-ordinate will be rotated to.

    new_observer_location : None, str, BaseCoordinateFrame, SkyCoord
        The solar-system body for which to calculate observer locations.  Note
        that spacecraft are not explicitly supported as yet.  Instruments in
        Earth orbit can be approximated by using the default setting.  If a
        BaseCoordinateFrame or SkyCoord are passed in, it is assumed that
        this observer location is correct for the observer time that was also
        passed in.

    **diff_rot_kwargs : keyword arguments
        Keyword arguments are passed on as keyword arguments to `~sunpy.physics.differential_rotation.diff_rot`.

    Returns
    -------
    `~sunpy.coordinates`
        The locations of the input coordinates after the application of
        solar rotation in the input coordinate frame.

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import frames
    >>> from sunpy.physics.differential_rotation import solar_rotate_coordinate
    >>> from sunpy.coordinates.ephemeris import get_earth
    >>> obstime = '2010-09-10 12:34:56'
    >>> c = SkyCoord(-570*u.arcsec, 120*u.arcsec, obstime=obstime, observer=get_earth(obstime), frame=frames.Helioprojective)
    >>> solar_rotate_coordinate(c, '2010-09-10 13:34:56')
    <SkyCoord (Helioprojective: obstime=2010-09-10 13:34:56, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2010-09-10 13:34:56): (lon, lat, radius) in (deg, deg, AU)
    ( 0.,  7.24822784,  1.00695436)>): (Tx, Ty, distance) in (arcsec, arcsec, km)
    (-562.37689548,  119.26840368,   1.50083152e+08)>
    """

    # Make sure we have enough time information to perform a solar differential
    # rotation.
    if coordinate.obstime is None:
        raise ValueError('Input coordinate(s) obstime must not be of type NoneType')

    # Calculate the interval between the start and end time
    interval = (parse_time(new_observer_time) - parse_time(coordinate.obstime)).total_seconds() * u.s

    # Compute Stonyhurst Heliographic co-ordinates - returns (longitude,
    # latitude). Points off the limb are returned as nan.
    heliographic_coordinate = coordinate.transform_to('heliographic_stonyhurst')

    # Compute the differential rotation
    drot = diff_rot(interval, heliographic_coordinate.lat.to(u.degree), **diff_rot_kwargs)

    # Rotate the input co-ordinate and update the observer
    heliographic_rotated = SkyCoord(heliographic_coordinate.lon + drot,
                                    heliographic_coordinate.lat,
                                    obstime=new_observer_time,
                                    observer=new_observer_location,
                                    frame=frames.HeliographicStonyhurst)

    # Return the rotated coordinates to the input coordinate frame
    return heliographic_rotated.transform_to(coordinate.frame.name)
