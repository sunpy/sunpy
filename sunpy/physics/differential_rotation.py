from __future__ import division

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

from sunpy.time import parse_time
from sunpy.coordinates.ephemeris import get_earth
from sunpy.coordinates import frames

__author__ = ["Jose Ivan Campos Rozo", "Stuart Mumford", "Jack Ireland"]
__all__ = ['diff_rot', 'rot_hpc']


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


def solar_rotate_coord_from_earth(start_coordinate, tend, **diff_rot_kwargs):
    """Given a location on the Sun as seen from the Earth, use the solar
    rotation profile to find that location at some later or earlier time.
    Note that this function assumes that the data was observed from the Earth or
    near Earth vicinity.  Specifically, data from SOHO and STEREO observatories
    are not supported.  Note also that the function does NOT use solar B0 and L0
    values provided in the input start co-ordinate -
    these quantities are calculated.

    Parameters
    ----------
    start_coordinate : `~sunpy.coordinates`
        a sunpy co-ordinate

    tend : `sunpy.time.time`
        date/time at which the input co-ordinate will be rotated to.

    **diff_rot_kwargs : keyword arguments
        Keyword arguments are passed on as keyword arguments to `~sunpy.physics.differential_rotation.diff_rot`.

    Returns
    -------
    `~sunpy.coordinates`
        The locations of the input co-ordinates after the application of
         solar rotation in the input co-ordiinate frame.

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import frames
    >>> from sunpy.physics.differential_rotation import solar_rotate_coord_from_earth
    >>> c = SkyCoord(-570 * u.arcsec, 120 * u.arcsec, dateobs='2010-09-10 12:34:56', frame=frame.Helioprojective)
    >>> solar_rotate_coord_from_earth(c, dateobs)
    <SkyCoord (Helioprojective: D0=150634662.59404698 km, dateobs=2010-09-10 13:34:56, L0=0d00m00s, B0=7d14m46.821s, rsun=695508.0 km): (Tx, Ty, distance) in (arcsec, arcsec, km)
    (-562.90765805,  119.31706625,   1.50079871e+08)>

    Notes
    -----
    SSWIDL code equivalent: http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/rot_xy.pro .
    The function rot_xy uses arcmin2hel.pro and hel2arcmin.pro to implement the
    same functionality as this function.  These two functions seem to perform
    inverse operations of each other to a high accuracy.  The corresponding
    equivalent functions here are convert_hpc_hg and convert_hg_hpc
    respectively. These two functions seem to perform inverse
    operations of each other to a high accuracy.  However, the values
    returned by arcmin2hel.pro are slightly different from those provided
    by convert_hpc_hg.  This leads to very slightly different results from
    rot_hpc compared to rot_xy.
    """

    # Make sure we have enough time information to perform a solar differential
    # rotation
    if start_coordinate.dateobs is None:
        raise ValueError('Input co-ordinate(s) must not be of type NoneType')

    # Calculate the interval between the start and end time
    interval = (parse_time(tend) - parse_time(start_coordinate.dateobs)).total_seconds() * u.s

    # Compute Stonyhurst Heliographic co-ordinates - returns (longitude,
    # latitude). Points off the limb are returned as nan.
    heliographic_coordinate = start_coordinate.transform_to('heliographic_stonyhurst')

    # Compute the differential rotation
    drot = diff_rot(interval, heliographic_coordinate.lat.to(u.degree),
                    **diff_rot_kwargs)

    # Rotate the input co-ordinate and update the observer
    heliographic_rotated = SkyCoord(heliographic_coordinate.lon + drot,
                                    heliographic_coordinate.lat,
                                    observer=get_earth(time=tend),
                                    frame=frames.HeliographicStonyhurst)

    # Return the rotated coordinates to the input coordinate frame
    return heliographic_rotated.transform_to(start_coordinate.frame)
