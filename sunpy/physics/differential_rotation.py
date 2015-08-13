from __future__ import division
from datetime import timedelta

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord, Longitude

from sunpy.time import parse_time
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

    return Longitude(rotation_deg * u.deg)


def solar_rotate_coordinate(coordinate, new_observer_time, new_observer_location="earth", **diff_rot_kwargs):
    """
    Given a coordinate on the Sun, calculate where that coordinate maps to
    at some later or earlier time, given the solar rotation profile.  Note that
    if the new observer location is defined using a BaseCoordinateFrame or
    SkyCoord, then it is assumed that the new observer location is correct for
    the new observer time that was also passed in.
    Parameters
    ----------
    coordinate : `~astropy.coordinates.SkyCoord`
        Any valid coordinate which is transformable to Heliographic Stonyhurst.
    new_observer_time : sunpy-compatible time
        date/time at which the input co-ordinate will be rotated to.
    new_observer_location : `str`, `~astropy.coordinates.BaseCoordinateFrame`, `~astropy.coordinates.SkyCoord`
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
    `~astropy.coordinates.SkyCoord``
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


def _to_norm(arr):
    '''
    Helpper function to normalise/scale an array.  This is needed for example
    for scikit-image which uses flotas between 0 and 1

    Parameters
    ----------
    arr : `~numpy.ndarray`
        Array of floats whised to normalise

    Returns
    -------
    arr : `~numpy.ndarray`
        Array with values between 0 (min) and 1 (max)

    Examples
    --------
    >>> import numpy as np
    >>> from sunpy.physics.transforms.differential_rotation import _to_norm
    >>> out = _to_norm(np.array([-1, 0, 1]))
    >>> out
    array([ 0. ,  0.5,  1. ])
    '''
    from skimage.util import img_as_float
    arr = np.array(arr, dtype='double')
    arr = img_as_float(arr, force_copy=True)
    if arr.min() < 0:
        arr += np.abs(arr.min())
    arr /= arr.max()
    return arr

def _un_norm(arr, original):
    '''
    Helpper function tu Un-normalises (or re-scale) an array based in
    the values of the original array.

    Parameters
    ----------
    arr : `~numpy.ndarray`
        Array of floats whised to un-normalise with values in [0,1]
    original : `~numpy.ndarray`
        Original array with the min and max values

    Returns
    -------
    arr : `~numpy.ndarray`
        Array with values between `original.min()` and `original.max()`

    Examples
    --------
    >>> import numpy as np
    >>> from sunpy.physics.transforms.differential_rotation import _un_norm
    >>> original = np.array([-1, 0, 1])
    >>> normalised = np.array([0., 0.5, 1.])
    >>> out = _un_norm(normalised, original)
    >>> out
    array([-1.,  0.,  1.])
    '''
    level = 0 if original.min() > 0 else np.abs(original.min())
    arr *= original.max() + level
    arr -= level
    return arr

def _warp_sun(xy, smap, timedelta):
    '''
    Function that returns a new list of coordinates for each input coord.
    This is an inverse function needed by the scikit-image `transform.warp`
    function.

    Parameters
    ----------
    xy :
        Array from `transform.warp`
    smap : `~sumpy.map`
        Original map that we want to transform
    timedelta : `~datetime.timedelta`
        Desired interval for the differential rotation

    Returns
    -------
    xy2 : `~numpy.ndarray`
        Array with the inverse transformation
    '''
    #Calculate the hpc coords
    x = np.arange(0, smap.dimensions.x.value)
    y = np.arange(0, smap.dimensions.y.value)
    xx, yy = np.meshgrid(x, y)
    hpc_coords = smap.pixel_to_data(xx * u.pix, yy * u.pix)

    #Do the diff rot
    rotted = rot_hpc(hpc_coords[1], hpc_coords[0], smap.date + timedelta,
                     smap.date, occultation=True)
    # `transform.warp` needs the inverse rotation, therefore it's needed
    # to provide the transform from the desired date to the original date

    #Go back to pixel coords
    x2,y2 = smap.data_to_pixel(rotted[0], rotted[1])

    #Restack the data to make it correct output form
    xy2 = np.column_stack([x2.value.flat,y2.value.flat])

    #Remove NaNs
    mask = np.isnan(xy2)
    xy2[mask] = 0.0

    return xy2

def difrot_map(smap, timedelta):
    '''
    Function to compensate for differential rotation a sunpy map

    Parameters
    ----------
    smap : `~sumpy.map`
        Original map that we want to transform
    timedelta : `~datetime.timedelta`
        Desired interval for the differential rotation

    Returns
    -------
    rotmap : `~sunpy.map`
        Array with the inverse transformation
    '''
    from skimage import transform

    out = transform.warp(_to_norm(smap.data), inverse_map=_warp_sun,
                         map_args={'smap':smap, 'timedelta': timedelta })
    out = _un_norm(out, smap.data)
    rotmap = sunpy.map.Map((out, smap.meta))
    rotmap_rotdate = rotmap.date + timedelta # FixMe update the map date (maybe a new keyword in Map?)
    return rotmap
