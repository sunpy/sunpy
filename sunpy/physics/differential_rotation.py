from __future__ import division
import datetime
from copy import deepcopy
import warnings
from itertools import product

import numpy as np
from skimage import transform
from astropy import units as u
from astropy.coordinates import SkyCoord, Longitude

import sunpy.map
from sunpy.time import parse_time
from sunpy.coordinates import frames, HeliographicStonyhurst
from sunpy.image.util import to_norm, un_norm

__all__ = ['diff_rot', 'solar_rotate_coordinate', 'diffrot_map']


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
    rot_type : `str`
        The differential rotation model to use.

        One of:

        | ``howard`` : Use values for small magnetic features from Howard et al.
        | ``snodgrass`` : Use Values from Snodgrass et. al
        | ``allen`` : Use values from Allen's Astrophysical Quantities, and simpler equation.

    frame_time : `str`
        One of : ``'sidereal'`` or  ``'synodic'``. Choose 'type of day' time reference frame.

    Returns
    -------
    longitude_delta : `~astropy.units.Quantity`
        The change in longitude over days (units=degrees)

    References
    ----------

    * `IDL code equivalent <https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/diff_rot.pro>`__
    * `Howard rotation <http://adsabs.harvard.edu/abs/1990SoPh..130..295H>`__
    * `A review of rotation parameters (including Snodgrass values) <https://doi.org/10.1023/A:1005226402796>`__

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

    sin2l = (np.sin(latitude))**2
    sin4l = sin2l**2

    rot_params = {'howard': [2.894, -0.428, -0.370] * u.urad / u.second,
                  'snodgrass': [2.851, -0.343, -0.474] * u.urad / u.second,
                  'allen': [14.44, -3.0, 0] * u.deg / u.day
                  }

    if rot_type not in ['howard', 'allen', 'snodgrass']:
        raise ValueError(("rot_type must equal one of "
                          "{{ {} }}".format(" | ".join(rot_params.keys()))))

    A, B, C = rot_params[rot_type]

    rotation = (A + B * sin2l + C * sin4l) * duration

    if frame_time == 'synodic':
        rotation -= 0.9856 * u.deg / u.day * duration

    return Longitude(rotation.to(u.deg))


def solar_rotate_coordinate(coordinate,
                            new_observer_time,
                            new_observer_location="earth",
                            **diff_rot_kwargs):
    """
    Given a coordinate on the Sun, calculate where that coordinate maps to
    at some later or earlier time, given the solar rotation profile.

    Note that if the new observer location is defined using a
    BaseCoordinateFrame or SkyCoord, then it is assumed that the new observer
    location is correct for the new observer time that was also passed in.

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
    coordinate : `~astropy.coordinates.SkyCoord``
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
        (0., 7.24822784, 1.00695436)>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (-562.37689548, 119.26840368, 1.50083152e+08)>

    """

    # Calculate the interval between the start and end time
    interval = (
        parse_time(new_observer_time) - parse_time(coordinate.obstime)).total_seconds() * u.s

    # Compute Stonyhurst Heliographic co-ordinates - returns (longitude,
    # latitude). Points off the limb are returned as nan.
    heliographic_coordinate = coordinate.transform_to('heliographic_stonyhurst')

    # Compute the differential rotation
    drot = diff_rot(interval, heliographic_coordinate.lat.to(u.degree), **diff_rot_kwargs)

    # Rotate the input co-ordinate and update the observer
    heliographic_rotated = SkyCoord(
        heliographic_coordinate.lon + drot,
        heliographic_coordinate.lat,
        obstime=new_observer_time,
        observer=new_observer_location,
        frame=frames.HeliographicStonyhurst)

    # Return the rotated coordinates to the input coordinate frame
    return heliographic_rotated.transform_to(coordinate.frame.name)


@u.quantity_input(dt=u.s)
def _warp_sun_coordinates(xy, smap, dt, **diffrot_kwargs):
    """
    Function that returns a new list of coordinates for each input coord.
    This is an inverse function needed by the scikit-image `transform.warp`
    function.

    Parameters
    ----------
    xy : `numpy.ndarray`
        Array from `transform.warp`
    smap : `~sunpy.map`
        Original map that we want to transform
    dt : `~astropy.units.Quantity`
        Desired interval to rotate the input map by solar differential rotation.

    Returns
    -------
    xy2 : `~numpy.ndarray`
        Array with the inverse transformation
    """
    # NOTE: The time is being subtracted - this is because this function
    # calculates the inverse of the transformation.
    rotated_time = smap.date - datetime.timedelta(seconds=dt.to(u.s).value)

    # Calculate the hpc coords
    x = np.arange(0, smap.dimensions.x.value)
    y = np.arange(0, smap.dimensions.y.value)
    xx, yy = np.meshgrid(x, y)
    # the xy input array would have the following shape
    # xy = np.dstack([xx.T.flat, yy.T.flat])[0]

    # We start by converting the pixel to world
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        hpc_coords = smap.pixel_to_world(xx * u.pix, yy * u.pix)

        # then diff-rotate the hpc coordinates to the desired time
        rotated_coord = solar_rotate_coordinate(hpc_coords, rotated_time, **diffrot_kwargs)

        # To find the values that are behind the sun we need to convert them
        # to HeliographicStonyhurst
        findOccult = rotated_coord.transform_to(HeliographicStonyhurst)

        with np.errstate(invalid='ignore'):
            # and find which ones are outside the [-90, 90] range.
            occult = np.logical_or(np.less(findOccult.lon, -90 * u.deg),
                                   np.greater(findOccult.lon, 90 * u.deg))

        # NaN-ing values that move to the other side of the sun
        rotated_coord.data.lon[occult] = np.nan * u.deg
        rotated_coord.data.lat[occult] = np.nan * u.deg
        rotated_coord.cache.clear()

        # Go back to pixel co-ordinates
        x2, y2 = smap.world_to_pixel(rotated_coord)

    # Re-stack the data to make it correct output form
    xy2 = np.dstack([x2.T.value.flat, y2.T.value.flat])[0]
    # Returned a masked array with the non-finite entries masked.
    xy2 = np.ma.array(xy2, mask=np.isnan(xy2))
    return xy2


@u.quantity_input(dt='time')
def diffrot_map(smap, time=None, dt=None, pad=False, **diffrot_kwargs):
    """
    Function to apply solar differential rotation to a sunpy map.

    Parameters
    ----------
    smap : `~sunpy.map`
        Original map that we want to transform.
    time : sunpy-compatible time
        date/time at which the input co-ordinate will be rotated to.
    dt : `~astropy.units.Quantity` or `datetime`
        Desired interval between the input map and returned map.
    pad : `bool`
        Whether to create a padded map for submaps to don't loose data

    Returns
    -------
    diffrot_map : `~sunpy.map`
        A map with the result of applying solar differential rotation to the
        input map.
    """
    if (time is not None) and (dt is not None):
        raise ValueError('Only a time or an interval is accepted')
    elif not (time or dt):
        raise ValueError('Either a time or an interval (`dt=`) needs to be provided')
    elif time:
        new_time = parse_time(time)
        dt = (new_time - smap.date).total_seconds() * u.s
    else:
        new_time = smap.date + datetime.timedelta(seconds=dt.to(u.s).value)

    # Check for masked maps
    if smap.mask is not None:
        smap_data = np.ma.array(smap.data, mask=smap.mask)
    else:
        smap_data = smap.data

    submap = False
    # Check whether the input is a submap
    if ((2 * smap.rsun_obs > smap.top_right_coord.Tx - smap.bottom_left_coord.Tx) or
        (2 * smap.rsun_obs > smap.top_right_coord.Ty - smap.bottom_left_coord.Ty)):

        submap = True
        if pad:
            # Calculating the largest distance between the corners and their rotation values
            deltax = deltay = 0
            for corner in product(*product([0 * u.pix], smap.dimensions)):
                corner_world = smap.pixel_to_world(*corner)
                corner_world_rotated = solar_rotate_coordinate(corner_world, new_time, **diffrot_kwargs)
                corner_px_rotated = smap.world_to_pixel(corner_world_rotated)
                dx = np.abs(corner_px_rotated.x - corner[0])
                dy = np.abs(corner_px_rotated.y - corner[1])
                deltax = dx if dx > deltax else deltax
                deltay = dy if dy > deltay else deltay

            deltax = np.int(np.ceil(deltax.value))
            deltay = np.int(np.ceil(deltay.value))
            # Create a new `smap` with the padding around it
            smap_data = np.pad(smap.data, ((deltay, deltay), (deltax, deltax)),
                               'constant', constant_values=0)
            smap_meta = deepcopy(smap.meta)
            smap_meta['naxis2'], smap_meta['naxis1'] = smap_data.shape
            smap_meta['crpix1'] += deltax
            smap_meta['crpix2'] += deltay
            smap = sunpy.map.Map(smap_data, smap_meta)

    warp_args = {'smap': smap, 'dt': dt}
    warp_args.update(diffrot_kwargs)
    # Apply solar differential rotation as a scikit-image warp
    out = transform.warp(to_norm(smap_data), inverse_map=_warp_sun_coordinates,
                         map_args=warp_args)

    # Recover the original intensity range.
    out = un_norm(out, smap.data)

    # Update the meta information with the new date and time, and reference pixel.
    out_meta = deepcopy(smap.meta)
    if out_meta.get('date_obs', False):
        del out_meta['date_obs']
    out_meta['date-obs'] = "{:%Y-%m-%dT%H:%M:%S}".format(new_time)

    if submap:
        crval_rotated = solar_rotate_coordinate(smap.reference_coordinate, new_time, **diffrot_kwargs)
        out_meta['crval1'] = crval_rotated.Tx.value
        out_meta['crval2'] = crval_rotated.Ty.value

    return sunpy.map.Map((out, out_meta))
