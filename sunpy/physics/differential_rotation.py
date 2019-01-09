from itertools import chain, product
import warnings
from copy import deepcopy

import numpy as np

from astropy import units as u
from astropy.time import TimeDelta
from astropy.coordinates import SkyCoord, Longitude, BaseCoordinateFrame, get_body

from sunpy.time import parse_time, is_time
from sunpy.coordinates import frames


__all__ = ['diff_rot', 'solar_rotate_coordinate', 'diffrot_map',
           'all_pixel_indices_from_map', 'all_coordinates_from_map',
           'find_pixel_radii', 'map_edges', 'contains_full_disk',
           'is_all_off_disk']


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

    # This calculation of the rotation assumes a sidereal frame time.
    rotation = (A + B * sin2l + C * sin4l) * duration

    # Applying this correction assumes that the observer is on the Earth,
    # and that the Earth is at the same distance from the Sun at all times
    # during the year.
    if frame_time == 'synodic':
        rotation -= 0.9856 * u.deg / u.day * duration

    return Longitude(rotation.to(u.deg))


def solar_rotate_coordinate(coordinate, observer=None, time=None, **diff_rot_kwargs):
    """
    Given a coordinate on the Sun, calculate where that coordinate maps to
    at as seen by a new observer at some later or earlier time, given that
    the input coordinate rotates according to the solar rotation profile.

    The amount of solar rotation is based on the amount of time between the
    observation time of the input coordinate and the observation time of the
    new observer.

    Parameters
    ----------
    coordinate : `~astropy.coordinates.SkyCoord`
        Any valid coordinate which is transformable to Heliographic Stonyhurst.

    observer : `~astropy.coordinates.BaseCoordinateFrame`, `~astropy.coordinates.SkyCoord`
        The location of the new observer.
        Instruments in Earth orbit can be approximated by using the position
        of the Earth at the observation time of the new observer.

    time : `~astropy.time.Time`


    **diff_rot_kwargs : keyword arguments
        Keyword arguments are passed on as keyword arguments to `~sunpy.physics.differential_rotation.diff_rot`.
        Note that the keyword "frame_time" is automatically set to the value
        "sidereal".

    Returns
    -------
    coordinate : `~astropy.coordinates.SkyCoord`
        The locations of the input coordinates after the application of
        solar rotation as seen from the point-of-view of the new observer.

    Example
    -------
    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord, get_body
    >>> from sunpy.coordinates import frames
    >>> from sunpy.physics.differential_rotation import solar_rotate_coordinate
    >>> from sunpy.time import parse_time
    >>> start_time = parse_time('2010-09-10 12:34:56')
    >>> end_time = parse_time('2010-09-10 13:34:56')
    >>> c = SkyCoord(-570*u.arcsec, 120*u.arcsec, obstime=start_time, frame=frames.Helioprojective)
    >>> solar_rotate_coordinate(c, time=end_time)
    <SkyCoord (Helioprojective: obstime=2010-09-10T13:34:56.000, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
    (-562.89877819, 119.3152842, 1.50085078e+08)>
    >>> new_observer = get_body("earth", end_time)
    >>> solar_rotate_coordinate(c, observer=new_observer)
    <SkyCoord (Helioprojective: obstime=2010-09-10T13:34:56.000, rsun=695508.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, km)
        (-562.89877819, 119.3152842, 1.50085078e+08)>
    >>> new_observer = get_body("mars", end_time)
    >>> solar_rotate_coordinate(c, observer=new_observer)
    ???
    """
    # Check the input and create the new observer
    if (observer is not None) and (time is not None):
        raise ValueError("Either the 'observer' or the 'time' keyword must be specified, but not both simultaneously.")

    if observer is not None:
        # Check that the new_observer is specified correctly.
        if not (isinstance(observer, (BaseCoordinateFrame, SkyCoord))):
            raise ValueError(
                "The 'observer' must be an astropy.coordinates.BaseCoordinateFrame or an astropy.coordinates.SkyCoord.")

        new_observer = observer

    if time is not None:
        warnings.warn("Using 'time' assumes an Earth-based observer.")
        if isinstance(time, TimeDelta) or isinstance(time, u.Quantity):
            new_observer_time = coordinate.obstime + time
        else:
            new_observer_time = parse_time(time)

        new_observer = get_body("earth", new_observer_time)

    # The keyword "frame_time" must be explicitly set to "sidereal"
    # when using this function.
    diff_rot_kwargs.update({"frame_time": "sidereal"})

    # Calculate the interval between the start and end time
    interval = (new_observer.obstime - coordinate.obstime).to(u.s)

    # Compute Stonyhurst Heliographic co-ordinates - returns (longitude,
    # latitude). Points off the limb are returned as nan.
    heliographic_coordinate = coordinate.transform_to(frames.HeliographicStonyhurst)

    # Compute the differential rotation
    drot = diff_rot(interval, heliographic_coordinate.lat.to(u.degree), **diff_rot_kwargs)

    # Rotate the input co-ordinate as seen by the original observer
    heliographic_rotated = SkyCoord(heliographic_coordinate.lon + drot, heliographic_coordinate.lat,
                                    obstime=coordinate.obstime, observer=coordinate.observer,
                                    frame=frames.HeliographicStonyhurst)

    # Calculate where the rotated co-ordinate appears as seen by new observer,
    # and then transform it into the co-ordinate system of the input
    # co-ordinate.
    return heliographic_rotated.transform_to(new_observer).transform_to(coordinate.frame.name)


def _warp_sun_coordinates(xy, smap, observer, **diffrot_kwargs):
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

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Get all the co-ordinates from the map
        hpc_coords = all_coordinates_from_map(smap)
        print(xy)

        # then diff-rotate the hpc coordinates to the desired time
        rotated_coord = solar_rotate_coordinate(hpc_coords, observer=observer, **diffrot_kwargs)

        # To find the values that are behind the sun we need to convert them
        # to HeliographicStonyhurst
        find_occult = rotated_coord.transform_to(frames.HeliographicStonyhurst)

        with np.errstate(invalid='ignore'):
            # and find which ones are outside the [-90, 90] range.
            occult = np.logical_or(np.less(find_occult.lon, -90 * u.deg),
                                   np.greater(find_occult.lon, 90 * u.deg))

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


def diffrot_map(smap, observer=None, time=None, **diffrot_kwargs):
    """
    Function to apply solar differential rotation to a sunpy map.

    Parameters
    ----------
    smap : `~sunpy.map`
        Original map that we want to transform.

    observer : `~astropy.coordinates.BaseCoordinateFrame`, `~astropy.coordinates.SkyCoord`
        The location of the new observer.
        Instruments in Earth orbit can be approximated by using the position
        of the Earth at the observation time of the new observer.

    time : sunpy-compatible time
        date/time at which the input co-ordinate will be rotated to.

    Returns
    -------
    diffrot_map : `~sunpy.map`
        A map with the result of applying solar differential rotation to the
        input map.
    """

    # If the entire map is off-disk, then there is nothing to do.
    if is_all_off_disk(smap):
        raise ValueError("The entire map is off disk. No data to differentially rotate.")

    # Only this function needs scikit image
    from skimage import transform
    from sunpy.image.util import to_norm, un_norm
    # Import map here for performance reasons.
    import sunpy.map

    # Check for masked maps
    if smap.mask is not None:
        smap_data = np.ma.array(smap.data, mask=smap.mask)
    else:
        smap_data = smap.data

    # At least part of the input map is on the disk.
    # Check whether the input contains the full disk of the Sun
    submap = not contains_full_disk(smap)
    if submap:
        # Get the edges of the map
        edges = map_edges(smap)

        # Calculate the size of the output array.
        # Calculate the difference between the top and the bottom.
        # Rotate the top and bottom edges
        rotated_top = solar_rotate_coordinate(smap.pixel_to_world(*edges["top"]), observer=observer, time=time, **diffrot_kwargs)
        rotated_bottom = solar_rotate_coordinate(smap.pixel_to_world(*edges["bottom"]), observer=observer, time=time, **diffrot_kwargs)

        # Calculate the difference between the rotated top and bottom
        difference_top_bottom_x = np.abs(rotated_top.Tx - rotated_bottom.Tx)
        difference_top_bottom_y = np.abs(rotated_top.Ty - rotated_bottom.Ty)

        # Calculate the difference between the left and right hand side.
        # Rotate the left and right hand edges
        rotated_lhs = solar_rotate_coordinate(smap.pixel_to_world(*edges["lhs"]), observer=observer, time=time, **diffrot_kwargs)
        rotated_rhs = solar_rotate_coordinate(smap.pixel_to_world(*edges["rhs"]), observer=observer, time=time, **diffrot_kwargs)

        # Calculate the difference between the rotated left and right hand sides.
        difference_lhs_rhs_x = np.abs(rotated_lhs.Tx - rotated_rhs.Tx)
        difference_lhs_rhs_y = np.abs(rotated_lhs.Ty - rotated_rhs.Ty)

        # Calculate the size of the bounding box
        rotated_nx = int(np.ceil(np.max([difference_top_bottom_x, difference_lhs_rhs_x]) / smap.scale.axis1).value)
        rotated_ny = int(np.ceil(np.max([difference_top_bottom_y, difference_lhs_rhs_y]) / smap.scale.axis2).value)

        # Change in size of the padded array relative to the input map
        deltax = np.abs(rotated_nx - smap.data.shape[0])
        deltay = np.abs(rotated_ny - smap.data.shape[1])

        # Create a new `smap` with the padding around it
        smap_data = np.pad(smap.data, ((deltay, deltay), (deltax, deltax)), 'constant', constant_values=0)
        smap_meta = deepcopy(smap.meta)
        smap_meta['naxis2'], smap_meta['naxis1'] = smap_data.shape

        smap_meta['crpix1'] += deltax
        smap_meta['crpix2'] += deltay

        # Create the padded map that will be used to create the rotated map.
        smap = sunpy.map.Map(smap_data, smap_meta)

    warp_args = {'smap': smap, 'new_observer': new_observer}
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
    out_meta['date-obs'] = "{:%Y-%m-%dT%H:%M:%S}".format(new_observer.obstime)

    if submap:
        # Put the reference pixel at (0, 0)
        out_meta['crpix1'] = 0  # a proper calculation of crpix1 and crpix2 is required
        out_meta['crpix2'] = 0

        # Calculate where the center of the field of view is
        crval_rotated = solar_rotate_coordinate(smap.pixel_to_world(0 * u.pix, 0 * u.pix), observer=observer, time=time, **diffrot_kwargs)

        # Calculate where the center of the field of view is
        crval_rotated = solar_rotate_coordinate(smap.reference_coordinate, observer=observer, time=time, **diffrot_kwargs)
        out_meta['crval1'] = crval_rotated.Tx.value
        out_meta['crval2'] = crval_rotated.Ty.value

    return sunpy.map.Map((out, out_meta))


# Functions that calculate useful quantities from maps. The functions
# all_pixel_indices_from_map, all_coordinates_from_map and find_pixel_radii
# were originally written for sunkit-image
def all_pixel_indices_from_map(smap):
    """
    Returns pixel pair indices of every pixel in a map.

    Parameters
    ----------
    smap : `~sunpy.map.Map`
        A SunPy map.

    Returns
    -------
    out : `~numpy.array`
        A numpy array with the all the pixel indices built from the
        dimensions of the map.
    """
    return np.meshgrid(*[np.arange(v.value) for v in smap.dimensions]) * u.pix


def all_coordinates_from_map(smap):
    """
    Returns the co-ordinates of every pixel in a map.

    Parameters
    ----------
    smap : `~sunpy.map.Map`
        A SunPy map.

    Returns
    -------
    out : `~astropy.coordinates.SkyCoord`
        An array of sky coordinates in the coordinate system "coordinate_system".
    """
    x, y = all_pixel_indices_from_map(smap)
    return smap.pixel_to_world(x, y)


def find_pixel_radii(smap, scale=None):
    """
    Find the distance of every pixel in a map from the center of the Sun.
    The answer is returned in units of solar radii.

    Parameters
    ----------
    smap : `~sunpy.map.Map`
        A SunPy map.

    scale : None | `~astropy.units.Quantity`
        The radius of the Sun expressed in map units.  For example, in typical
        helioprojective Cartesian maps the solar radius is expressed in units
        of arcseconds.  If None then the map is queried for the scale.

    Returns
    -------
    radii : `~astropy.units.Quantity`
        An array the same shape as the input map.  Each entry in the array
        gives the distance in solar radii of the pixel in the corresponding
        entry in the input map data.
    """

    # Calculate the helioprojective Cartesian co-ordinates of every pixel.
    coords = all_coordinates_from_map(smap).transform_to(frames.Helioprojective)

    # Calculate the radii of every pixel in helioprojective Cartesian
    # co-ordinate distance units.
    radii = np.sqrt(coords.Tx ** 2 + coords.Ty ** 2)

    # Re-scale the output to solar radii
    if scale is None:
        return u.R_sun * (radii / smap.rsun_obs)
    else:
        return u.R_sun * (radii / scale)


def map_edges(smap):
    """
    Returns the pixel locations of the edges of a rectangular input map.

    Parameters
    ----------
    smap : `~sunpy.map`
        The input map

    Returns
    -------
    maps_edges : `~dict`
        Returns the pixels of edge of the map
    """
    # Calculate all the edge pixels
    nx, ny = smap.dimensions.x.value, smap.dimensions.y.value
    top = list(product([0.0], np.arange(nx))) * u.pix
    bottom = list(product([ny - 1], np.arange(nx))) * u.pix
    lhs = list(product(np.arange(ny), [0])) * u.pix
    rhs = list(product(np.arange(ny), [nx - 1])) * u.pix
    return {"top": top, "bottom": bottom, "lhs": lhs, "rhs": rhs}


def contains_full_disk(smap):
    """
    Checks if a map contains the full disk of the Sun.  A map contains
    the full disk of the Sun if the following two conditions are met:
    (1) all the pixels at the edge of the map are more than 1 solar
    radius from the center of the Sun and, (2) the map is not all off
    disk.  If both these conditions are met, the function returns True.
    Otherwise, the function returns False.  Note that the function
    assumes that the input map is rectangular.  Note also that in the
    case of coronagraph images the disk itself need not be observed.

    Parameters
    ----------
    smap : `~sunpy.map`
        The input map

    Returns
    -------
    contains_full_disk : `~bool`
        Returns False if any of the edge pixels are less than one solar radius
        away from the center of the Sun.
    """
    # Calculate all the edge pixels
    edges = map_edges(smap)
    edge_pixels = list(chain.from_iterable([edges["lhs"], edges["rhs"], edges["top"], edges["bottom"]]))
    x = [p[0] for p in edge_pixels] * u.pix
    y = [p[1] for p in edge_pixels] * u.pix

    # Calculate the edge of the world
    edge_of_world = smap.pixel_to_world(x, y).transform_to(frames.Helioprojective)

    # Calculate the distance of the edge of the world in solar radii
    distance = u.R_sun * np.sqrt(edge_of_world.Tx ** 2 + edge_of_world.Ty ** 2) / smap.rsun_obs

    # Test if all the edge pixels are more than one solar radius distant
    # and that the whole map is not all off disk.
    return np.all(distance >= 1*u.R_sun) and ~is_all_off_disk(smap)


def is_all_off_disk(smap):
    """
    Checks to see if the entire map is off the solar disk.  The check is
    performed by calculating the distance of every pixel from the center of
    the Sun.  If they are all off-disk, then the function returns True.
    Otherwise, the function returns False.

    Parameters
    ----------
    smap : `~sunpy.map`
        The input map.

    Returns
    -------
    is_all_off_disk : `~bool`
        Returns True if all map pixels are strictly more than one solar radius
        away from the center of the Sun.
    """
    return np.all(find_pixel_radii(smap) > 1 * u.R_sun)


def is_all_on_disk(smap):
    """
    Checks to see if the entire map is on the solar disk.  The check is
    performed by calculating the distance of every pixel from the center of
    the Sun.  If they are all on-disk, then the function returns True.
    Otherwise, the function returns False.

    Parameters
    ----------
    smap : `~sunpy.map`
        The input map.

    Returns
    -------
    is_all_on_disk : `~bool`
        Returns True if all map pixels are strictly less than one solar radius
        away from the center of the Sun.
    """
    return np.all(find_pixel_radii(smap) < 1 * u.R_sun)


def contains_limb(smap):
    """
    Checks to see if a map contains a portion of the solar limb.  The check is
    performed by calculating the distance of every pixel from the center of
    the Sun.  If at least one pixel is on disk and at least one pixel is off
    disk, the function returns True.  Otherwise, the function returns False.
    Note that this function will also true if the entire solar limb is within
    the field of view of the map.  Note also that in the case of coronagraph
    images the limb itself need not be observed.

    Parameters
    ----------
    smap : `~sunpy.map`
        The input map.

    Returns
    -------
    contains_limb : `~bool`
        Returns True if all map pixels are strictly less than one solar radius
        away from the center of the Sun.
    """
    pixel_radii = find_pixel_radii(smap)
    return np.logical_and(np.any(pixel_radii < 1 * u.R_sun), np.any(pixel_radii > 1 * u.R_sun))
