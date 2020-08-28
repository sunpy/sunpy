"""
This submodule provides utility functions to act on `sunpy.map.GenericMap` instances.
"""
from itertools import product

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.coordinates import Helioprojective

__all__ = ['all_pixel_indices_from_map', 'all_coordinates_from_map',
           'map_edges', 'solar_angular_radius', 'sample_at_coords',
           'contains_full_disk', 'is_all_off_disk', 'is_all_on_disk',
           'contains_limb', 'coordinate_is_on_solar_disk',
           'on_disk_bounding_coordinates', '_convert_pixels']


def all_pixel_indices_from_map(smap):
    """
    Returns pixel pair indices of every pixel in a map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.

    Returns
    -------
    `~numpy.array`
        A `numpy.array` with the all the pixel indices built from the
        dimensions of the map.
    """
    return np.meshgrid(*[np.arange(v.value) for v in smap.dimensions]) * u.pix


def all_coordinates_from_map(smap):
    """
    Returns the coordinates of every pixel in a map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.

    Returns
    -------
    `~astropy.coordinates.SkyCoord`
        An two-dimensional array of sky coordinates in the coordinate
        system "coordinate_system".
    """
    x, y = all_pixel_indices_from_map(smap)
    return smap.pixel_to_world(x, y)


def map_edges(smap):
    """
    Returns the pixel locations of the edges of an input map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.

    Returns
    -------
    top, bottom, left_hand_side, right_hand_side : `~astropy.units.Quantity`
        Returns the pixel locations at the edge of the map;
        the zeroth, first, second and third tuple values
        return the top, bottom, left hand side and right
        hand side pixel locations respectively of the input map.
    """
    # Calculate all the edge pixels
    nx, ny = smap.dimensions.x.value, smap.dimensions.y.value
    top = list(product(np.arange(nx), [ny - 1])) * u.pix
    bottom = list(product(np.arange(nx), [0])) * u.pix
    left_hand_side = list(product([0], np.arange(ny))) * u.pix
    right_hand_side = list(product([nx - 1], np.arange(ny))) * u.pix
    return top, bottom, left_hand_side, right_hand_side


@u.quantity_input
def solar_angular_radius(coordinates):
    """
    Calculates the solar angular radius as seen by the observer.

    The tangent of the angular size of the Sun is equal to the radius
    of the Sun divided by the distance between the observer and the
    center of the Sun.

    Parameters
    ----------
    coordinates : `~astropy.coordinates.SkyCoord`, `~sunpy.coordinates.frames.Helioprojective`
        The input coordinate. The coordinate frame must be
        `~sunpy.coordinates.Helioprojective`.

    Returns
    -------
    angle : `~astropy.units.Quantity`
        The solar angular radius.
    """
    return np.arctan(coordinates.rsun / coordinates.observer.radius)


def sample_at_coords(smap, coordinates):
    """
    Samples the data in a map at given series of coordinates.
    Uses nearest-neighbor interpolation of coordinates in map, as
    it effectively uses array indexing.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A SunPy map.
    coordinates : `~astropy.coordinates.SkyCoord`
        Input coordinates.
    Returns
    -------
    `numpy.array`
        A `numpy.array` corresponding to the data obtained from the map,
        at the input coordinates.
    """
    return smap.data[smap.wcs.world_to_array_index(coordinates)]


def _convert_pixels(edge):
    """
    Helper function to convert a list of edge pixels of the form [(x0, y0), (x1, y1),..., (xn,yn)]
    into the form ([x0, x1,...,xn], [y0, y1,...,yn])

    Parameter
    ---------
    edge : list
        A list of pairs of tuples.

    Returns
    -------
    tuple
        A tuple containing two lists. The first entry is a list containing the first coordinate of
        each pixel, and similarly for the second entry.
    """
    x = [p[0] for p in edge.value] * u.pix
    y = [p[1] for p in edge.value] * u.pix
    return x, y


def contains_full_disk(smap):
    """
    Checks if a map contains the full disk of the Sun.

    A map contains the full disk of the Sun if the following two
    conditions are met: (1) all the coordinates at the edge of the map are
    more than solar angular radius from the center of the Sun and, (2) the
    location of the boundaries extend across the disk. If both these
    conditions are met, the function returns `True`. Otherwise, the
    function returns `False`.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A map in helioprojective Cartesian coordinates.

    Returns
    -------
    `~bool`
        Returns `False` if any of the coordinates at the edge of the map
         are less than one solar radius away from the center of the Sun.

    Notes
    -----
    This function checks if the image coordinates include the solar disk.
    Therefore this function would return `True` for a coronagraph image
    such as from LASCO/C3 or STEREO/SECCHI COR1 since the solar disk is
    within the field of the view of the instrument (although no emission
    from the disk itself is present in the data.)
    """
    # Get the edge pixels
    top_, bottom, left_hand_side, right_hand_side = map_edges(smap)

    # Coordinates of the edges
    horizontal1 = smap.pixel_to_world(_convert_pixels(top_))
    horizontal2 = smap.pixel_to_world(_convert_pixels(bottom))
    vertical1 = smap.pixel_to_world(_convert_pixels(left_hand_side))
    vertical2 = smap.pixel_to_world(_convert_pixels(right_hand_side))

    # The radius of the Sun
    radius = smap.rsun_obs

    # Determine which of the inputs are top and bottom edges of the map
    top = None
    bot = None
    h1 = horizontal1.Ty
    if np.all(h1 > radius):
        top = h1
    elif np.all(h1 < -radius):
        bot = h1

    h2 = horizontal2.Ty
    if np.all(h2 > radius):
        top = h2
    elif np.all(h2 < -radius):
        bot = h2

    # If neither edge lies further away than a radius then the map
    # cannot contain the full disk.
    if top is None or bot is None:
        return False

    # Determine which of the inputs are the left and right hand edges of the map
    lhs = None
    rhs = None
    v1 = vertical1.Tx
    if np.all(v1 > radius):
        rhs = v1
    elif np.all(v1 < -radius):
        lhs = v1

    v2 = vertical2.Tx
    if np.all(v2 > radius):
        rhs = v2
    elif np.all(v2 < -radius):
        lhs = v2

    # If neither edge lies further away than a radius then the map
    # cannot contain the full disk.
    if lhs is None or rhs is None:
        return False

    return np.all(top > radius) and np.all(bot < -radius) and np.all(lhs < -radius) and np.all(rhs > radius)


@u.quantity_input
def coordinate_is_on_solar_disk(coordinates):
    """
    Checks if the helioprojective Cartesian coordinates are on the solar disk.

    The check is performed by comparing the coordinate's angular distance
    to the angular size of the solar radius.  The solar disk is assumed to be
    a circle i.e., solar oblateness and other effects that cause the solar disk to
    be non-circular are not taken in to account.

    Parameters
    ----------
    coordinates : `~astropy.coordinates.SkyCoord`, `~sunpy.coordinates.frames.Helioprojective`
        The input coordinate. The coordinate frame must be
        `~sunpy.coordinates.Helioprojective`.

    Returns
    -------
    `~bool`
        Returns `True` if the coordinate is on disk, `False` otherwise.
    """

    if not isinstance(coordinates.frame, Helioprojective):
        raise ValueError('The input coordinate(s) must be in the Helioprojective Cartesian frame.')
    # Calculate the angle of every pixel from the center of the Sun and compare it the angular
    # radius of the Sun.
    return np.sqrt(coordinates.Tx ** 2 + coordinates.Ty ** 2) < solar_angular_radius(coordinates)


def is_all_off_disk(smap):
    """
    Checks if none of the coordinates in the `~sunpy.map.GenericMap` are on the solar disk.

    The check is performed by calculating the angle of every pixel from
    the center of the Sun. If they are all greater than the angular
    radius of the Sun, then the function returns `True`. Otherwise, the function
    returns `False`.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A map in helioprojective Cartesian coordinates.

    Returns
    -------
    `~bool`
        Returns `True` if all map pixels have an angular radius greater than
        the angular radius of the Sun.

    Notes
    -----
    For coronagraph images such as those from LASCO C2 and C3 the full disk is
    within the field of view of the instrument, but the solar disk itself is not imaged.
    For such images this function will return `False`.
    """
    return np.all(~coordinate_is_on_solar_disk(all_coordinates_from_map(smap)))


def is_all_on_disk(smap):
    """
    Checks if all of the coordinates in the `~sunpy.map.GenericMap` are on the solar disk.

    The check is performed by calculating the angle of every pixel from
    the center of the Sun. If they are all less than the angular
    radius of the Sun, then the function returns `True`. Otherwise, the function
    returns `False`.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A map in helioprojective Cartesian coordinates.

    Returns
    -------
    `~bool`
        Returns `True` if all map coordinates have an angular radius less than
        the angular radius of the Sun.
    """
    return np.all(coordinate_is_on_solar_disk(all_coordinates_from_map(smap)))


def contains_limb(smap):
    """
    Checks if a map contains any part of the solar limb or equivalently whether
    the map contains both on-disk and off-disk pixels.

    The check is performed by calculating the angular distance of every pixel from
    the center of the Sun.  If at least one pixel is on disk (less than the solar
    angular radius) and at least one pixel is off disk (greater than the solar
    angular distance), the function returns `True`. Otherwise, the function
    returns `False`.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A map in helioprojective Cartesian coordinates.

    Returns
    -------
    `~bool`
        Returns `True` If at least one coordinate of the map is on disk and at
        least one coordinate of the map is off disk.

    Notes
    -----
    For coronagraph images such as those from LASCO C2 and C3 the full disk is
    within the field of view of the instrument, but the solar disk itself is not imaged.
    For such images this function will return `True`.
    """
    on_disk = coordinate_is_on_solar_disk(all_coordinates_from_map(smap))
    return np.logical_and(np.any(on_disk), np.any(~on_disk))


def on_disk_bounding_coordinates(smap):
    """
    Returns the the bottom left and top right coordinates of the smallest
    rectangular region that contains all the on disk coordinates of the input map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A map in helioprojective Cartesian coordinates.

    Returns
    -------
    `~astropy.coordinates.SkyCoord`
        A `~astropy.coordinates.SkyCoord` of length 2 such that the
        first entry is the bottom left coordinate and the second entry is the
        top right coordinate of the smallest rectangular region that contains
        all the on-disk pixels in the input map.
    """
    # Check that the input map is not all off disk.
    if is_all_off_disk(smap):
        raise ValueError("The entire map is off disk.")

    # Get all the coordinates from the input map
    coordinates = all_coordinates_from_map(smap)

    # Find which coordinates are on the disk
    on_disk = coordinate_is_on_solar_disk(coordinates)
    on_disk_coordinates = coordinates[on_disk]

    # The bottom left and top right coordinates that contain
    # the on disk coordinates.
    tx = on_disk_coordinates.Tx.value
    ty = on_disk_coordinates.Ty.value
    return SkyCoord([np.nanmin(tx), np.nanmax(tx)] * u.arcsec,
                    [np.nanmin(ty), np.nanmax(ty)] * u.arcsec,
                    frame=smap.coordinate_frame)
