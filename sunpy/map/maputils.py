"""
This submodule provides utility functions to act on `sunpy.map.GenericMap` instances.
"""
from itertools import chain, product

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.map import GenericMap
from sunpy.coordinates import Helioprojective

__all__ = ['all_pixel_indices_from_map', 'all_coordinates_from_map',
           'map_edges', 'contains_full_disk',
           'is_all_off_disk', 'is_all_on_disk', 'contains_limb',
           'coordinate_is_on_disk', 'on_disk_bounding_coordinates']


def all_pixel_indices_from_map(smap):
    """
    Returns pixel pair indices of every pixel in a map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
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
    smap : `~sunpy.map.GenericMap`
        A SunPy map.

    Returns
    -------
    out : `~astropy.coordinates.SkyCoord`
        An array of sky coordinates in the coordinate system "coordinate_system".
    """
    x, y = all_pixel_indices_from_map(smap)
    return smap.pixel_to_world(x, y)


def map_edges(smap):
    """
    Returns the pixel locations of the edges of a rectangular input map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        The input map.

    Returns
    -------
    maps_edges : `~dict`
        Returns the pixels of edge of the map
    """
    # Calculate all the edge pixels
    nx, ny = smap.dimensions.x.value, smap.dimensions.y.value
    bottom = list(product([0.0], np.arange(nx))) * u.pix
    top = list(product([ny - 1], np.arange(nx))) * u.pix
    lhs = list(product(np.arange(ny), [0])) * u.pix
    rhs = list(product(np.arange(ny), [nx - 1])) * u.pix
    return {"top": top, "bottom": bottom, "lhs": lhs, "rhs": rhs}


def contains_full_disk(smap):
    """
    Checks if a map contains the full disk of the Sun.

    A map contains the full disk of the Sun if the following two
    conditions are met: (1) all the pixels at the edge of the map are
    more than 1 solar radius from the center of the Sun and, (2) the
    map is not all off disk.  If both these conditions are met, the
    function returns True. Otherwise, the function returns False.
    Note that the function assumes that the input map is rectangular.
    Note also that in the case of coronagraph images the disk itself
    need not be observed.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        The input map.

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
    edge_of_world = smap.pixel_to_world(x, y)

    # Calculate the distance of the edge of the world in solar radii
    distance = u.R_sun * np.sqrt(edge_of_world.Tx ** 2 + edge_of_world.Ty ** 2) / smap.rsun_obs

    # Test if all the edge pixels are more than one solar radius distant
    # and that the whole map is not all off disk.
    return np.all(distance > 1*u.R_sun) and ~is_all_off_disk(smap)


def is_all_off_disk(smap):
    """
    Checks if the entire map is off the solar disk.

    The check is performed by calculating the distance of every
    pixel from the center of the Sun.  If they are all off-disk,
    then the function returns True. Otherwise, the function
    returns False.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        The input map.

    Returns
    -------
    is_all_off_disk : `~bool`
        Returns True if all map pixels are strictly more than one solar radius
        away from the center of the Sun.
    """
    return np.all(~coordinate_is_on_disk(all_coordinates_from_map(smap), smap.rsun_obs))


def is_all_on_disk(smap):
    """
    Checks if the entire map is on the solar disk.

    The check is performed by calculating the distance of every pixel
    from the center of the Sun.  If they are all on-disk, then the
    function returns True. Otherwise, the function returns False.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        The input map.

    Returns
    -------
    is_all_on_disk : `~bool`
        Returns True if all map pixels are strictly less than one solar radius
        away from the center of the Sun.
    """
    return np.all(coordinate_is_on_disk(all_coordinates_from_map(smap), smap.rsun_obs))


def contains_limb(smap):
    """
    Checks if a map contains a portion of the solar limb.

    The check is performed by calculating the distance of every pixel from
    the center of the Sun.  If at least one pixel is on disk and at least one
    pixel is off disk, the function returns True.  Otherwise, the function
    returns False. Note that this function will also true if the entire solar
    limb is within the field of view of the map.  Note also that in the case
    of coronagraph images the limb itself need not be observed.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        The input map.

    Returns
    -------
    contains_limb : `~bool`
        Returns True If at least one pixel is on disk and at least one pixel
        is off disk.
    """
    on_disk = coordinate_is_on_disk(all_coordinates_from_map(smap), smap.rsun_obs)
    return np.logical_and(np.any(on_disk), np.any(~on_disk))


@u.quantity_input
def coordinate_is_on_disk(coordinate, scale: u.arcsecond):
    """
    Checks if the helioprojective Cartesian coordinate is on disk.

    The check is performed by comparing the coordinate's distance
    from the center of the Sun to the solar radius.

    Parameters
    ----------
    coordinate : `~astropy.coordinates.SkyCoord`, `~sunpy.coordinates.frames.Helioprojective`
        The input coordinate.  The frame property of the input
        coordinate must be helioprojective.

    scale : `~astropy.units.Quantity`
        The pixel scale size in units of arcseconds.

    Returns
    -------
    is_on_disk : `~bool`
        Returns True if the coordinate is on disk, False otherwise.
    """
    # Calculate the radii of every pixel in helioprojective Cartesian
    # co-ordinate distance units.
    return u.R_sun * (np.sqrt(coordinate.Tx ** 2 + coordinate.Ty ** 2) / scale) < 1 * u.R_sun


def on_disk_bounding_coordinates(smap):
    """
    Returns the the bottom left and top-right coordinates of the smallest
    rectangular region that contains all the on-disk pixels in the input map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        The input map.

    Returns
    -------
    bl, tr : `~astropy.coordinates.SkyCoord`
        Returns a SkyCoord of length 2 such that the first entry is the bottom left
        coordinate and the second entry is the top right coordinate of the smallest
        rectangular region that contains all the on-disk pixels in the input map.
    """
    # Check that the input map is not all off disk.
    if is_all_off_disk(smap):
        raise ValueError("The entire map is off disk.")

    # Get all the coordinates from the input map
    coordinates = all_coordinates_from_map(smap)

    # Find which coordinates are on the disk
    on_disk = coordinate_is_on_disk(coordinates, smap.rsun_obs)
    on_disk_coordinates = coordinates[on_disk]

    # The bottom left and top right coordinates that contain
    # the on disk coordinates.
    tx = on_disk_coordinates.Tx.value
    ty = on_disk_coordinates.Ty.value
    return SkyCoord([np.nanmin(tx), np.nanmax(tx)] * u.arcsec,
                    [np.nanmin(ty), np.nanmax(ty)] * u.arcsec,
                    frame=Helioprojective, observer=smap.observer_coordinate)
