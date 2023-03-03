"""
This submodule provides utility functions to act on `sunpy.map.GenericMap` instances.
"""
import numbers
from itertools import product

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.coordinates import Helioprojective, sun

__all__ = ['all_pixel_indices_from_map', 'all_coordinates_from_map',
           'all_corner_coords_from_map',
           'map_edges', 'solar_angular_radius', 'sample_at_coords',
           'contains_full_disk', 'is_all_off_disk', 'is_all_on_disk',
           'contains_limb', 'coordinate_is_on_solar_disk',
           'on_disk_bounding_coordinates',
           'contains_coordinate', 'contains_solar_center',
           'extract_along_coord']


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
    y, x = np.indices(smap.data.shape)
    return [x, y] * u.pix


def all_coordinates_from_map(smap):
    """
    Returns the coordinates of the center of every pixel in a map.

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


def all_corner_coords_from_map(smap):
    """
    Returns the coordinates of the pixel corners in a map.
    """
    ny, nx = smap.data.shape
    y, x = np.indices((ny + 1, nx + 1))
    return smap.pixel_to_world((x - 0.5) * u.pix, (y - 0.5) * u.pix)


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


def _verify_coordinate_helioprojective(coordinates):
    """
    Raises an error if the coordinate is not in the
    `~sunpy.coordinates.frames.Helioprojective` frame.

    Parameters
    ----------
    coordinates : `~astropy.coordinates.SkyCoord`, `~astropy.coordinates.BaseCoordinateFrame`
    """
    frame = coordinates.frame if hasattr(coordinates, 'frame') else coordinates
    if not isinstance(frame, Helioprojective):
        raise ValueError(f"The input coordinate(s) is of type {type(frame).__name__}, "
                         "but must be in the Helioprojective frame.")


def solar_angular_radius(coordinates):
    """
    Calculates the solar angular radius as seen by the observer.

    The tangent vector from the observer to the edge of the Sun forms a
    right-angle triangle with the radius of the Sun as the far side and the
    Sun-observer distance as the hypotenuse.  Thus, the sine of the angular
    radius of the Sun is ratio of these two distances.

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
    _verify_coordinate_helioprojective(coordinates)
    return sun._angular_radius(coordinates.rsun, coordinates.observer.radius)


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


def _edge_coordinates(smap):
    # Calculate all the edge pixels
    edges = map_edges(smap)
    # We need to strip the units from edges before handing it to np.concatenate,
    # as the .unit attribute is not propagated in np.concatenate for numpy<1.17
    # When sunpy depends on numpy>=1.17 this unit replacing code can be removed
    edge_pixels = u.Quantity(np.concatenate(edges).value, unit=u.pix, copy=False)
    # Calculate the edge of the world
    return smap.pixel_to_world(edge_pixels[:, 0], edge_pixels[:, 1])


def contains_full_disk(smap):
    """
    Checks if a map contains the full disk of the Sun.

    A map contains the full disk of the Sun if the following two
    conditions are met: (1) all the coordinates at the edge of the map are
    more than solar angular radius from the center of the Sun and, (2) the
    map is not all off disk. If both these conditions are met, the
    function returns `True`. Otherwise, the function returns `False`.

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
    _verify_coordinate_helioprojective(smap.coordinate_frame)
    edge_of_world = _edge_coordinates(smap)
    # Calculate the distance of the edge of the world in solar radii
    coordinate_angles = np.sqrt(edge_of_world.Tx ** 2 + edge_of_world.Ty ** 2)

    # Test if all the edge pixels are more than one solar radius distant
    # and that the whole map is not all off disk.
    return np.all(coordinate_angles > solar_angular_radius(edge_of_world)) and contains_solar_center(smap)


def contains_solar_center(smap):
    """
    Returns `True` if smap contains the solar center.

    This is the case if and only if the solar center is inside or on the edges of the map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A map in helioprojective Cartesian coordinates.

    Returns
    -------
    bool
        True if the map contains the solar center.
    """
    _verify_coordinate_helioprojective(smap.coordinate_frame)
    return contains_coordinate(smap, SkyCoord(0*u.arcsec, 0*u.arcsec, frame=smap.coordinate_frame))


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
    _verify_coordinate_helioprojective(coordinates)
    # Calculate the angle of every pixel from the center of the Sun and compare it the angular
    # radius of the Sun.
    return np.sqrt(coordinates.Tx ** 2 + coordinates.Ty ** 2) < solar_angular_radius(coordinates)


def is_all_off_disk(smap):
    """
    Checks if none of the coordinates in the `~sunpy.map.GenericMap` are on the solar disk.

    This is done by checking if the edges of the map do not contain the solar limb, and
    checking that the solar center is not in the map.

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
    _verify_coordinate_helioprojective(smap.coordinate_frame)
    edge_of_world = _edge_coordinates(smap)
    # Calculate the distance of the edge of the world in solar radii
    coordinate_angles = np.sqrt(edge_of_world.Tx ** 2 + edge_of_world.Ty ** 2)

    # Test if all the edge pixels are more than one solar radius distant
    # and that the solar center is
    return np.all(coordinate_angles > solar_angular_radius(edge_of_world)) and ~contains_solar_center(smap)


def is_all_on_disk(smap):
    """
    Checks if all of the coordinates in the `~sunpy.map.GenericMap` are on the solar disk.

    The check is performed by calculating the angle of the edges of the map from
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
    _verify_coordinate_helioprojective(smap.coordinate_frame)
    edge_of_world = _edge_coordinates(smap)
    return np.all(coordinate_is_on_solar_disk(edge_of_world))


def contains_limb(smap):
    """
    Checks if a map contains any part of the solar limb or equivalently whether
    the map contains both on-disk and off-disk pixels.

    The check is performed by calculating the angular distance of the edge pixels from
    the center of the Sun. If at least one edge pixel is on disk (less than the solar
    angular radius) and at least one edge pixel is off disk (greater than the solar
    angular distance), or the map contains the full disk, the function returns `True`.
    Otherwise, the function returns `False`.

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
    _verify_coordinate_helioprojective(smap.coordinate_frame)
    if contains_full_disk(smap):
        return True
    on_disk = coordinate_is_on_solar_disk(_edge_coordinates(smap))
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
    _verify_coordinate_helioprojective(smap.coordinate_frame)
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


def contains_coordinate(smap, coordinates):
    """
    Checks whether a coordinate falls within the bounds of a map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        The input map.
    coordinates : `~astropy.coordinates.SkyCoord`
        The input coordinate.

    Returns
    -------
    bool
        `True` if ``coordinates`` falls within the bounds of ``smap``.
        This includes the edges of the map. If multiple coordinates are input,
        returns a boolean array.
    """
    # Dimensions of smap
    ys, xs = smap.wcs.array_shape
    # Converting coordinates to pixels
    xc, yc = smap.wcs.world_to_pixel(coordinates)
    return ((xc >= -0.5) &
            (xc <= xs - 0.5) &
            (yc >= -0.5) &
            (yc <= ys - 0.5))


def _bresenham(*, x1, y1, x2, y2):
    """
    Returns an array of all pixel coordinates which the line defined by `x1, y1` and
    `x2, y2` crosses. Uses Bresenham's line algorithm to enumerate the pixels along
    a line. This was adapted from ginga.

    Parameters
    ----------
    x1, y1, x2, y2 :`int`

    References
    ----------
    * https://github.com/ejeschke/ginga/blob/c8ceaf8e559acc547bf25661842a53ed44a7b36f/ginga/BaseImage.py#L503
    * http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
    """
    for x in [x1, y1, x2, y2]:
        if not isinstance(x, numbers.Integral):
            raise TypeError('All pixel coordinates must be of type int')
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    err = dx - dy
    res = []
    x, y = x1, y1
    while True:
        res.append((x, y))
        if (x == x2) and (y == y2):
            break
        e2 = 2 * err
        if e2 > -dy:
            err = err - dy
            x += sx
        if e2 < dx:
            err = err + dx
            y += sy
    return np.array(res)


def extract_along_coord(smap, coord):
    """
    Extract pixel values from a map along a path that approximates a coordinate path.

    Each pair of consecutive coordinates in the provided coordinate array defines a
    line segment in pixel space.  The approximate pixel path for each line segment is
    determined by applying
    `Bresenham's line algorithm <http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm>`_.
    The resulting pixel path does not necessarily include every pixel that is
    intersected by each line segment.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        The sunpy map.
    coord : `~astropy.coordinates.SkyCoord`
        The coordinate array that defines the path for extraction.

    Returns
    -------
    pixel_values : `~astropy.units.Quantity`
         The pixel values along a path that approximates the coordinate path.
    pixel_coords : `~astropy.coordinates.SkyCoord`
         The coordinates for the pixels from which the values were extracted.

    Notes
    -----
    The provided coordinates are first rounded to the nearest corresponding pixel,
    which means that the coordinates used for calculations may be shifted relative
    to the provided coordinates by up to half a pixel.

    Examples
    --------
    .. minigallery:: sunpy.map.extract_along_coord
    """
    if not len(coord.shape) or coord.shape[0] < 2:
        raise ValueError('At least two points are required for extracting intensity along a '
                         'line. To extract points at single coordinates, use '
                         'sunpy.map.maputils.sample_at_coords.')
    if not all(contains_coordinate(smap, coord)):
        raise ValueError('At least one coordinate is not within the bounds of the map.'
                         'To extract the intensity along a coordinate, all points must fall within '
                         'the bounds of the map.')
    # Find pixels between each loop segment
    px, py = smap.wcs.world_to_array_index(coord)
    pix = []
    for i in range(len(px)-1):
        b = _bresenham(x1=px[i], y1=py[i], x2=px[i+1], y2=py[i+1])
        # Pop the last one, unless this is the final entry because the first point
        # of the next section will be the same
        if i < (len(px) - 2):
            b = b[:-1]
        pix.append(b)
    pix = np.vstack(pix)

    pixel_values = u.Quantity(smap.data[pix[:, 0], pix[:, 1]], smap.unit)
    pixel_coords = smap.wcs.pixel_to_world(pix[:, 1], pix[:, 0])

    return pixel_values, pixel_coords
