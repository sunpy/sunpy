#
# Calculates the co-ordinates along great arcs between two specified points
# which are assumed to be on disk.
#
from __future__ import absolute_import, division, print_function

import numpy as np
from astropy.coordinates import SkyCoord

__all__ = ['great_arc']


def great_arc(start_point, end_point, center=None, number_points=100):
    """
    Calculates a user-specified number of points on an arc of a great circle
    connecting two points.

    Parameters
    ----------
    start_point : `~astropy.coordinates.SkyCoord`
        Start point.

    end_point : `~astropy.coordinates.SkyCoord`
        End point.

    center : `~astropy.coordinates.SkyCoord`
        Center of the sphere.

    number_points : int
        Number of points along the great arc.

    Returns
    -------
    arc : `~astropy.coordinates.SkyCoord`
        Co-ordinates along the great arc in the co-ordinate frame of the
        start point.

    Example
    -------
    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import great_arc
    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE
    >>> m = sunpy.map.Map(AIA_171_IMAGE)
    >>> a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=m.coordinate_frame)
    >>> b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=m.coordinate_frame)
    >>> v = great_arc(a, b)
    """
    start_point_unit = start_point.transform_to('heliocentric').cartesian.xyz.unit
    if center is None:
        c = SkyCoord(0*start_point_unit, 0*start_point_unit, 0*start_point_unit, frame='heliocentric')
    a_cartesian = start_point.transform_to('heliocentric').cartesian.xyz.to(start_point_unit).value
    b_cartesian = end_point.transform_to('heliocentric').cartesian.xyz.to(start_point_unit).value
    c_cartesian = c.transform_to('heliocentric').cartesian.xyz.to(start_point_unit).value

    # Calculate the points along the great arc.
    v_cartesian = _calculate_great_arc(a_cartesian, b_cartesian, c_cartesian, number_points)*start_point_unit

    # Transform the great arc back into the input frame.
    return SkyCoord(v_cartesian[:, 0], v_cartesian[:, 1], v_cartesian[:, 2],
                    frame='heliocentric', observer=start_point.observer).transform_to(start_point.frame)


def _calculate_great_arc(start, end, center, number_points):
    """
    Calculate a user-specified number of points on a great arc between a start
    and end point on a sphere where the start and end points are assumed to be
    x,y,z Cartesian triples on a sphere relative to a center.  See the
    references below for a description of the algorithm

    Parameters
    ----------
    start : `~numpy.ndarray`
        Start point expressed as a Cartesian xyz triple.

    end : `~numpy.ndarray`
        End point expressed as a Cartesian xyz triple.

    center : `~numpy.ndarray`
        Center of the sphere expressed as a Cartesian xyz triple

    number_points : int
        Number of points along the great arc.

    Returns
    -------
    arc : `~numpy.ndarray`
        Co-ordinates along the great arc expressed as Cartesian xyz triples.
        The shape of the array is (number_points, 3).

    References
    ----------
    [1] https://www.mathworks.com/matlabcentral/newsreader/view_thread/277881
    [2] https://en.wikipedia.org/wiki/Great-circle_distance#Vector_version

    """
    # Vector from center to first point
    v1 = start - center

    # Distance of the first point from the center
    r = np.linalg.norm(v1)

    # Vector from center to second point
    v2 = end - center

    # The v3 lies in plane of v1 & v2 and is orthogonal to v1
    v3 = np.cross(np.cross(v1, v2), v1)

    # Ensure that the vector has length r
    v3 = r * v3 / np.linalg.norm(v3)

    # Range through the inner angle between v1 and v2
    inner_angles = np.linspace(0, np.arctan2(np.linalg.norm(np.cross(v1, v2)),
                                             np.dot(v1, v2)), num=number_points).reshape(number_points, 1)

    # Calculate the Cartesian locations from the first to second points
    return v1[np.newaxis, :] * np.cos(inner_angles) + \
           v3[np.newaxis, :] * np.sin(inner_angles) + center
