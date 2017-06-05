#
# Calculates the co-ordinates along great arcs between two specified points
# which are assumed to be on disk.
#
from __future__ import absolute_import, division, print_function
import numpy as np
import sunpy.coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u


def great_arc(a, b, center=SkyCoord(0*u.km, 0*u.km, 0*u.km, frame='heliocentric'), num=100):
    """
    Calculate a user-specified number of points on a great arc between a start
    and end point on a sphere.

    Parameters
    ----------
    :param a: `~astropy.coordinates.SkyCoord`
        Start point.

    :param b: `~astropy.coordinates.SkyCoord`
        End point.

    :param center: `~astropy.coordinates.SkyCoord`
        Center of the sphere.

    :param num: int
        Number of points along the great arc.

    :return: `~astropy.coordinates.SkyCoord`
        Co-ordinates along the great arc in the co-ordinate frame of the
        start point.

    Example
    -------
    >>> import sunpy.coordinates
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE
    >>> m = sunpy.map.Map(AIA_171_IMAGE)
    >>> a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=m.coordinate_frame)
    >>> b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=m.coordinate_frame)
    >>> v = great_arc(a, b)
    """
    input_frame = a.frame
    a_unit = a.transform_to('heliocentric').cartesian.xyz.unit
    a_xyz = a.transform_to('heliocentric').cartesian.xyz.value
    b_xyz = b.transform_to('heliocentric').cartesian.xyz.to(a_unit).value
    c_xyz = center.transform_to('heliocentric').cartesian.xyz.to(a_unit).value

    # Calculate the points along the great arc.
    v_xyz = calculate_great_arc(a_xyz, b_xyz, c_xyz, num=num)*a_unit

    # Transform the great arc back into the input frame.
    return SkyCoord(v_xyz[:, 0], v_xyz[:, 1], v_xyz[:, 2],
                    frame='heliocentric', observer=a.observer).transform_to(input_frame)


def calculate_great_arc(a, b, c, num):
    """
    Calculate a user-specified number of points on a great arc between a start
    and end point on a sphere where the start and end points are assumed to be
    x,y,z Cartesian triples on a sphere relative to a center.  See the
    references below for a description of the algorithm

    Parameters
    ----------
    :param a: `~numpy.ndarray`
        Start point expressed as a Cartesian xyz triple.

    :param b: `~numpy.ndarray`
        End point expressed as a Cartesian xyz triple.

    :param c: `~numpy.ndarray`
        Center of the sphere expressed as a Cartesian xyz triple

    :param num: int
        Number of points along the great arc.

    :return: `~numpy.ndarray`
        Co-ordinates along the great arc expressed as Cartesian xyz triples.
        The shape of the array is (num, 3).

    References
    ----------
    [1] https://www.mathworks.com/matlabcentral/newsreader/view_thread/277881
    [2] https://en.wikipedia.org/wiki/Great-circle_distance#Vector_version

    """
    x0 = c[0]
    y0 = c[1]
    z0 = c[2]

    # Vector from center to first point
    v1 = np.asarray([a[0] - x0, a[1] - y0, a[2] - z0])

    # Distance of the first point from the center
    r = np.linalg.norm(v1)

    # Vector from center to second point
    v2 = np.asarray([b[0] - x0, b[1] - y0, b[2] - z0])

    # The v3 lies in plane of v1 & v2 and is orthogonal to v1
    v3 = np.cross(np.cross(v1, v2), v1)

    # Ensure that the vector has length r
    v3 = r * v3 / np.linalg.norm(v3)

    # Range through the inner angle between v1 and v2
    inner_angles = np.linspace(0, np.arctan2(np.linalg.norm(np.cross(v1, v2)),
                                             np.dot(v1, v2)), num=num)

    # Calculate the Cartesian locations from the first to second points
    v_xyz = np.zeros(shape=(num, 3))
    for k in range(0, num):
        inner_angle = inner_angles[k]
        v_xyz[k, :] = v1 * np.cos(inner_angle) + v3 * np.sin(inner_angle) + c

    return v_xyz
