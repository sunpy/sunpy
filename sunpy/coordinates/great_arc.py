#
# Calculates the co-ordinates along great arcs between two specified points
# which are assumed to be on disk.
#
from __future__ import absolute_import, division, print_function

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.coordinates


def great_arc(start, end, center=None, number_points=100):
    """
    Calculate a user-specified number of points on a great arc between a start
    and end point on a sphere.

    Parameters
    ----------
    start : `~astropy.coordinates.SkyCoord`
        Start point.

    end : `~astropy.coordinates.SkyCoord`
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
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> import sunpy.coordinates
    >>> from sunpy.coordinates import great_arc
    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE
    >>> m = sunpy.map.Map(AIA_171_IMAGE)
    >>> a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=m.coordinate_frame)
    >>> b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=m.coordinate_frame)
    >>> v = great_arc(a, b)
    """

    # Create a helper object that contains all the information we need.
    gc = GreatArcConvertToCartesian(start, end, center)

    # Calculate the points along the great arc.
    great_arc_points_cartesian = calculate_great_arc(gc.start_cartesian, gc.end_cartesian, gc.center_cartesian, number_points)*gc.start_unit

    # Transform the great arc back into the input frame.
    return SkyCoord(great_arc_points_cartesian[:, 0],
                    great_arc_points_cartesian[:, 1],
                    great_arc_points_cartesian[:, 2],
                    frame='heliocentric', observer=gc.observer).transform_to(gc.start_frame)


def great_arc_distance(start, end, center=None):
    """
    Calculate the distance between the start point and end point on a sphere.

    Parameters
    ----------
    start : `~astropy.coordinates.SkyCoord`
        Start point.

    end : `~astropy.coordinates.SkyCoord`
        End point.

    center : `~astropy.coordinates.SkyCoord`
        Center of the sphere.

    Returns
    -------
    distance : `astropy.units.Quantity`
        The distance between the two points on the sphere.
    """
    gc = GreatArcConvertToCartesian(start, end, center)
    this_great_arc = CalculateGreatArcCartesian(gc.start_cartesian, gc.end_cartesian, gc.center_cartesian)
    return this_great_arc.distance() * gc.start_unit


def great_arc_angular_distance(start, end, center=None):
    """
    Calculate the angular distance between the start point and end point on a
    sphere.

    Parameters
    ----------
    start : `~astropy.coordinates.SkyCoord`
        Start point.

    end : `~astropy.coordinates.SkyCoord`
        End point.

    center : `~astropy.coordinates.SkyCoord`
        Center of the sphere.

    Returns
    -------

    """
    gc = GreatArcConvertToCartesian(start, end, center)
    this_great_arc = CalculateGreatArcCartesian(gc.start_cartesian, gc.end_cartesian, gc.center_cartesian)
    return np.rad2deg(this_great_arc.inner_angle) * u.degree


def calculate_great_arc(start_cartesian, end_cartesian, center_cartesian, number_points):
    """
    Calculate a user-specified number of points on a great arc between a start
    and end point on a sphere where the start and end points are assumed to be
    x,y,z Cartesian triples on a sphere relative to a center.

    Parameters
    ----------
    start_cartesian : `~numpy.ndarray`
        Start point expressed as a Cartesian xyz triple.

    end_cartesian : `~numpy.ndarray`
        End point expressed as a Cartesian xyz triple.

    center_cartesian : `~numpy.ndarray`
        Center of the sphere expressed as a Cartesian xyz triple

    number_points : int
        Number of points along the great arc.

    Returns
    -------
    arc : `~numpy.ndarray`
        Co-ordinates along the great arc expressed as Cartesian xyz triples.
        The shape of the array is (num, 3).
    """
    this_great_arc = CalculateGreatArcCartesian(start_cartesian, end_cartesian, center_cartesian)

    # Range through the inner angle between v1 and v2
    inner_angles = np.linspace(0, this_great_arc.inner_angle, num=number_points).reshape(number_points, 1)

    # Calculate the Cartesian locations from the first to second points
    return this_great_arc.v1[np.newaxis, :] * np.cos(inner_angles) + \
           this_great_arc.v3[np.newaxis, :] * np.sin(inner_angles) + \
           this_great_arc.center


class CalculateGreatArcCartesian:
    def __init__(self, start_cartesian, end_cartesian, center_cartesian):
        """
        Calculate the properties of a great arc between a start point and an
        end point on a sphere.  See the references below for a description of
        the algorithm.

        Parameters
        ----------
        start_cartesian : `~numpy.ndarray`
            Start point expressed as a Cartesian xyz triple.

        end_cartesian : `~numpy.ndarray`
            End point expressed as a Cartesian xyz triple.

        center_cartesian : `~numpy.ndarray`
            Center of the sphere expressed as a Cartesian xyz triple

        References
        ----------
        [1] https://www.mathworks.com/matlabcentral/newsreader/view_thread/277881
        [2] https://en.wikipedia.org/wiki/Great-circle_distance#Vector_version

        """
        self.start_cartesian = start_cartesian
        self.end_cartesian = end_cartesian
        self.center_cartesian = center_cartesian
        self.center = np.asarray(self.center_cartesian)

        # Vector from center to first point
        self.v1 = np.asarray(self.start_cartesian) - self.center

        # Distance of the first point from the center
        self.r = np.linalg.norm(self.v1)

        # Vector from center to second point
        self.v2 = np.asarray(self.end_cartesian) - self.center

        # The v3 vector lies in plane of v1 & v2 and is orthogonal to v1
        self.v3 = np.cross(np.cross(self.v1, self.v2), self.v1)
        self.v3 = self.r * self.v3 / np.linalg.norm(self.v3)

        # Inner angle between v1 and v2 in radians
        self.inner_angle = np.arctan2(np.linalg.norm(np.cross(self.v1, self.v2)),
                                      np.dot(self.v1, self.v2))

    def distance(self):
        """
        Calculate the distance on the sphere between the start point and the end
        point.

        Returns
        -------
        distance : float
            The distance on the sphere between the start point and the end
            point.
        """
        return self.r * self.inner_angle


class GreatArcConvertToCartesian:
    def __init__(self, start, end, center=None):
        """
        A helper class that takes the inputs required to compute a great arc
        and returns Cartesian triples, co-ordinate frame and unit information
        required by other functions.

        Parameters
        ----------
        start : `~astropy.coordinates.SkyCoord`
            Start point.

        end : `~astropy.coordinates.SkyCoord`
            End point.

        center : `~astropy.coordinates.SkyCoord`
            Center of the sphere.
        """
        self.start = start
        self.end = end
        self.center = center

        # Units of the start point
        self.start_unit = self.start.transform_to('heliocentric').cartesian.xyz.unit

        # Co-ordinate frame
        self.start_frame = self.start.frame

        # Observer details
        self.observer = self.start.observer

        if self.center is None:
            self.c = SkyCoord(0*self.start_unit,
                              0*self.start_unit,
                              0*self.start_unit, frame='heliocentric')

        self.start_cartesian = self.start.transform_to('heliocentric').cartesian.xyz.to(self.start_unit).value
        self.end_cartesian = self.end.transform_to('heliocentric').cartesian.xyz.to(self.start_unit).value
        self.center_cartesian = self.c.transform_to('heliocentric').cartesian.xyz.to(self.start_unit).value
