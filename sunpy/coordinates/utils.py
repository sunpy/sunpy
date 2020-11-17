"""
Miscellaneous utilities related to coordinates
"""

import numpy as np

import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, SkyCoord

from sunpy.coordinates import Heliocentric, get_body_heliographic_stonyhurst

__all__ = ['GreatArc', 'get_rectangle_coordinates', 'solar_angle_equivalency']


class GreatArc:
    """
    Calculate the properties of a great arc at user-specified points between a
    start and end point on a sphere.

    The coordinates of the great arc are returned with the observation time
    and coordinate frame of the starting point of the arc.

    Parameters
    ----------
    start : `~astropy.coordinates.SkyCoord`
        Start point.

    end : `~astropy.coordinates.SkyCoord`
        End point.

    center : `~astropy.coordinates.SkyCoord`
        Center of the sphere.

    points : `None`, `int`, `numpy.ndarray`
        Number of points along the great arc.  If None, the arc is calculated
        at 100 equally spaced points from start to end.  If int, the arc is
        calculated at "points" equally spaced points from start to end.  If a
        numpy.ndarray is passed, it must be one dimensional and have values
        >=0 and <=1.  The values in this array correspond to parameterized
        locations along the great arc from zero, denoting the start of the arc,
        to 1, denoting the end of the arc.  Setting this keyword on initializing
        a GreatArc object sets the locations of the default points along the
        great arc.

    Methods
    -------
    inner_angles : `~astropy.units.Quantity`
        Angles of the points along the great arc from the start to end
        co-ordinate.

    distances : `~astropy.units.Quantity`
        Distances of the points along the great arc from the start to end
        co-ordinate.  The units are defined as those returned after transforming
        the co-ordinate system of the start co-ordinate into its Cartesian
        equivalent.

    coordinates : `~astropy.coordinates.SkyCoord`
        Co-ordinates along the great arc in the co-ordinate frame of the
        start point.

    References
    ----------
    [1] https://www.mathworks.com/matlabcentral/newsreader/view_thread/277881
    [2] https://en.wikipedia.org/wiki/Great-circle_distance#Vector_version

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    >>> from sunpy.coordinates.utils import GreatArc
    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_IMAGE  # doctest: +REMOTE_DATA
    >>> m = sunpy.map.Map(AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=m.coordinate_frame)  # doctest: +REMOTE_DATA
    >>> b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=m.coordinate_frame)  # doctest: +REMOTE_DATA
    >>> great_arc = GreatArc(a, b)  # doctest: +REMOTE_DATA
    >>> ax = plt.subplot(projection=m)  # doctest: +SKIP
    >>> m.plot(axes=ax)  # doctest: +SKIP
    >>> ax.plot_coord(great_arc.coordinates(), color='c')  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

    """

    def __init__(self, start, end, center=None, points=None):

        # Observer
        self.observer = start.observer

        # Co-ordinate frame of the starting point
        self.start_frame = start.frame

        # Observation time
        self.obstime = start.obstime

        # Start point of the great arc
        self.start = start.transform_to(Heliocentric)

        # End point of the great arc
        self.end = end.transform_to(self.start_frame).transform_to(Heliocentric)

        # Parameterized location of points between the start and the end of the
        # great arc.
        # Default parameterized points location.
        self.default_points = np.linspace(0, 1, 100)

        # If the user requests a different set of default parameterized points
        # on initiation of the object, then these become the default.  This
        # allows the user to use the methods without having to specify their
        # choice of points over and over again, while also allowing the
        # flexibility in the methods to calculate other values.
        self.default_points = self._points_handler(points)

        # Units of the start point
        self.distance_unit = self.start.cartesian.xyz.unit

        # Set the center of the sphere
        if center is None:
            self.center = SkyCoord(0 * self.distance_unit,
                                   0 * self.distance_unit,
                                   0 * self.distance_unit,
                                   obstime=self.obstime,
                                   observer=self.observer,
                                   frame=Heliocentric)
        else:
            self.center = center.transform_to(self.start_frame).transform_to(Heliocentric)

        # Convert the start, end and center points to their Cartesian values
        self.start_cartesian = self.start.cartesian.xyz.to(self.distance_unit).value
        self.end_cartesian = self.end.cartesian.xyz.to(self.distance_unit).value
        self.center_cartesian = self.center.cartesian.xyz.to(self.distance_unit).value

        # Great arc properties calculation
        # Vector from center to first point
        self.v1 = self.start_cartesian - self.center_cartesian

        # Distance of the first point from the center
        self._r = np.linalg.norm(self.v1)

        # Vector from center to second point
        self.v2 = self.end_cartesian - self.center_cartesian

        # The v3 vector lies in plane of v1 & v2 and is orthogonal to v1
        self.v3 = np.cross(np.cross(self.v1, self.v2), self.v1)
        self.v3 = self._r * self.v3 / np.linalg.norm(self.v3)

        # Inner angle between v1 and v2 in radians
        self.inner_angle = np.arctan2(np.linalg.norm(np.cross(self.v1, self.v2)),
                                      np.dot(self.v1, self.v2)) * u.rad

        # Radius of the sphere
        self.radius = self._r * self.distance_unit

        # Distance on the sphere between the start point and the end point.
        self.distance = self.radius * self.inner_angle.value

    def _points_handler(self, points):
        """
        Interprets the points keyword.
        """
        if points is None:
            return self.default_points
        elif isinstance(points, int):
            return np.linspace(0, 1, points)
        elif isinstance(points, np.ndarray):
            if points.ndim > 1:
                raise ValueError('One dimensional numpy ndarrays only.')
            if np.any(points < 0) or np.any(points > 1):
                raise ValueError('All value in points array must be strictly >=0 and <=1.')
            return points
        else:
            raise ValueError('Incorrectly specified "points" keyword value.')

    def inner_angles(self, points=None):
        """
        Calculates the inner angles for the parameterized points along the arc
        and returns the value in radians, from the start co-ordinate to the
        end.

        Parameters
        ----------
        points : `None`, `int`, `numpy.ndarray`
            If None, use the default locations of parameterized points along the
            arc.  If int, the arc is calculated at "points" equally spaced
            points from start to end.  If a numpy.ndarray is passed, it must be
            one dimensional and have values >=0 and <=1.  The values in this
            array correspond to parameterized locations along the great arc from
            zero, denoting the start of the arc, to 1, denoting the end of the
            arc.

        Returns
        -------
        inner_angles : `~astropy.units.Quantity`
            Angles of the points along the great arc from the start to
            end co-ordinate.

        """
        these_points = self._points_handler(points)
        return these_points.reshape(len(these_points), 1)*self.inner_angle

    def distances(self, points=None):
        """
        Calculates the distance from the start co-ordinate to the end
        co-ordinate on the sphere for all the parameterized points.

        Parameters
        ----------
        points : `None`, `int`, `numpy.ndarray`
            If None, use the default locations of parameterized points along the
            arc.  If int, the arc is calculated at "points" equally spaced
            points from start to end.  If a numpy.ndarray is passed, it must be
            one dimensional and have values >=0 and <=1.  The values in this
            array correspond to parameterized locations along the great arc from
            zero, denoting the start of the arc, to 1, denoting the end of the
            arc.

        Returns
        -------
        distances : `~astropy.units`
            Distances of the points along the great arc from the start to end
            co-ordinate.  The units are defined as those returned after
            transforming the co-ordinate system of the start co-ordinate into
            its Cartesian equivalent.
        """
        return self.radius * self.inner_angles(points=points).value

    def coordinates(self, points=None):
        """
        Calculates the co-ordinates on the sphere from the start to the end
        co-ordinate for all the parameterized points.  Co-ordinates are
        returned in the frame of the start coordinate.

        Parameters
        ----------
        points : `None`, `int`, `numpy.ndarray`
            If None, use the default locations of parameterized points along the
            arc.  If int, the arc is calculated at "points" equally spaced
            points from start to end.  If a numpy.ndarray is passed, it must be
            one dimensional and have values >=0 and <=1.  The values in this
            array correspond to parameterized locations along the great arc from
            zero, denoting the start of the arc, to 1, denoting the end of the
            arc.

        Returns
        -------
        arc : `~astropy.coordinates.SkyCoord`
            Co-ordinates along the great arc in the co-ordinate frame of the
            start point.

        """
        # Calculate the inner angles
        these_inner_angles = self.inner_angles(points=points)

        # Calculate the Cartesian locations from the first to second points
        great_arc_points_cartesian = (self.v1[np.newaxis, :] * np.cos(these_inner_angles) +
                                      self.v3[np.newaxis, :] * np.sin(these_inner_angles) +
                                      self.center_cartesian) * self.distance_unit

        # Return the coordinates of the great arc between the start and end
        # points
        return SkyCoord(great_arc_points_cartesian[:, 0],
                        great_arc_points_cartesian[:, 1],
                        great_arc_points_cartesian[:, 2],
                        obstime=self.obstime,
                        observer=self.observer,
                        frame=Heliocentric).transform_to(self.start_frame)


def get_rectangle_coordinates(bottom_left, *, top_right=None,
                              width: u.deg = None, height: u.deg = None):
    """
    Specify a rectangular region of interest in longitude and latitude in a given coordinate frame.

    Parameters
    ----------
    bottom_left : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The bottom-left coordinate of the rectangle. Supports passing both the
        bottom left and top right coordinates by passing with a shape of ``(2,)``.
    top_right : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The top-right coordinate of the rectangle.
        If in a different frame than ``bottom_left`` and all required metadata
        for frame conversion is present, ``top_right`` will be transformed to
        ``bottom_left`` frame.
    width : `~astropy.units.Quantity`
        The width of the rectangle.
        Must be omitted if the coordinates of both corners have been specified.
    height : `~astropy.units.Quantity`
        The height of the rectangle.
        Must be omitted if the coordinates of both corners have been specified.

    Returns
    -------
    `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The bottom left coordinate of the rectangular region of interest.
    `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
        The top right coordinate of the rectangular region of interest.

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates.frames import HeliographicStonyhurst
    >>> from sunpy.coordinates.utils import get_rectangle_coordinates

    >>> # With bottom left as a SkyCoord, width and height
    >>> bottom_left = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame='heliographic_stonyhurst')
    >>> width = 10 * u.arcsec
    >>> height = 10 * u.arcsec
    >>> bottom_left, top_right = get_rectangle_coordinates(bottom_left, width=width, height=height)

    >>> # With bottom left of shape (2,)
    >>> bottom_left_vector = SkyCoord([0 * u.arcsec, 10 * u.arcsec], [0 * u.arcsec, 10 * u.arcsec], frame='heliographic_stonyhurst')
    >>> bottom_left, top_right = get_rectangle_coordinates(bottom_left_vector)

    >>> # With bottom left as a BaseCoordinateFrame instance, width and height
    >>> bottom_left =  HeliographicStonyhurst(0 * u.arcsec, 0 * u.arcsec)
    >>> width = 10 * u.arcsec
    >>> height = 10 * u.arcsec
    >>> bottom_left, top_right = get_rectangle_coordinates(bottom_left, width=width, height=height)

    Notes
    -----
    ``width`` is always treated as an increase in longitude, but ``bottom_left`` may have a higher
    value of longitude than ``top_right`` due to the wrapping of the longitude angle.  Appropriate
    care should be taken when using this function's output to define a range of longitudes.

    ``height`` is always treated as an increase in latitude, but this function does not enforce
    that ``bottom_left`` has a lower value of latitude than ``top_right``, in case that orientation
    is valid for the intended use.
    """
    if not (hasattr(bottom_left, 'transform_to') and
            hasattr(bottom_left, 'shape') and
            hasattr(bottom_left, 'spherical')):
        raise TypeError(
            "Invalid input, bottom_left must be of type SkyCoord or BaseCoordinateFrame.")

    if (top_right is not None and not ((hasattr(top_right, 'transform_to') and
                                        hasattr(top_right, 'shape') and
                                        hasattr(top_right, 'spherical')))):
        raise TypeError("Invalid input, top_right must be of type SkyCoord or BaseCoordinateFrame.")

    if bottom_left.shape == (2,) and any((x is not None for x in (width, height, top_right))):
        raise ValueError("Invalid input, if bottom_left.shape == (2,) "
                         "other parameters should not be passed.")

    if all(x is not None for x in (width, height, top_right)):
        raise ValueError("Invalid input, width, height and top_right "
                         "parameters should not be passed simultaneously.")

    if top_right is None and bottom_left.shape != (2,) and (width is None or height is None):
        raise ValueError("Invalid input, either bottom_left and top_right "
                         "or bottom_left and height and width should be provided.")

    if width is not None:
        if width < 0*u.deg:
            raise ValueError("The specified width cannot be negative.")
        if width > 360*u.deg:
            raise ValueError("The specified width cannot be greater than 360 degrees.")

    if height is not None:
        if height < 0*u.deg:
            raise ValueError("The specified height cannot be negative.")
        if bottom_left.spherical.lat + height > 90*u.deg:
            raise ValueError("The specified height exceeds the maximum latitude.")

    if bottom_left.shape == (2,):
        top_right = bottom_left[1]
        bottom_left = bottom_left[0]

    elif top_right is not None:
        top_right = top_right.transform_to(bottom_left)

    else:
        # If bottom left is a ``SkyCoord``, top right is constructed
        # as a SkyCoord using width and height. If bottom left is a
        # ``frame``, the top right is typecasted to its respective
        # frame. This is done to ensure that the output coordinates
        # are of the same type.
        top_right = SkyCoord(bottom_left.spherical.lon + width,
                             bottom_left.spherical.lat + height,
                             frame=bottom_left)

        if isinstance(bottom_left, BaseCoordinateFrame):
            top_right = top_right.frame

    return bottom_left, top_right


def solar_angle_equivalency(observer):
    """
    Return the equivalency to convert between a physical distance on the Sun
    and an angular separation as seen by a specified observer.

    .. note::
        This equivalency assumes that the physical distance is perpendicular to
        the Sun-observer line.  That is, the tangent of the angular separation
        is equal to the ratio of the physical distance to the Sun-observer
        distance.  For large physical distances, a different assumption may be
        more appropriate.

    Parameters
    ----------
    observer : `~astropy.coordinates.SkyCoord`
        Observer location for which the equivalency is calculated.

    Returns
    -------
    equiv : equivalency function that can be used as keyword `equivalencies` for astropy unit conversion.

    Examples
    --------
    >>> import astropy.units as u
    >>> from sunpy.coordinates import get_body_heliographic_stonyhurst
    >>> earth_observer = get_body_heliographic_stonyhurst("earth", "2013-10-28")
    >>> distance_in_km = 725*u.km
    >>> distance_in_km.to(u.arcsec, equivalencies=solar_angle_equivalency(earth_observer))
    INFO: Apparent body location accounts for 495.82 seconds of light travel time [sunpy.coordinates.ephemeris]
    <Quantity 1.00603718 arcsec>
    """

    if not isinstance(observer, (SkyCoord, BaseCoordinateFrame)):
        raise TypeError(
            "Invalid input, observer must be of type SkyCoord or BaseCoordinateFrame.")
    if observer.obstime is None:
        raise ValueError(
            "Observer must have an observation time, `obstime`.")

    obstime = observer.obstime
    sun_coord = get_body_heliographic_stonyhurst("sun", time=obstime, observer=observer)
    sun_observer_distance = sun_coord.separation_3d(observer).to_value(u.m)

    equiv = [(u.radian,
              u.meter,
              lambda x: np.tan(x)*sun_observer_distance,
              lambda x: np.arctan(x/sun_observer_distance))]

    return equiv
