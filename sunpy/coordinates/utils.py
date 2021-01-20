"""
Miscellaneous utilities related to coordinates
"""

import numpy as np

import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, SkyCoord

from sunpy.coordinates import Heliocentric, get_body_heliographic_stonyhurst

__all__ = ['GreatArc', 'CoordinateVisibility', 'get_rectangle_coordinates', 'solar_angle_equivalency']


class GreatArc:
    """
    Calculate the properties of points on a great arc defined by two
    user-specified coordinates on a sphere.

    Parameters
    ----------
    initial : `~astropy.coordinates.SkyCoord`
        The first coordinate.

    target : `~astropy.coordinates.SkyCoord`
        The second coordinate.

    center : `~astropy.coordinates.SkyCoord`
        Center of the sphere.

    points : `None`, `int`, `numpy.ndarray`
        Number of points along the great arc.  If `None`, coordinates are
        calculated at 100 equally spaced points along the arc.  If `int`,
        coordinates are calculated at "points" equally spaced points along
        the arc.  If a numpy.ndarray is passed, it must be one dimensional
        and have values >=0 and <=1.  The values in this array correspond to
        parameterized locations along the arc, with zero corresponding to
        the initial coordinate and 1 corresponding to the last point of the
        arc. Setting this keyword on initializing a GreatArc object sets the
        locations of the default points along the great arc.

    great_circle : `bool`
        If True, calculate a great circle that passes through the initial and
        target coordinates.  If False, calculate points that lie along an arc
        between the initial and target coordinate.

    use_inner_angle_direction : `bool`
        Defines the direction of the great arc on the sphere. If True, then
        the great arc is directed along the inner angle from the initial to the
        target coordinate. If False, then the great arc is directed along the
        outer angle from the initial to the target coordinate.

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
    >>> ax.plot_coord(great_arc.coordinates, color='c')  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP
    """
    def __init__(self, initial, target, center=None, points=100, great_circle=False, use_inner_angle_direction=True):
        # Initial and target coordinates
        self.initial = initial
        self.target = target

        # Observer
        self._output_observer = self.initial.observer

        # Co-ordinate frame of the initial coordinate is the frame used in the output
        self._output_frame = self.initial.frame

        # Observation time of the initial coordinate is the frame used in the output
        self._output_obstime = self.initial.obstime

        # Initial point of the great arc
        self._initial = self.initial.transform_to(Heliocentric)

        # Target point of the great arc
        self._target = self.target.transform_to(self._output_frame).transform_to(Heliocentric)

        # Units of the initial point
        self._distance_unit = self._initial.cartesian.xyz.unit

        # Set the center of the sphere
        if center is None:
            self._center = SkyCoord(0 * self._distance_unit,
                                    0 * self._distance_unit,
                                    0 * self._distance_unit,
                                    obstime=self._output_obstime,
                                    observer=self._output_observer,
                                    frame=Heliocentric)
        else:
            self._center = center.transform_to(self._output_frame).transform_to(Heliocentric)

        # Interpret the points keyword
        self.points = points
        if isinstance(self.points, int):
            self._points = np.linspace(0, 1, self.points)
        elif isinstance(self.points, np.ndarray):
            if self.points.ndim > 1:
                raise ValueError('One dimensional numpy ndarrays only.')
            if np.any(self.points < 0) or np.any(self.points > 1):
                raise ValueError('All value in points array must be strictly >=0 and <=1.')
            self._points = self.points
        else:
            raise ValueError('Incorrectly specified "points" keyword value.')

        # Did the user ask for a great circle?
        self.great_circle = great_circle

        # Which direction between the initial and target points?
        self.use_inner_angle_direction = use_inner_angle_direction

        # Convert the initial, target and center points to their Cartesian values
        self._initial_cartesian = self._initial.cartesian.xyz.to(self._distance_unit).value
        self._target_cartesian = self._target.cartesian.xyz.to(self._distance_unit).value
        self._center_cartesian = self._center.cartesian.xyz.to(self._distance_unit).value

        # Great arc properties calculation
        # Vector from center to first point
        self._v1 = self._initial_cartesian - self._center_cartesian

        # Distance of the first point from the center
        self._r = np.linalg.norm(self._v1)

        # Vector from center to second point
        self._v2 = self._target_cartesian - self._center_cartesian

        # The v3 vector lies in plane of v1 & v2 and is orthogonal to v1
        self._v3 = np.cross(np.cross(self._v1, self._v2), self._v1)
        self._v3 = self._r * self._v3 / np.linalg.norm(self._v3)

        # Radius of the sphere
        self._radius = self._r * self._distance_unit

        # Calculate the angle subtended by the requested arc
        if self.great_circle:
            full_circle = 2 * np.pi * u.rad
            if self.use_inner_angle_direction:
                self._angle = full_circle
            else:
                self._angle = -full_circle
        else:
            # Inner angle between v1 and v2 in radians
            inner_angle = np.arctan2(np.linalg.norm(np.cross(self._v1, self._v2)), np.dot(self._v1, self._v2)) * u.rad
            if self.use_inner_angle_direction:
                self._angle = inner_angle
            else:
                self._angle = inner_angle - 2 * np.pi * u.rad

        # Distance on the sphere between the initial point and the target point.
        self._distance = self._radius * self._angle.value

        # Interpret the points keyword
        self.points = points
        if isinstance(self.points, int):
            self._points = np.linspace(0, 1, self.points)
        elif isinstance(self.points, np.ndarray):
            if self.points.ndim > 1:
                raise ValueError('One dimensional numpy ndarrays only.')
            if np.any(self.points < 0) or np.any(self.points > 1):
                raise ValueError('All value in points array must be strictly >=0 and <=1.')
            self._points = self.points
        else:
            raise ValueError('Incorrectly specified "points" keyword value.')

        # Angles needed for the coordinate calculation
        self._angles = self._points.reshape(len(self._points), 1)*self._angle.value

    @property
    def angles(self):
        """
        The angles subtended by the arc coordinates relative to the initial
        coordinate.
        """
        return self._points * self._angle

    @property
    def distances(self):
        """
        The distance along the sphere for all the coordinates relative to the
        initial coordinate.
        """
        return self._radius * self.angles.value

    @property
    def coordinates(self):
        """
        Calculates the co-ordinates of the arc, returned in the frame of
        the start coordinate.

        Returns
        -------
        coordinates : `~astropy.coordinates.SkyCoord`
            Co-ordinates along the great arc in the co-ordinate frame of the
            initial coordinate.
        """

        # Calculate the Cartesian locations from the first to second points
        great_arc_points_cartesian = (self._v1[np.newaxis, :] * np.cos(self._angles) +
                                      self._v3[np.newaxis, :] * np.sin(self._angles) +
                                      self._center_cartesian) * self._distance_unit

        # Return the coordinates of the great arc between the start and end
        # points
        return SkyCoord(great_arc_points_cartesian[:, 0],
                        great_arc_points_cartesian[:, 1],
                        great_arc_points_cartesian[:, 2],
                        obstime=self._output_obstime,
                        observer=self._output_observer,
                        frame=Heliocentric).transform_to(self._output_frame)


class CoordinateVisibility:
    def __init__(self, coordinates):
        """
        Calculates the visibility properties of a SkyCoord with respect to the solar disk.
        Although this object can be used with any input SkyCoord its primary purpose is
        calculate the visibility properties of great arcs.

        Parameters
        ----------
        coordinates : `~astropy.coordinates.SkyCoord`
            Visibility properties are calculated for these coordinates.
        """
        self.coordinates = coordinates
        self._visible = self.coordinates.transform_to(Heliocentric).z.value > 0
        self._front = self.visible.astype(np.int)
        self._change = self._front[1:] - self._front[0:-1]

    @property
    def visible(self):
        """
        The visibility of each of the arc coordinates as seen from the arc's
        observer.  If a value is True, then the coordinate is in front of the
        plane of the sky.  If False, then the coordinate is behind the plane of
        the sky.
        """
        return self._visible

    @property
    def all_on_front(self):
        """
        Returns True if every arc coordinate is visible, False otherwise.  When
        True, every arc coordinate is visible from the arc's observer.
        """
        return np.all(self.visible)

    @property
    def all_on_back(self):
        """
        Returns True if every arc coordinate is not visible, False otherwise. When
        True, every coordinate is on the back of the Sun as seen from the arc's
        observer.
        """
        return np.all(~self.visible)

    @property
    def from_front_to_back(self):
        """
        Returns the index at which the visibility of the arc switches from being
        visible (on the disk of the Sun) to not visible (on the back of the Sun),
        as seen by the arc's observer.  The index indicates the first coordinate
        along the direction of the arc that is not visible from the arc's observer.
        If all the arc coordinates are on the front or back of the visible disk of
        the Sun as seen from the arc's viewpoint, None is returned.
        """
        if self.all_on_front or self.all_on_back:
            return None
        else:
            test = self._change == -1
            if np.any(test):
                return np.where(test)[0][0]
            else:
                return None

    @property
    def from_back_to_front(self):
        """
        Returns the index at which the visibility of the arc switches from being
        not visible (on the back of the Sun) to visible (on the front of the Sun),
        as seen by the arc's observer.  The index indicates the first coordinate
        along the direction of the arc that is visible from the arc's observer.
        If all the arc coordinates are on the front or back of the visible disk of
        the Sun as seen from the arc's viewpoint, None is returned.
        """
        if self.all_on_front or self.all_on_back:
            return None
        else:
            test = self._change == 1
            if np.any(test):
                return np.where(test)[0][0]
            else:
                return None


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
    equiv : equivalency function that can be used as keyword ``equivalencies`` for astropy unit conversion.

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
