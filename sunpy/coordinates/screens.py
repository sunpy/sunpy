"""
Screen class definitions for making assumptions about off-disk emission
"""
import abc

import numpy as np

import astropy.units as u
from astropy.coordinates.representation import CartesianRepresentation, UnitSphericalRepresentation

from sunpy.coordinates import HeliographicStonyhurst, Helioprojective
from sunpy.util.decorators import ACTIVE_CONTEXTS

__all__ = ['BaseScreen', 'SphericalScreen', 'PlanarScreen']


class BaseScreen(abc.ABC):

    def __init__(self, only_off_disk=False):
        self.only_off_disk = only_off_disk

    @property
    def _context_name(self):
        # Used for setting active context
        return f'assume_{self.screen_type}_screen'

    @abc.abstractmethod
    @u.quantity_input
    def calculate_distance(self) -> u.cm:
        ...

    def __enter__(self):
        ACTIVE_CONTEXTS[self._context_name] = True
        self._old_assumed_screen = Helioprojective._assumed_screen  # nominally None
        Helioprojective._assumed_screen = self

    def __exit__(self, exc_type, exc_value, exc_tb):
        ACTIVE_CONTEXTS[self._context_name] = False
        Helioprojective._assumed_screen = self._old_assumed_screen


class SphericalScreen(BaseScreen):
    """
    Context manager to interpret 2D coordinates as being on the inside of a spherical screen.

    The radius of the screen is the distance between the specified ``center`` and Sun center.
    This ``center`` does not have to be the same as the observer location for the coordinate
    frame.  If they are the same, then this context manager is equivalent to assuming that the
    helioprojective "zeta" component is zero.
    This replaces the default assumption where 2D coordinates are mapped onto the surface of the
    Sun.

    .. note:: This applies only to coordinates in a `~sunpy.coordinates.Helioprojective` frame.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The center of the spherical screen
    radius : `~astropy.units.Quantity`, optional
        The radius of the spherical screen. The default sets the radius to the distance from the
        screen center to the Sun.
    only_off_disk : `bool`, optional
        If `True`, apply this assumption only to off-disk coordinates, with on-disk coordinates
        still mapped onto the surface of the Sun.  Defaults to `False`.

    See Also
    --------
    sunpy.coordinates.PlanarScreen

    Examples
    --------
    .. minigallery:: sunpy.coordinates.SphericalScreen

    >>> import astropy.units as u
    >>> from sunpy.coordinates import Helioprojective, SphericalScreen
    >>> h = Helioprojective(range(7)*u.arcsec*319, [0]*7*u.arcsec,
    ...                     observer='earth', obstime='2020-04-08')
    >>> print(h.make_3d())
    <Helioprojective Coordinate (obstime=2020-04-08T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        [(   0., 0., 0.99660825), ( 319., 0., 0.99687244),
            ( 638., 0., 0.99778472), ( 957., 0., 1.00103285),
            (1276., 0.,        nan), (1595., 0.,        nan),
            (1914., 0.,        nan)]>

    >>> with SphericalScreen(h.observer):
    ...     print(h.make_3d())
    <Helioprojective Coordinate (obstime=2020-04-08T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        [(   0., 0., 1.00125872), ( 319., 0., 1.00125872),
            ( 638., 0., 1.00125872), ( 957., 0., 1.00125872),
            (1276., 0., 1.00125872), (1595., 0., 1.00125872),
            (1914., 0., 1.00125872)]>

    >>> with SphericalScreen(h.observer, only_off_disk=True):
    ...     print(h.make_3d())
    <Helioprojective Coordinate (obstime=2020-04-08T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        [(   0., 0., 0.99660825), ( 319., 0., 0.99687244),
            ( 638., 0., 0.99778472), ( 957., 0., 1.00103285),
            (1276., 0., 1.00125872), (1595., 0., 1.00125872),
            (1914., 0., 1.00125872)]>
    """
    screen_type = 'spherical'

    @u.quantity_input
    def __init__(self, center, *, radius: u.m=None, **kwargs):
        self._center = center
        if radius is not None:
            self._radius = radius
        else:
            hgs_frame = HeliographicStonyhurst(obstime=self._center.obstime)
            center_hgs = self._center.transform_to(hgs_frame)
            self._radius = center_hgs.radius
        super().__init__(**kwargs)

    def calculate_distance(self, frame):
        sphere_center = self._center.transform_to(frame).cartesian
        c = sphere_center.norm()**2 - self._radius**2
        rep = frame.represent_as(UnitSphericalRepresentation)
        b = -2 * sphere_center.dot(rep)
        # Ignore sqrt of NaNs
        with np.errstate(invalid='ignore'):
            distance = ((-1*b) + np.sqrt(b**2 - 4*c)) / 2  # use the "far" solution
        return distance


class PlanarScreen(BaseScreen):
    """
    Context manager to interpret 2D coordinates as being on the inside of a planar screen.

    The plane goes through Sun center and is perpendicular to the vector between the
    specified vantage point and Sun center.
    This replaces the default assumption where 2D coordinates are mapped onto the surface of the
    Sun.

    .. note:: This applies only to coordinates in a `~sunpy.coordinates.Helioprojective` frame.

    Parameters
    ----------
    vantage_point : `~astropy.coordinates.SkyCoord`
        The vantage point that defines the orientation of the plane.
    distance_from_center : `~astropy.units.Quantity`, optional
        Distance from Sun center of the planar screen. Defaults to 0 such
        that the plane goes through Sun center
    only_off_disk : `bool`, optional
        If `True`, apply this assumption only to off-disk coordinates, with on-disk coordinates
        still mapped onto the surface of the Sun.  Defaults to `False`.

    See Also
    --------
    sunpy.coordinates.SphericalScreen

    Examples
    --------

    >>> import astropy.units as u
    >>> from sunpy.coordinates import Helioprojective, PlanarScreen
    >>> h = Helioprojective(range(7)*u.arcsec*319, [0]*7*u.arcsec,
    ...                     observer='earth', obstime='2020-04-08')
    >>> print(h.make_3d())
    <Helioprojective Coordinate (obstime=2020-04-08T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        [(   0., 0., 0.99660825), ( 319., 0., 0.99687244),
            ( 638., 0., 0.99778472), ( 957., 0., 1.00103285),
            (1276., 0.,        nan), (1595., 0.,        nan),
            (1914., 0.,        nan)]>

    >>> with PlanarScreen(h.observer):
    ...     print(h.make_3d())
    <Helioprojective Coordinate (obstime=2020-04-08T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        [(   0., 0., 1.00125872), ( 319., 0., 1.00125992),
            ( 638., 0., 1.00126351), ( 957., 0., 1.0012695 ),
            (1276., 0., 1.00127788), (1595., 0., 1.00128866),
            (1914., 0., 1.00130183)]>

    >>> with PlanarScreen(h.observer, distance_from_center=1*u.R_sun):
    ...     print(h.make_3d())
    <Helioprojective Coordinate (obstime=2020-04-08T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        [(   0., 0., 0.99660825), ( 319., 0., 0.99660944),
            ( 638., 0., 0.99661302), ( 957., 0., 0.99661898),
            (1276., 0., 0.99662732), (1595., 0., 0.99663805),
            (1914., 0., 0.99665116)]>
    """
    screen_type = 'planar'

    @u.quantity_input
    def __init__(self, vantage_point, *, distance_from_center: u.m=0*u.m, **kwargs):
        self._vantage_point = vantage_point
        self._distance_from_center = distance_from_center
        super().__init__(**kwargs)

    def calculate_distance(self, frame):
        direction = self._vantage_point.transform_to(frame).cartesian
        direction = CartesianRepresentation(1, 0, 0) * frame.observer.radius - direction
        direction /= direction.norm()
        d_from_plane = (frame.observer.radius - self._distance_from_center) * direction.x
        rep = frame.represent_as(UnitSphericalRepresentation)
        distance = d_from_plane / rep.dot(direction)
        return distance
