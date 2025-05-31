"""
Screen class definitions for making assumptions about off-disk emission
"""
import abc
import logging

import numpy as np

import astropy.units as u
from astropy.coordinates.representation import CartesianRepresentation, UnitSphericalRepresentation

from sunpy import log
from sunpy.coordinates import HeliographicStonyhurst, Helioprojective, _transformations
from sunpy.util.decorators import ACTIVE_CONTEXTS
from sunpy.util.exceptions import warn_user

__all__ = ['BaseScreen', 'SphericalScreen', 'PlanarScreen']


class BaseScreen(abc.ABC):

    def __init__(self, only_off_disk=False):
        self.only_off_disk = only_off_disk

    @property
    def _context_name(self):
        # Used for setting active context
        return f"{self.__module__}.{self.__class__.__qualname__}"

    @abc.abstractmethod
    @u.quantity_input
    def calculate_distance(self) -> u.cm:
        ...

    def __enter__(self):
        ACTIVE_CONTEXTS.append(self._context_name)
        self._old_assumed_screen = Helioprojective._assumed_screen  # nominally None
        Helioprojective._assumed_screen = self

    def __exit__(self, exc_type, exc_value, exc_tb):
        if (removed := ACTIVE_CONTEXTS.pop()) != self._context_name:
            raise RuntimeError(f"Cannot remove {self._context_name} from tracking stack because {removed} is last active.")
        Helioprojective._assumed_screen = self._old_assumed_screen

    def _iterate_calculate_distance(self, coord, distance, screen_frame):
        """
        Numerically calculates the distance component to a screen to promote a
        coordinate from 2D to 3D.

        This method should be used only if the screen has been differentially rotated,
        since otherwise promoting a coordinate from 2D to 3D is straightforward.

        The approach used for numerical solving is to iterate to the solution starting
        with an initial guess. For each guess, we construct a 3D coordinate and perform
        a 3D->2D->3D transformation in the native frame of the screen. The amount that
        the coordinate shifts in the 3D->2D->3D transformation is related to how far off
        the guess is: when the shift is zero, the guess is correct. To form the next
        guess, we linearly interpolate/extrapolate between the latest guess and the
        closest guess over past iterations to estimate where the shift will be zero.

        The maximum number of iterations is hard-coded to 20, but convergence will
        typically take no more than 8 iterations. When working with an array coordinate,
        elements that have converged are skipped in subsequent iterations.
        """
        # Check if the logging level is at least DEBUG (for performance reasons)
        debug_output = log.getEffectiveLevel() <= logging.DEBUG

        log.debug("Differentially rotating the screen")
        use = np.ones(distance.shape, dtype=bool)
        delta = np.empty_like(distance)
        for niter in range(20):
            # If this is not the first iteration, update best distance and calculate new distance
            if niter > 0:
                if niter == 1:
                    # For the second iteration, use a simple guess
                    best_distance, distance = distance, distance + delta
                    best_delta = delta.copy()
                else:
                    last_distance = distance.copy()
                    # Starting with the third iteration, linearly interpolate between the last distance and the best distance
                    with np.errstate(invalid='ignore'):
                        distance[use] = last_distance[use] - delta[use] * np.nan_to_num((best_distance[use] - last_distance[use]) / (best_delta[use] - delta[use]))

                    # Update the best distance (not including the most recent calculation)
                    best_distance = np.where(np.abs(best_delta) < np.abs(delta), best_distance, last_distance)
                    best_delta = np.where(np.abs(best_delta) < np.abs(delta), best_delta, delta)

            # Calculate the corresponding delta of a 3D->2D->3D transformation in the native frame of the screen
            # If delta is zero in that frame, then the 3D point is exactly on the screen, so the guessed distance is correct
            other_3d = coord[use].realize_frame(coord[use].represent_as(UnitSphericalRepresentation) * distance[use])
            native_3d = other_3d.transform_to(screen_frame)
            native_2d = native_3d.realize_frame(native_3d.represent_as(UnitSphericalRepresentation))
            delta[use] = self.calculate_distance(native_2d) - native_3d.spherical.distance

            # Continue computations on only those elements that do not meet the tolerance
            # A tolerance of 1e-11 is larger than numerical-precision errors (~1e-12), and only 1.5 meters at 1 AU
            with np.errstate(invalid='ignore'):
                use = np.abs(delta / distance) >= 1e-11*u.one

            if debug_output:
                argmax = np.unravel_index(np.argmax(np.abs(delta)), delta.shape)
                log.debug("%d points remaining; largest delta is %s at distance %s", np.sum(use), delta[argmax], distance[argmax])

            # Return the solution if all elements have converged
            if np.sum(use) == 0:
                log.debug("Solved for the differentially rotated screen after %d iterations", niter + 1)
                return distance

            # Squash small deltas to avoid warnings from later computations
            delta = np.where(use, delta, 0)

        warn_user(f"Failed to solve for the differentially rotated screen after {niter + 1} iterations. "
                  "Using the best guess.")
        return distance


class SphericalScreen(BaseScreen):
    """
    Context manager to interpret 2D coordinates as being on the inside of a spherical screen.

    The radius of the screen is the distance between the specified ``center`` and Sun center.
    This ``center`` does not have to be the same as the observer location for the coordinate
    frame. If they are the same, then this context manager is equivalent to assuming that the
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
        still mapped onto the surface of the Sun. Defaults to `False`.

    See Also
    --------
    sunpy.coordinates.PlanarScreen

    Notes
    -----
    If this context manager is combined with the :func:`~sunpy.coordinates.propagate_with_solar_surface`
    context manager, significantly more computations are required to numerically solve for the
    distance to the differentially rotated screen.

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

        # Iterate the calculation if differential rotation is being applied
        if _transformations._autoapply_diffrot is not None and np.any(frame.obstime != self._center.obstime):
            screen_frame = Helioprojective(observer=self._center, obstime=self._center.obstime)
            distance = self._iterate_calculate_distance(frame, distance, screen_frame)

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
        still mapped onto the surface of the Sun. Defaults to `False`.

    See Also
    --------
    sunpy.coordinates.SphericalScreen

    Notes
    -----
    If this context manager is combined with the :func:`~sunpy.coordinates.propagate_with_solar_surface`
    context manager, significantly more computations are required to numerically solve for the
    distance to the differentially rotated screen.

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
        obs_to_vantage = self._vantage_point.transform_to(frame).cartesian
        vantage_to_sun = CartesianRepresentation(1, 0, 0) * frame.observer.radius - obs_to_vantage
        direction = vantage_to_sun / vantage_to_sun.norm()
        d_from_plane = frame.observer.radius * direction.x - self._distance_from_center
        rep = frame.represent_as(UnitSphericalRepresentation)
        distance = d_from_plane / rep.dot(direction)

        # Iterate the calculation if differential rotation is being applied
        if _transformations._autoapply_diffrot is not None and np.any(frame.obstime != self._vantage_point.obstime):
            screen_frame = Helioprojective(observer=self._vantage_point, obstime=self._vantage_point.obstime)
            distance = self._iterate_calculate_distance(frame, distance, screen_frame)

        return distance
