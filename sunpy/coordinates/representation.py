import inspect
from collections import OrderedDict

import numpy as np

import astropy.units as u
from astropy.coordinates.representation import (BaseRepresentation,
                                                CartesianRepresentation,
                                                SphericalRepresentation,
                                                UnitSphericalRepresentation)
from astropy.coordinates.angles import Angle, Longitude
from astropy.coordinates.distances import Distance

__all__ = ['UnitSouthPoleSphericalRepresentation', 'SouthPoleSphericalRepresentation']


class UnitSouthPoleSphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates where the latitude
    coordinate is measured as inclination from the *south* pole. This

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity` or str
        The azimuth and inclination of the point(s), in angular units. The
        inclination should be between 0 and 180 degrees, and the azimuth will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`.  If ``copy`` is False, `lon`
        will be changed inplace if it is not between 0 and 360 degrees.

    differentials : dict, `PhysicsSphericalDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single
        `PhysicsSphericalDifferential` instance, or a dictionary of of
        differential instances with keys set to a string representation of the
        SI unit with which the differential (derivative) is taken. For example,
        for a velocity differential on a positional representation, the key
        would be ``'s'`` for seconds, indicating that the derivative is a time
        derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Angle)])

    def __init__(self, lon, lat, differentials=None, copy=True):
        super().__init__(lon, lat, copy=copy, differentials=differentials)

        if np.any(self._lat < 0.*u.deg) or np.any(self._lat > 180.*u.deg):
            raise ValueError('Inclination angle(s) must be within '
                             '0 deg <= angle <= 180 deg, '
                             'got {0}'.format(lat.to(u.degree)))

    @property
    def lon(self):
        """
        The azimuth of the point(s).
        """
        return self._lon

    @property
    def lat(self):
        """
        The elevation of the point(s).
        """
        return self._lat

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is a spherical representation

        # TODO: this could be optimized to shortcut even if a differential_class
        # is passed in, using the ._re_represent_differentials() method
        if inspect.isclass(other_class) and not differential_class:
            if issubclass(other_class, SphericalRepresentation):
                return other_class(lon=self.lon, lat=self.lat - 90 * u.deg,
                                   distance=1.0)
            elif issubclass(other_class, UnitSphericalRepresentation):
                return other_class(lon=self.lon, lat=self.lat - 90 * u.deg)

        return super().represent_as(other_class, differential_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """

        x = np.sin(self.lat) * np.cos(self.lon)
        y = np.sin(self.lat) * np.sin(self.lon)
        z = -1 * np.cos(self.lat)

        return CartesianRepresentation(x=x, y=y, z=z, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        z = -1 * cart.z
        s = np.hypot(cart.x, cart.y)

        lon = np.arctan2(cart.y, cart.x)
        lat = np.arctan2(s, z)

        return cls(lon=lon, lat=lat, copy=False)

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.  For
        spherical coordinates, this is just the absolute value of the radius.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return u.Quantity(np.ones(self.shape), u.dimensionless_unscaled,
                          copy=False)


class SouthPoleSphericalRepresentation(BaseRepresentation):
    """
    Representation of points in 3D spherical coordinates where the latitude
    coordinate is measured as inclination from the *south* pole. This

    Parameters
    ----------
    lon, lat : `~astropy.units.Quantity` or str
        The azimuth and inclination of the point(s), in angular units. The
        inclination should be between 0 and 180 degrees, and the azimuth will
        be wrapped to an angle between 0 and 360 degrees. These can also be
        instances of `~astropy.coordinates.Angle`.  If ``copy`` is False, `lon`
        will be changed inplace if it is not between 0 and 360 degrees.

    r : `~astropy.units.Quantity`
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, otherwise
        it is passed to the :class:`~astropy.units.Quantity` class.

    differentials : dict, `PhysicsSphericalDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single
        `PhysicsSphericalDifferential` instance, or a dictionary of of
        differential instances with keys set to a string representation of the
        SI unit with which the differential (derivative) is taken. For example,
        for a velocity differential on a positional representation, the key
        would be ``'s'`` for seconds, indicating that the derivative is a time
        derivative.

    copy : bool, optional
        If `True` (default), arrays will be copied rather than referenced.
    """

    _unit_representation = UnitSouthPoleSphericalRepresentation

    attr_classes = OrderedDict([('lon', Longitude),
                                ('lat', Angle),
                                ('distance', u.Quantity)])

    def __init__(self, lon, lat, distance, differentials=None, copy=True):
        super().__init__(lon, lat, distance, copy=copy, differentials=differentials)

        if np.any(self._lat < 0.*u.deg) or np.any(self._lat > 180.*u.deg):
            raise ValueError('Inclination angle(s) must be within '
                             '0 deg <= angle <= 180 deg, '
                             'got {0}'.format(lat.to(u.degree)))

        if self._distance.unit.physical_type == 'length':
            self._distance = self.distance.view(Distance)

    @property
    def lon(self):
        """
        The azimuth of the point(s).
        """
        return self._lon

    @property
    def lat(self):
        """
        The elevation of the point(s).
        """
        return self._lat

    @property
    def distance(self):
        """
        The distance from the origin to the point(s).
        """
        return self._distance

    # def unit_vectors(self):
    #     sinlon, coslon = np.sin(self.lon), np.cos(self.lon)
    #     sinlat, coslat = np.sin(self.lat), np.cos(self.lat)
    #     return OrderedDict(
    #         (('lon', CartesianRepresentation(-sinlon, coslon, 0., copy=False)),
    #          ('lat', CartesianRepresentation(coslat*coslon,
    #                                            coslat*sinlon,
    #                                            -sinlat, copy=False)),
    #          ('r', CartesianRepresentation(sinlat*coslon, sinlat*sinlon,
    #                                        coslat, copy=False))))

    # def scale_factors(self):
    #     r = self.distance / u.radian
    #     sinlat = np.sin(self.lat)
    #     l = np.broadcast_to(1.*u.one, self.shape, subok=True)
    #     return OrderedDict((('lon', r * sinlat),
    #                         ('lat', r),
    #                         ('r', l)))

    def represent_as(self, other_class, differential_class=None):
        # Take a short cut if the other class is a spherical representation

        # TODO: this could be optimized to shortcut even if a differential_class
        # is passed in, using the ._re_represent_differentials() method
        if inspect.isclass(other_class) and not differential_class:
            if issubclass(other_class, SphericalRepresentation):
                return other_class(lon=self.lon, lat=self.lat - 90 * u.deg,
                                   distance=1.0)
            elif issubclass(other_class, UnitSphericalRepresentation):
                return other_class(lon=self.lon, lat=self.lat - 90 * u.deg)

        return super().represent_as(other_class, differential_class)

    def to_cartesian(self):
        """
        Converts spherical polar coordinates to 3D rectangular cartesian
        coordinates.
        """

        # We need to convert Distance to Quantity to allow negative values.
        if isinstance(self.distance, Distance):
            d = self.distance.view(u.Quantity)
        else:
            d = self.distance

        x = d * np.sin(self.lat) * np.cos(self.lon)
        y = d * np.sin(self.lat) * np.sin(self.lon)
        z = -d * np.cos(self.lat)

        return CartesianRepresentation(x=x, y=y, z=z, copy=False)

    @classmethod
    def from_cartesian(cls, cart):
        """
        Converts 3D rectangular cartesian coordinates to spherical polar
        coordinates.
        """

        z = -1 * cart.z
        s = np.hypot(cart.x, cart.y)
        d = np.hypot(s, z)

        lon = np.arctan2(cart.y, cart.x)
        lat = np.arctan2(s, z)

        return cls(lon=lon, lat=lat, distance=d, copy=False)

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.  For
        spherical coordinates, this is just the absolute value of the radius.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return np.abs(self.distance)
