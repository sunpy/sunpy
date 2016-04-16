"""
SunPy specific representations.

This submodule extends `astropy.coordinates.representations` with two Spherical
representation classes, primarily for use with the Helioprojective coordinate
system, due to the convention of Longitude going from -180 to 180 degrees.
"""

from __future__ import absolute_import, division

from astropy.utils.compat.odict import OrderedDict
from astropy import units as u
from astropy.coordinates.representation import (SphericalRepresentation,
                                                UnitSphericalRepresentation)
from astropy.coordinates import Longitude, Latitude

__all__ = ['Longitude180', 'SphericalWrap180Representation',
           'UnitSphericalWrap180Representation']


class Longitude180(Longitude):
    """
    Quantity class that represents Longitude.
    It sets the default wrap_angle to 180 degrees.

    Parameters
    ----------
    angle : array, list, scalar, `~astropy.units.Quantity`,
       :class:`~astropy.coordinates.Angle` The angle value(s). If a tuple,
       will be interpreted as ``(h, m s)`` or ``(d, m, s)`` depending
       on ``unit``. If a string, it will be interpreted following the
       rules described for :class:`~astropy.coordinates.Angle`.

       If ``angle`` is a sequence or array of strings, the resulting
       values will be in the given ``unit``, or if `None` is provided,
       the unit will be taken from the first given value.

    unit : :class:`~astropy.units.UnitBase`, str, optional
       The unit of the value specified for the angle.  This may be
       any string that `~astropy.units.Unit` understands, but it is
       better to give an actual unit object.  Must be an angular
       unit.

    wrap_angle : :class:`~astropy.coordinates.Angle` or equivalent, or None
       Angle at which to wrap back to ``wrap_angle - 180 deg``.
       If ``None`` (default), it will be taken to be 180 deg unless ``angle``
       has a ``wrap_angle`` attribute already (i.e., is a ``Longitude``),
       in which case it will be taken from there.
    """

    def __new__(cls, angle, unit=None, wrap_angle=180 * u.deg, **kwargs):
        self = super(Longitude180, cls).__new__(cls,
                                                angle,
                                                unit=unit,
                                                wrap_angle=wrap_angle,
                                                **kwargs)
        return self


class UnitSphericalWrap180Representation(UnitSphericalRepresentation):
    """
    Representation of points in 3D Spherical coordinates.
    This representation allows for a negative Longitude.
    It does so by setting wrap_angle=180 degrees.

    Parameters
    ----------
    lon, lat: `~astropy.units.Quantity`
        The longitude and latitude of the point(s) in angular units. The
        latitude should be between -90 and +90 degrees, and the longitude
        is allowed to have any value between -180 to 180 degrees. These
        can also be instances of `~astropy.units.Angle`,
        `~astropy.coordinates.Longitude`, or `~astropy.coordinates.Latitude`.

    distance: `~astropy.units.Quantity`
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, else
        it is passed to the :class:`~astropy.units.Quantity` class.

    copy: bool, optional
        If True, arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude180), ('lat', Latitude)])
    recommended_units = {'lon': u.deg, 'lat': u.deg}


class SphericalWrap180Representation(SphericalRepresentation):
    """
    Representation of points in 3D Spherical coordinates.
    This representation allows for a negative Longitude.
    It does so by setting wrap_angle=180 degrees.

    Parameters
    ----------
    lon, lat: `~astropy.units.Quantity`
        The longitude and latitude of the point(s) in angular units. The
        latitude should be between -90 and +90 degrees, and the longitude
        is allowed to have any value between -180 to 180 degrees. These
        can also be instances of `~astropy.units.Angle`,
        `~astropy.coordinates.Longitude`, or `~astropy.coordinates.Latitude`.

    distance: `~astropy.units.Quantity`
        The distance to the point(s). If the distance is a length, it is
        passed to the :class:`~astropy.coordinates.Distance` class, else
        it is passed to the :class:`~astropy.units.Quantity` class.

    copy: bool, optional
        If True, arrays will be copied rather than referenced.
    """

    attr_classes = OrderedDict([('lon', Longitude180), ('lat', Latitude),
                                ('distance', u.Quantity)])
    recommended_units = {'lon': u.deg, 'lat': u.deg}

    _unitrep = UnitSphericalWrap180Representation
