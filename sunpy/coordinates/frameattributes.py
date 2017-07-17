# -*- coding: utf-8 -*-
from __future__ import absolute_import, division

import inspect
import datetime
import warnings

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import TimeAttribute, CoordinateAttribute

import sunpy.sun
from sunpy.extern import six
from sunpy.time import parse_time
from sunpy.util.exceptions import SunpyUserWarning

__all__ = ['TimeFrameAttributeSunPy', 'ObserverCoordinateAttribute']


class TimeFrameAttributeSunPy(TimeAttribute):
    """
    Frame attribute descriptor for quantities that are Time objects.
    See the `~astropy.coordinates.Attribute` API doc for further
    information.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value
        if ``default`` is ``None`` and no value was supplied during initialization.

    Returns
    -------
    frame_attr : descriptor
        A new data descriptor to hold a frame attribute
    """

    def convert_input(self, value):
        """
        Convert input value to a Time object and validate by running through the
        Time constructor.  Also check that the input was a scalar.

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """
        if value is None:
            return None, False

        elif value == 'now':
            return Time(datetime.datetime.now()), True

        elif isinstance(value, Time):
            out = value
            converted = False

        elif isinstance(value, six.string_types):
            try:
                out = Time(parse_time(value))
            except Exception as err:
                raise ValueError('Invalid time input {0}={1!r}\n{2}'.format(self.name, value, err))
            converted = True
        else:
            try:
                out = Time(value)
            except Exception as err:
                raise ValueError('Invalid time input {0}={1!r}\n{2}'.format(self.name, value, err))
            converted = True

        if not out.isscalar:
            raise ValueError('Time input {0}={1!r} must be a single (scalar) value'
                             .format(self.name, value))

        return out, converted


class ObserverCoordinateAttribute(CoordinateAttribute):
    """
    An Attribute to describe the location of the observer in the solar system.
    The observer location can be given as a string of a known observer, a
    low-level frame class *or* a `~astropy.coordinates.SkyCoord`, but will
    always be converted to the low-level frame class when accessed.

    Parameters
    ----------
    frame : a coordinate frame class
        The type of frame this attribute can be
    default : object
        Default value for the attribute if not provided
    fallback_coordinate: `astropy.coordinates.BaseCoordinateFrame`
        The coordinate value to return if one can not be calculated.
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    def __init__(self, frame, default=None, fallback_coordinate=None, secondary_attribute=''):
        if fallback_coordinate is None:
            # Import here to prevent circular import
            from .frames import HeliographicStonyhurst

            self.fallback_coordinate = HeliographicStonyhurst(0 * u.deg,
                                                              0 * u.deg,
                                                              1 * u.AU)

        super(ObserverCoordinateAttribute, self).__init__(frame, default=default,
                                                          secondary_attribute=secondary_attribute)

    def convert_input(self, value):
        # If we are reaching here, we do not have an instance, so we fall back
        # to the default location when obstime is not set.
        if isinstance(value, six.string_types):
            return self.fallback_coordinate, True
        else:
            return super(ObserverCoordinateAttribute, self).convert_input(value)

    def _convert_string_to_coord(self, out, instance):
        """
        Given a value and and frame instance calculate the position of the
        object given as a string.
        """

        # Import here to prevent circular import
        from .frames import HeliographicStonyhurst

        obstime = getattr(instance, 'obstime')

        # If no time assume nothing
        if obstime is None:
            warnings.warn("Can not compute location of object '{}'"
                          " without obstime attribute being set.".format(out),
                          SunpyUserWarning, stacklevel=4)

            # If obstime is not set, we can't work out where an object is.
            return self.fallback_coordinate

        if out == 'earth':
            # Import here to prevent circular import
            from .frames import HeliographicStonyhurst

            distance = sunpy.sun.sunearth_distance(obstime)
            lon = 0 * u.deg
            lat = sunpy.sun.heliographic_solar_center(obstime)[1]
            return HeliographicStonyhurst(lon=lon, lat=lat, radius=distance)

        else:
            raise ValueError("Only Earth is currently supported as a known observer location.")

    def __get__(self, instance, frame_cls=None):
        # If instance is None then we can't get obstime so it doesn't matter.
        if instance is not None:
            out = getattr(instance, '_' + self.name, self.default)
            print(out)

            # Convert strings to coordinates
            if isinstance(out, six.string_types):
                out = out.lower()
                out = self._convert_string_to_coord(out, instance)
                setattr(instance, '_' + self.name, out)

        return super(ObserverCoordinateAttribute, self).__get__(instance, frame_cls=frame_cls)
