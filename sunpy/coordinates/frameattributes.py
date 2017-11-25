# -*- coding: utf-8 -*-
from __future__ import absolute_import, division

import datetime
import warnings

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import TimeAttribute, CoordinateAttribute, get_body_barycentric, ICRS

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

        return out, converted


class ObserverCoordinateAttribute(CoordinateAttribute):
    """
    An Attribute to describe the location of the observer in the solar system.
    The observer location can be given as a string of a known observer, which
    will be converted to a coordinate as long as the ``obstime`` attribute is
    valid on the instance of the frame. Alternatively a low-level frame class
    *or* a `~astropy.coordinates.SkyCoord` can be provided to specify the
    location of the observer. If a `~astropy.coordinates.SkyCoord` is passed it
    will always be converted to the low-level frame class when accessed.

    Parameters
    ----------
    frame : a coordinate frame class
        The type of frame this attribute can be
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    def convert_input(self, value):
        # Keep string here.
        if isinstance(value, six.string_types):
            return value, False
        else:
            return super(ObserverCoordinateAttribute, self).convert_input(value)

    def _convert_string_to_coord(self, out, obstime):
        """
        Given a value and and frame instance calculate the position of the
        object given as a string.
        """

        # Import here to prevent circular import
        from .frames import HeliographicStonyhurst
        from .ephemeris import get_body_heliographic_stonyhurst

        obscoord = get_body_heliographic_stonyhurst(out, obstime)

        if out == "earth":
            rep = obscoord.spherical
            rep.lon[()] = 0*u.deg
            obscoord = obscoord.realize_frame(rep)

        return obscoord

    def __get__(self, instance, frame_cls=None):
        # If instance is None then we can't get obstime so it doesn't matter.
        if instance is not None:
            # Get observer if the instance has one, or the default.
            observer = getattr(instance, '_' + self.name, self.default)

            # We have an instance of a frame, so get obstime
            obstime = getattr(instance, 'obstime', None)

            # If the observer is a string and we have obstime then calculate
            # the position of the observer.
            if isinstance(observer, six.string_types):
                if obstime is not None:
                    observer = self._convert_string_to_coord(observer.lower(), obstime)
                    setattr(instance, '_' + self.name, observer)
                else:
                    return observer

        return super(ObserverCoordinateAttribute, self).__get__(instance, frame_cls=frame_cls)
