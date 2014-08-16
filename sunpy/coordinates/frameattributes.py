# -*- coding: utf-8 -*-

import datetime

from astropy.time import Time
from astropy.coordinates.baseframe import TimeFrameAttribute

from sunpy.time import parse_time

__all__ = [TimeFrameAttributeSunpy]

class TimeFrameAttributeSunPy(TimeFrameAttribute):
    """
    Frame attribute descriptor for quantities that are Time objects.
    See the `~astropy.coordinates.FrameAttribute` API doc for further
    information.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value
        if ``default is None`` and no value was supplied during initialization.

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
        
        if value == 'now':
            return datetime.datetime.now(), True

        if isinstance(value, Time):
            out = value
            converted = False
        
        if isinstance(value, basestring):
            out = Time(parse_time(value))
            converted = True
        else:
            try:
                out = Time(value)
            except Exception as err:
                raise ValueError('Invalid time input {0}={1!r}\n{2}'
                                 .format(self.name, value, err))
            converted = True

        if not out.isscalar:
            raise ValueError('Time input {0}={1!r} must be a single (scalar) value'
                             .format(self.name, value))

        return out, converted