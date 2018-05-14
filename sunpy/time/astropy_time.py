import astropy.time
import astropy.units as u

import numpy as np

from time import strftime, strptime
from datetime import date, datetime

__all__ = ['Time', 'is_time_equal']

if hasattr(astropy.time.Time, 'strptime'):
    from astropy.time import Time
else:
    class Time(astropy.time.Time):
        @classmethod
        def strptime(cls, time_string, format_string, **kwargs):
            """
            Parse a string to a Time according to a format specification.
            See `time.strptime` documentation for format specification.

            >>> Time.strptime('2012-Jun-30 23:59:60', '%Y-%b-%d %H:%M:%S')
            <Time object: scale='utc' format='isot' value=2012-06-30T23:59:60.000>

            Parameters
            ----------
            time_string : string, sequence, ndarray
                Objects containing time data of type string
            format_string : string
                String specifying format of time_string.
            kwargs : dict
                Any keyword arguments for ``Time``.  If the ``format`` keyword
                argument is present, this will be used as the Time format.

            Returns
            -------
            time_obj : `~astropy.time.Time`
                A new `~astropy.time.Time` object corresponding to the input
                ``time_string``.

            """
            time_array = np.asarray(time_string)

            if time_array.dtype.kind not in ('U', 'S'):
                err = "Expected type is string, a bytes-like object or a sequence"\
                    " of these. Got dtype '{}'".format(time_array.dtype.kind)
                raise TypeError(err)

            to_string = (str
                        if time_array.dtype.kind == 'U' else lambda x: str(x.item(), encoding='ascii'))
            iterator = np.nditer([time_array, None], op_dtypes=[time_array.dtype, 'U30'])

            for time, formatted in iterator:
                # TODO: Make this same as astropy version
                if '%f' in format_string:
                    fstring = str(datetime.strptime(to_string(time), format_string))
                    formatted[...] = fstring.replace(' ', 'T')
                else:
                    time_tuple = strptime(to_string(time), format_string)
                    formatted[...] = '{:04}-{}-{}T{}:{}:{}'.format(*time_tuple)

            format = kwargs.pop('format', None)
            out = cls(*iterator.operands[1:], format='isot', **kwargs)
            if format is not None:
                out.format = format

            return out

        def strftime(self, format_spec):
            """
            Convert Time to a string or a numpy.array of strings according to a
            format specification.
            See `time.strftime` documentation for format specification.

            Parameters
            ----------
            format_spec : string
                Format definition of return string.

            Returns
            -------
            formatted : string, numpy.array
                String or numpy.array of strings formatted according to the given
                format string.

            """
            formatted_strings = []
            for sk in self.replicate('iso')._time.str_kwargs():
                date_tuple = date(sk['year'], sk['mon'], sk['day']).timetuple()
                datetime_tuple = (sk['year'], sk['mon'], sk['day'], sk['hour'], sk['min'], sk['sec'],
                                date_tuple[6], date_tuple[7], -1)
                formatted_strings.append(strftime(format_spec, datetime_tuple))

            if self.isscalar:
                return formatted_strings[0]
            else:
                return np.array(formatted_strings).reshape(self.shape)


def is_time_equal(t1, t2):
    if abs(t1 - t2) < 1 * u.nanosecond:
        return True

    return False
