import astropy.time
from astropy.time import TimeDeltaFormat, TimeUnique
from astropy import _erfa as erfa

import numpy as np

from time import strftime, strptime
from datetime import date, datetime, timedelta

__all__ = ['Time', 'TimeDeltaDatetime']

STANDARD_TIME_SCALES = ('tai', 'tcb', 'tcg', 'tdb', 'tt', 'ut1', 'utc')
LOCAL_SCALES = ('local', )
TIME_TYPES = dict(
    (scale, scales) for scales in (STANDARD_TIME_SCALES, LOCAL_SCALES) for scale in scales)

if hasattr(astropy.time.Time, 'strptime'):
    from astropy.time import Time
else:
    class Time(astropy.time.Time):
        def __init__(self, val, val2=None, format=None, scale=None,
                     precision=None, in_subfmt=None, out_subfmt=None,
                     location=None, copy=False):

            if location is not None:
                from ..coordinates import EarthLocation
                if isinstance(location, EarthLocation):
                    self.location = location
                else:
                    self.location = EarthLocation(*location)
            else:
                self.location = None

            if isinstance(val, self.__class__):
                # Update _time formatting parameters if explicitly specified
                if precision is not None:
                    self._time.precision = precision
                if in_subfmt is not None:
                    self._time.in_subfmt = in_subfmt
                if out_subfmt is not None:
                    self._time.out_subfmt = out_subfmt
                self.SCALES = TIME_TYPES[self.scale]
                if scale is not None:
                    self._set_scale(scale)
            else:
                self._init_from_vals(val, val2, format, scale, copy,
                                     precision, in_subfmt, out_subfmt)
                self.SCALES = TIME_TYPES[self.scale]

            if self.location is not None and (self.location.size > 1 and
                                              self.location.shape != self.shape):
                try:
                    # check the location can be broadcast to self's shape.
                    self.location = np.broadcast_to(self.location, self.shape,
                                                    subok=True)
                except Exception:
                    raise ValueError('The location with shape {0} cannot be '
                                     'broadcast against time with shape {1}. '
                                     'Typically, either give a single location or '
                                     'one for each time.'
                                     .format(self.location.shape, self.shape))

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
            format_string : str
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

            to_string = (str if time_array.dtype.kind == 'U' else
                         lambda x: str(x.item(), encoding='ascii'))
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
            format_spec : str
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
                datetime_tuple = (sk['year'], sk['mon'], sk['day'], sk['hour'], sk['min'],
                                  sk['sec'], date_tuple[6], date_tuple[7], -1)
                formatted_strings.append(strftime(format_spec, datetime_tuple))

            if self.isscalar:
                return formatted_strings[0]
            else:
                return np.array(formatted_strings).reshape(self.shape)

        def __hash__(self):
            return hash((self.jd1, self.jd2, self.scale))


class TimeDeltaDatetime(TimeDeltaFormat, TimeUnique):
    """Time delta in datetime.timedelta"""
    name = 'datetime'

    def _check_val_type(self, val1, val2):
        # Note: don't care about val2 for this class
        if not all(isinstance(val, timedelta) for val in val1.flat):
            raise TypeError('Input values for {0} class must be '
                            'datetime.timedelta objects'.format(self.name))
        return val1, None

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        iterator = np.nditer([val1, None], flags=['refs_ok'], op_dtypes=[object] + [np.double])

        for val, sec in iterator:
            sec[...] = val.item().total_seconds()

        self.jd1, self.jd2 = astropy.time.utils.day_frac(
            iterator.operands[-1], 0.0, divisor=erfa.DAYSEC)

    @property
    def value(self):
        iterator = np.nditer(
            [self.jd1 + self.jd2, None], flags=['refs_ok'], op_dtypes=[self.jd1.dtype] + [object])

        for jd, out in iterator:
            out[...] = timedelta(days=jd.item())

        return iterator.operands[-1]
