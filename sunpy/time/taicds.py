"""
This module provides support for the "TAI" time format used in
the CDS software
"""
from astropy.time.formats import TimeFromEpoch, erfa

__all__ = ['TimeTAICDS']


class TimeTAICDS(TimeFromEpoch):
    """
    Number of seconds, including leap seconds, since 1958-01-01 00:00:00.

    This format originates from the CDS time conversion software and is the
    format returned by the often-used `anytim2tai` routine in SSW. This format
    and the `anytim2tai` function do vary in their treatment of leap seconds
    before 1972-01-01 00:00:00. For dates between 1958-01-01 00:00:00 UTC
    and 1972-01-01 00:00:00 UTC, `anytim2tai` returns a constant value of 9 s
    difference between UTC and TAI time while this implementation returns a value
    of 0 s difference between UTC and TAI on 1958-01-01 00:00:00, with this difference
    increasing approximately linearly to 10 s on 1972-01-01 00:00:00, where the two
    approaches agree.

    Examples
    --------
    >>> from astropy.time import Time
    >>> t = Time('1958-01-01T00:00:00', format='isot', scale='tai')
    >>> t.tai_cds
    0.0
    >>> t2 = Time('2015-10-25T05:24:08', format='isot', scale='tai')
    >>> t3 = Time(t2.tai_cds, format='tai_cds', scale='tai')
    >>> t3.isot
    '2015-10-25T05:24:08.000'

    References
    ----------
    * `CDS Time Conversion Software README <https://hesperia.gsfc.nasa.gov/ssw/gen/idl/time/aaareadme.txt>`_
    * `anytim2tai routine in SSW <https://hesperia.gsfc.nasa.gov/ssw/gen/idl/time/anytim2tai.pro>`_
    """
    name = 'tai_cds'
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = '1958-01-01 00:00:00'
    epoch_val2 = None
    epoch_scale = 'tai'
    epoch_format = 'iso'
