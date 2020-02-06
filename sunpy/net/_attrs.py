"""
Implementation of global attrs.

These are defined in here to keep the `sunpy.net.attrs` namespace clean, and to
prevent circular imports.
"""
import collections

import astropy.units as u

from sunpy.time import parse_time, TimeRange

from .attr import Range, SimpleAttr

from sunpy.util.decorators import add_common_docstring
from sunpy.time.time import _variables_for_parse_time_docstring

__all__ = ['Resolution', 'Detector', 'Sample', 'Level', 'Instrument', 'Wavelength', 'Time']


@add_common_docstring(**_variables_for_parse_time_docstring())
class Time(Range):
    """
    Specify the time range of the query.

    Parameters
    ----------
    start : {parse_time_types}
        The start time in a format parseable by `~sunpy.time.parse_time` or
        a `sunpy.time.TimeRange` object.
    end : {parse_time_types}
        The end time of the range.
    near : {parse_time_types}
        Return a singular record closest in time to this value as possible,
        inside the start and end window. Note: not all providers support this
        functionality.

    """
    def __init__(self, start, end=None, near=None):
        if end is None and not isinstance(start, TimeRange):
            raise ValueError("Specify start and end or start has to be a TimeRange")
        if isinstance(start, TimeRange):
            self.start = start.start
            self.end = start.end
        else:
            self.start = parse_time(start)
            self.end = parse_time(end)

        if self.start > self.end:
            raise ValueError("End time must be after start time.")
        self.near = None if near is None else parse_time(near)

        super().__init__(self.start, self.end)

    def __hash__(self):
        if not (isinstance(self.start, collections.Hashable) and
                isinstance(self.end, collections.Hashable)):
            # The hash is the hash of the start and end time
            return hash((self.start.jd1, self.start.jd2, self.start.scale,
                         self.end.jd1, self.end.jd2, self.end.scale))
        else:
            return super().__hash__()

    def collides(self, other):
        return isinstance(other, self.__class__)

    def __xor__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError
        if self.near is not None or other.near is not None:
            raise TypeError
        return Range.__xor__(self, other)

    def pad(self, timedelta):
        return Time(self.start - timedelta, self.start + timedelta)

    def __repr__(self):
        return '<Time({s.start!r}, {s.end!r}, {s.near!r})>'.format(s=self)


class Wavelength(Range):
    def __init__(self, wavemin, wavemax=None):
        """
        Specifies the wavelength or spectral energy range of the detector.

        Parameters
        ----------
        wavemin : `~astropy.units.Quantity`
            The lower bounds of the range.
        wavemax : `~astropy.units.Quantity`
            The upper bound of the range, if not specified it will default to
            the lower bound.

        Notes
        -----
        The VSO understands the 'wavelength' in one of three units, Angstroms,
        kHz or keV. Therefore any unit which is directly convertible to these
        units is valid input.
        """

        if wavemax is None:
            wavemax = wavemin

        if not all(isinstance(var, u.Quantity) for var in [wavemin, wavemax]):
            raise TypeError("Wave inputs must be astropy Quantities")

        if not all([wavemin.isscalar, wavemax.isscalar]):
            raise ValueError("Both wavemin and wavemax must be scalar values")

        # VSO just accept inputs as Angstroms, kHz or keV, the following
        # converts to any of these units depending on the spectral inputs
        # Note: the website asks for GHz, however it seems that using GHz
        # produces weird responses on VSO.
        supported_units = [u.AA, u.kHz, u.keV]
        for unit in supported_units:
            if wavemin.unit.is_equivalent(unit):
                break
            else:
                unit = None
        if unit is None:
            raise u.UnitsError(f"This unit is not convertable to any of {supported_units}")

        wavemin, wavemax = sorted([wavemin.to(unit), wavemax.to(unit)])
        self.unit = unit

        super().__init__(wavemin, wavemax)

    def collides(self, other):
        return isinstance(other, self.__class__)

    def __repr__(self):
        return "<Wavelength({!r}, {!r}, '{!s}')>".format(self.min.value,
                                                            self.max.value,
                                                            self.unit)


class Instrument(SimpleAttr):
    """
    Specifies the Instrument name for the search.

    Parameters
    ----------
    value : `str`

    Notes
    -----
    More information about each instrument supported by the VSO may be found
    within the VSO Registry. For a list of instruments see
    https://sdac.virtualsolar.org/cgi/show_details?keyword=INSTRUMENT.
    """
    def __init__(self, value):
        if not isinstance(value, str):
            raise ValueError("Instrument names must be strings")

        super().__init__(value)


class Level(SimpleAttr):
    """
    Specifies the data processing level to search for.  The data processing
    level is specified by the instrument PI.  May not work with all archives.

    Parameters
    ----------
    value : `float` or `str`

        The value can be entered in of three ways:

        #. May be entered as a string or any numeric type for equality matching
        #. May be a string of the format '(min) - (max)' for range matching
        #. May be a string of the form '(operator) (number)' where operator is\
        one of: lt gt le ge < > <= >=

    """
    pass


class Sample(SimpleAttr):
    """
    Time interval for data sampling.

    Parameters
    ----------
    value : `astropy.units.Quantity`
        A sampling rate convertible to seconds.
    """
    @u.quantity_input
    def __init__(self, value: u.s):
        super().__init__(value)
        self.value = value.to_value(u.s)


class Detector(SimpleAttr):
    """
    The detector from which the data comes from.

    Parameters
    ----------
    value : `str`
    """
    pass


class Resolution(SimpleAttr):
    """
    Resolution level of the data.

    Parameters
    ----------
    value : `float` or `str`

        The value can be entered in of three ways:

        #. May be entered as a string or any numeric type for equality matching
        #. May be a string of the format '(min) - (max)' for range matching
        #. May be a string of the form '(operator) (number)' where operator is\
        one of: lt gt le ge < > <= >=

        This attribute is currently implemented for SDO/AIA and HMI only.
        The "resolution" is a function of the highest level of data available.
        If the CCD is 2048x2048, but is binned to 512x512 before downlink,
        the 512x512 product is designated as '1'.  If a 2048x2048 and 512x512
        product are both available, the 512x512 product is designated '0.25'.

    References
    ----------
    Documentation in SSWIDL routine vso_search.pro.
    """
    pass
