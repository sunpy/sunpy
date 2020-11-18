"""
Implementation of global attrs.

These are defined in here to keep the `sunpy.net.attrs` namespace clean, and to
prevent circular imports.
"""
import collections.abc

import astropy.units as u

from sunpy.time import TimeRange, parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring
from .attr import Range, SimpleAttr

__all__ = ['Physobs', 'Resolution', 'Detector', 'Sample',
           'Level', 'Instrument', 'Wavelength', 'Time', 'Source', 'Provider']


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
    type_name = "time"

    def __init__(self, start, end=None, near=None):
        if end is None and not isinstance(start, TimeRange):
            raise ValueError("Specify start and end or start has to be a TimeRange")
        if isinstance(start, TimeRange):
            self.start, self.end = start.start, start.end
        else:
            self.start, self.end = parse_time(start), parse_time(end)

        if self.start > self.end:
            raise ValueError("End time must be after start time.")
        self.near = parse_time(near) if near else None

        super().__init__(self.start, self.end)

    def __hash__(self):
        if not isinstance(self.start, collections.abc.Hashable) or \
           not isinstance(self.end, collections.abc.Hashable):
            # The hash is the hash of the start and end time
            return hash((self.start.jd1, self.start.jd2, self.start.scale,
                         self.end.jd1, self.end.jd2, self.end.scale))
        else:
            return super().__hash__()

    def collides(self, other):
        # Use exact type checking here, because otherwise it collides with all
        # subclasses of itself which can have completely different search
        # meanings.
        return type(other) is type(self)

    def __xor__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError
        if self.near is not None or other.near is not None:
            raise TypeError
        return super().__xor__(other)

    def pad(self, timedelta):
        return type(self)(self.start - timedelta, self.start + timedelta)

    def __repr__(self):
        start = self.start.iso
        end = self.end.iso
        iso = self.near.iso if self.near else None
        str_repr = ", ".join(str(param) for param in [start, end, iso] if param)
        return f'<sunpy.net.attrs.Time({str_repr})>'


class Wavelength(Range):
    type_name = 'wave'

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
            raise u.UnitsError(f"This unit is not convertable to any of {supported_units}")

        wavemin, wavemax = sorted([wavemin.to(unit), wavemax.to(unit)])
        self.unit = unit

        super().__init__(wavemin, wavemax)

    def collides(self, other):
        return isinstance(other, self.__class__)

    def __repr__(self):
        return f"<sunpy.net.attrs.Wavelength({self.min.value}, {self.max.value}, '{self.unit}')>"


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

        # . May be entered as a string or any numeric type for equality matching
        # . May be a string of the format '(min) - (max)' for range matching
        # . May be a string of the form '(operator) (number)' where operator is\
        one of: lt gt le ge < > <= >=

    """


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


class Physobs(SimpleAttr):
    """
    Specifies the physical observable the VSO can search for.

    Parameters
    ----------
    value : `str`
        A keyword describing the observable in the data.

    Notes
    -----
    More information about the values of physobs used by the VSO
    registry can be found at
    https://sdac.virtualsolar.org/cgi/show_details?keyword=PHYSOBS.
    """


class Provider(SimpleAttr):
    """
    Specifies the data provider to search for data using Fido.

    Parameters
    ----------
    value : str
        A keyword describing the Provider for the data.

    Notes
    -----
    For VSO, more information about each provider may be found within in the VSO Registry.
    See `VSO providers <https://sdac.virtualsolar.org/cgi/show_details?keyword=PROVIDER>`__.
    """


class Source(SimpleAttr):
    """
    Data sources that Fido can search with.

    Parameters
    ----------
    value : str
        A keyword describing the Data Source.

    Notes
    -----
    For VSO, more information about each source may be found within in the VSO Registry.
    See `VSO sources <https://sdac.virtualsolar.org/cgi/show_details?keyword=SOURCE>`__.
    Please note that 'Source' is used internally by VSO to represent
    what the VSO Data Model refers to as 'Observatory'.
    """


class ExtentType(SimpleAttr):
    """
    The type of Extent; for example, "FULLDISK", "SYNOPTIC", "LIMB", etc.
    """
