"""
TimeSeries is a generic time series class from which all other TimeSeries
classes inherit from.
"""

from __future__ import absolute_import, division, print_function
__authors__ = ["Alex Hamilton, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

import warnings
import inspect
from abc import ABCMeta
from collections import OrderedDict

import matplotlib.pyplot as plt
import pandas as pd

from sunpy import config
from sunpy.time import TimeRange, parse_time
from sunpy.extern import six
from sunpy.sun import sun

import astropy
import astropy.units as u
from astropy.table import Table
from astropy.table import Column

TIME_FORMAT = config.get("general", "time_format")


# pylint: disable=E1101,E1121,W0404,W0612,W0613
__authors__ = ["Alex Hamilton"]
__email__ = "####"


# GenericTimeSeries subclass registry.
TIMESERIES_CLASSES = OrderedDict()


class GenericTimeSeriesMeta(ABCMeta):
    """
    Registration metaclass for `~sunpy.timeseries.GenericTimeSeries`.
    This class checks for the existance of a method named ``is_datasource_for``
    when a subclass of `GenericTimeSeries` is defined. If it exists it will add
    that class to the registry.
    """

    _registry = TIMESERIES_CLASSES

    def __new__(mcls, name, bases, members):
        cls = super(GenericTimeSeriesMeta, mcls).__new__(mcls, name, bases, members)

        # The registry contains the class as the key and the validation method
        # as the item.
        if 'is_datasource_for' in members:
            mcls._registry[cls] = cls.is_datasource_for

        return cls


@six.add_metaclass(GenericTimeSeriesMeta)
class GenericTimeSeries:
    """
    A generic time series object.

    Attributes
    ----------
    meta : `str` or `dict`
        The comment string or header associated with the data.
    data : `~pandas.DataFrame`
        An pandas DataFrame prepresenting one or more fields as a function of time.

    Parameters
    ----------------
    filename: `str` or File
        A file to read.
    source: `str`
        A string identifier for a registered subclass, matched by that
         subclasses `_is_source_for` method.
    concatenate :  `bool`
        Concatenate all files into one TimeSeries object if True, or return
        one TimeSeries for each file if False.

    All other keywords are passed to _is_source_for and then __init__.

    Examples ########
    --------
    >>> from sunpy.timeseries import TimeSeries
    >>> import datetime
    >>> import numpy as np
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, step=(12 * np.pi) / (24 * 60)))
    >>> ts = Timeseries({"param1": intensity}, index=dates)
    >>> df = 
    >>> ts = Timeseries({"param1": intensity}, index=dates)
    >>> ts.peek()   # doctest: +SKIP

    References
    ----------
    * `Pandas Documentation <http://pandas.pydata.org/pandas-docs/dev/dsintro.html>`_

    """

    def __init__(self, data, meta=None, units=None, **kwargs):
        self.data = data
        if meta is None:
            self.meta = OrderedDict()
            self.meta.update({'name':None})
        else:
            self.meta = OrderedDict(meta)

        if units is None:
            self.units = {}
        else:
            self.units = units

        # Validate input data
        #self._validate_meta()
        #self._validate_units()

        # Setup some attributes
        self._nickname = self.detector


# #### Keyword attribute and other attribute definitions #### #
    @property
    def name(self):
        """Human-readable description of time-series-type"""
        return "{obs} {detector} {instrument} {measurement} {date:{tmf}}".format(obs=self.observatory,
                                                                detector=self.detector,
                                                                instrument=self.instrument,
                                                                measurement=self.measurement,
                                                                date=parse_time(self.date),
                                                                tmf=TIME_FORMAT)

    @property
    def nickname(self):
        """An abbreviated human-readable description of the time-series-type"""
        return self._nickname

    @nickname.setter
    def nickname(self, n):
        self._nickname = n

    @property
    def date(self):
        """Image observation time"""
        time = parse_time(self.meta.get('date-obs', self.data.index.min()))
        if time is None:
            warnings.warn_explicit("Missing metadata for observation time. Using first measurement index.",
                                       Warning, __file__, inspect.currentframe().f_back.f_lineno)
        return parse_time(time)

    @property
    def detector(self):
        """Detector name"""
        return self.meta.get('detector', "")

    @property
    def dsun(self):
        """The observer distance from the Sun."""
        dsun = self.meta.get('dsun_obs', None)

        if dsun is None:
            warnings.warn_explicit("Missing metadata for Sun-spacecraft separation: assuming Sun-Earth distance",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
            dsun = sun.sunearth_distance(self.date).to(u.m)

        return u.Quantity(dsun, 'm')

    @property
    def exposure_time(self):
        """Exposure time of the image in seconds."""
        return self.meta.get('exptime', 0.0) * u.s

    @property
    def instrument(self):
        """Instrument name"""
        return self.meta.get('instrume', "").replace("_", " ")

    @property
    def measurement(self):
        """Measurement name, defaults to the wavelength of image"""
        #return u.Quantity(self.meta.get('wavelnth', 0), self.meta.get('waveunit', ""))
        return 'measurement'

    @property
    def observatory(self):
        """Observatory or Telescope name"""
        return self.meta.get('obsrvtry', self.meta.get('telescop', "")).replace("_", " ")

# #### From Pandas #### #
    def sort_index(self, **kwargs):
        """Returns a sorted version of the TimeSeries object.
        Generally this shouldn't be necessary as most TimeSeries operations sort
        the data anyway to ensure consistent behaviour when truncating.

        Returns
        -------
        newts : `~sunpy.timeseries.TimeSeries`
            A new time series in ascending chronological order.
        """
        return GenericTimeSeries(self.data.sort_index(**kwargs), self.meta.copy(), self.units.copy())

# #### From Generic Spec #### #
    def resample(self, rule, method='sum', **kwargs):
        """Returns a resampled version of the TimeSeries object.

        Parameters
        ----------
        rule : `str`
            The offset string or object representing target conversion

        method : `str`
            The mothod used to combine values.
            Defines the function called in Pandas.
            downsampling: 'sum', 'mean', 'std'
            upsampling: 'pad', 'bfill', 'ffill'

        Returns
        -------
        newts :
            A new time series with the data resampled.
        """
        # how parameter indicates using multiple methods to Pandas
        how = kwargs.get('how', None)

        # Create the resample using the given rule.
        if not how:
            method = method.lower()
            resampled_data = self.data.resample(rule, **kwargs)
            if hasattr(resampled_data, method):
                resampled_data = getattr(resampled_data, method)()
            else:
                raise ValueError("Resample method not found")
        else:
            # If the how kwarg was given then we simply pass it to the resample function.
            resampled_data = self.data.resample(rule, **kwargs)
        
        # ToDo: consider re-evaluating the metadata.

        # Return the same type of timeseries object.
        return self.__class__(resampled_data.sort_index(), self.meta.copy(), self.units)


# #### From OLD LightCurve #### #
    @property
    def time_range(self):
        """Returns the start and end times of the TimeSeries as a `~sunpy.time.TimeRange`
        object"""
        return TimeRange(self.data.index.min(), self.data.index.max())

    def plot(self, axes=None, **plot_args):
        """Plot a plot of the time series

        Parameters
        ----------
        axes : `~matplotlib.axes.Axes` or None
            If provided the image will be plotted on the given axes. Otherwise
            the current axes will be used.

        **plot_args : `dict`
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        axes : `~matplotlib.axes.Axes`
            The plot axes.
        """

        # Get current axes
        if axes is None:
            axes = plt.gca()

        axes = self.data.plot(ax=axes, **plot_args)

        return axes

    def peek(self, **kwargs):
        """Displays the time series in a new figure.

        Parameters
        ----------
        **kwargs : `dict`
            Any additional plot arguments that should be used when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """

        figure = plt.figure()
        self.plot(**kwargs)
        figure.show()

        return figure

    def truncate(self, a, b=None, int=None):
        """Returns a truncated version of the TimeSeries object.

        Parameters
        ----------
        a : `sunpy.time.TimeRange`, `str` or `int`
            Either a time range to truncate to, or a start time in some format
            recognised by pandas, or a index integer.

        b : `str` or `int`
            If specified, the end time of the time range in some format
            recognised by pandas, or a index integer.

        int : `int`
            If specified, the interger indicating the slicing intervals.

        Returns
        -------
        newts : `~sunpy.timeseries.TimeSeries`
            A new time series with only the selected times.
        """
        # Debate: consider adding slice notation instead of calling this function.

        # Evaluate inputs
        # If given strings, then use to create a sunpy.time.timerange.TimeRange
        # for the SunPy text date parser.
        if isinstance(a, str) and isinstance(b, str):
            a = TimeRange(a, b)
        if isinstance(a, TimeRange):
            # If we have a TimeRange, extract the values
            start = a.start
            end   = a.end
        else:
            # Otherwise we already have the values
            start = a
            end   = b

        # If an interval integer was given then use in truncation.
        truncated = self.data.sort_index()[start:end:int]

        # Return the same type of timeseries object.
        return self.__class__(truncated.sort_index(), self.meta.copy(), self.units)

    def extract(self, column_name):
        """Returns a new time series with the chosen column.

        Parameters
        ----------
        column_name : `str`
            A valid column name.

        Returns
        -------
        newts : `~sunpy.timeseries.TimeSeries`
            A new time series with only the selected column.
        """
        """
        # TODO allow the extract function to pick more than one column
        if isinstance(self, pandas.Series):
            return self
        else:
            return GenericTimeSeries(self.data[column_name], self.meta.copy())
        """
        # Extract column and remove empty rows
        data = self.data[[column_name]].dropna()
        
        # Return the same type of timeseries object.
        return self.__class__(data.sort_index(), self.meta.copy(), self.units.copy())
        # Debate: should this return the same TimeSeries type as the original?

    def concatenate(self, otherts, **kwargs):
        """Concatenate with another time series. This function will check and
        remove any duplicate times. It will keep the column values from the
        original time series to which the new time series is being added.

        Parameters
        ----------
        otherts : `~sunpy.timeseries.TimeSeries`
            Another time series.

        same_source : `bool` Optional
            Set to true to check if the sources of the time series match.

        Returns
        -------
        newts : `~sunpy.timeseries.TimeSeries`
            A new time series.
        
        Debate: decide if we want to be able to concatenate multiple time series
        at once.
        """
        
        # Check the sources match if specified.
        same_source = kwargs.get('same_source', False)
        if same_source and not isinstance(otherts, self.__class__):
            raise TypeError("TimeSeries classes must match.")


        # Metadata is implemented as with the legacy LightCurve class.
        # ToDo: Debate: we want to consider how to do this in future.
        meta = OrderedDict()
        meta.update({str(self.data.index[0]):self.meta.copy()})
        meta.update({str(otherts.data.index[0]):otherts.meta.copy()})

        data = pd.concat([self.data.copy(), otherts.data])
        
        # Add all the new units to the dictionary.
        units = OrderedDict()
        units.update(self.units)
        units.update(otherts.units)
        
        # Return the same type of timeseries object.
        return self.__class__(data.sort_index(), meta, units)

# #### Miscellaneous #### #
    def _validate_meta(self):
        """
        Validates the meta-information associated with a TimeSeries.

        This method includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.timeseries.sources package.

        Allows for default unit assignment for:
            COL_UNITS

        """

        warnings.simplefilter('always', Warning)

        #
        for meta_property in ('cunit1', 'cunit2', 'waveunit'):
            if (self.meta.get(meta_property) and
                u.Unit(self.meta.get(meta_property),
                       parse_strict='silent').physical_type == 'unknown'):

                warnings.warn("Unknown value for "+meta_property.upper(), Warning)

    def _validate_units(self, units, **kwargs):
        """
        Validates the astropy unit-information associated with a TimeSeries.

        This method includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.timeseries.sources package.

        Allows for default unit assignment for:
            COL_UNITS

        """

        warnings.simplefilter('always', Warning)

        result = True
        for key in units:
            if not isinstance(units[key], astropy.units.UnitBase):
                # If this is not a unit then this can't be a valid units dict.
                result = False
                warnings.warn("Invalid unit given for \""+str(key)+"\"", Warning)
        
        return result

    def _sanitize_units(self, **kwargs):
        """
        Sanitises the collections.OrderedDict used to store the units.
        Primarily this method will:
        
        Remove entries that don't match up to a column,
        Add unitless entries for columns with no units defined.        
        Re-arrange the order of the dictionary to match the columns.
        """
        warnings.simplefilter('always', Warning)

        # Populate unspecified units:
        for column in set(self.data.columns.tolist()) - set(self.units.keys()):
            # For all columns not present in the units dictionary.
            self.units[column] = u.dimensionless_unscaled
            warnings.warn("Unknown units for \""+str(column)+"\"", Warning)
            
        # Re-arrange so it's in the same order as the columns and removed unused.
        units = OrderedDict()
        for column in self.data.columns.tolist():
            units.update({column:self.units[column]})
            
        # Now use the amended units Ordered Dictionary
        self.units = units
            

# #### New Methods #### #
    def to_table(self, **kwargs):
        """
        Return an Astropy Table of the give TimeSeries object.

        Returns
        -------
        newtable : `~astrpy.table`
            A new astropy table containing the data from the time series.
            The table will include units where relevant.
        """
        # ToDo: Table.from_pandas(df) doesn't include the index column. Add request?
        # Get data columns
        table = Table.from_pandas(self.data)

        # Get index column and add to table.
        index_col = Column(self.data.index.values, name='date')
        table.add_column(index_col, index=0)

        # Add in units.
        for key in self.units:
            table[key].unit = self.units[key]

        # Output the table
        return table

    def to_dataframe(self, **kwargs):
        """
        Return a Pandas DataFrame of the give TimeSeries object.

        Returns
        -------
        newdf : `~pandas.core.frame.DataFrame`
            A Pandas Dataframe containing the data.
        """
        return self.data

    def to_array(self, **kwargs):
        """
        Return a numpy array of the give TimeSeries object.

        Parameters
        ----------
        columns: `list`, optional, default:None
            If None, return all columns minus the index, otherwise, returns
            specified columns.

        Returns
        -------
        values : `~numpy.ndarray`
            If the caller is heterogeneous and contains booleans or objects,
            the result will be of dtype=object. See Notes.
        """
        return self.data.as_matrix(**kwargs)

    def quantity(self, colname, **kwargs):
        """
        Return a `~astropy.units.quantity.Quantity` for the given column.

        Parameters
        ----------
        colname: `str`
            The heading of the column you want output.

        Returns
        -------
        quantity : `~astropy.units.quantity.Quantity`
        """
        values = self.data[colname].values
        unit   = self.units[colname]
        return u.Quantity(values, unit)

    def add_column(self, colname, quantity, **kwargs):
        """
        Return an new TimeSeries with the given column added or updated.

        Parameters
        ----------
        colname: `str`
            The heading of the column you want output.

        quantity : `~astropy.units.quantity.Quantity` or `~numpy.ndarray`
            The values to be placed within the column.
            If updating values only then a numpy array is permitted.
            
        overwrite : `bool`, optional, default:True
            Set to true to allow the method to overwrite a column already present
            in the TimeSeries.
            
        Returns
        -------
        newts : TimeSeries
        
        """
        unit      = kwargs.get('unit', None)
        overwrite = kwargs.get('overwrite', True)
        values    = quantity

        # Make a copy of all the TimeSeries components.
        data  = self.data.copy()
        meta  = self.meta.copy()
        units = self.units.copy()
        
        # Check if we are updating a current column or adding a new one.
        if (colname in self.data.columns) and unit and overwrite:
            # Updating, so change the unit if appicable.
            units[colname] = unit
        
        # Convert the given values into the column units.
        if isinstance(values, astropy.units.quantity.Quantity) and overwrite:
            values = values.to(units[colname]).value
        
        # Update or add the data.
        if not (colname in self.data.columns) or overwrite:
            data[colname] = values
            
        # Return a new TimeSeries with the given updated/added column.
        return self.__class__(data, meta, units)

    @classmethod
    def _parse_file(cls, filepath):
        """Parses a file - to be implemented in any subclass that may use files"""
        return NotImplemented
