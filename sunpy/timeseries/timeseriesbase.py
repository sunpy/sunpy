"""
TimeSeries is a generic time series class from which all other TimeSeries
classes inherit from.
"""

from __future__ import absolute_import, division, print_function
__authors__ = ["Alex Hamilton, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

import warnings
from abc import ABCMeta
from collections import OrderedDict
import copy

import matplotlib.pyplot as plt
import pandas as pd

from sunpy import config
from sunpy.time import TimeRange
from sunpy.extern import six
from sunpy.timeseries import TimeSeriesMetaData
from sunpy.util.metadata import MetaDict

import astropy
import astropy.units as u
from astropy.table import Table
from astropy.table import Column

# define and register a new unit, needed for RHESSI
det = u.def_unit('detector')
u.add_enabled_units([det])

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

    Parameters
    ----------
    data : `~pandas.DataFrame`
        A pandas DataFrame representing one or more fields as a function
        of time.
    meta : `~sunpy.timeseries.metadata.TimeSeriesMetaData`, optional
        The metadata giving details about the time series data/instrument.
    units : dict, optional
        A mapping from column names in *data* to the physical units of
        that column.

    Attributes
    ----------
    data : `~pandas.DataFrame`
        A pandas DataFrame representing one or more fields as a function
        of time.
    meta : `~sunpy.timeseries.metadata.TimeSeriesMetaData`
        The metadata giving details about the time series data/instrument.
    units : dict
        A mapping from column names in *data* to the physical units of
        that column.

    Examples
    --------
    >>> from sunpy.timeseries import TimeSeries
    >>> import datetime
    >>> import numpy as np
    >>> import pandas as pd
    >>> base = datetime.datetime.today()
    >>> times = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, step=(12 * np.pi) / (24 * 60)))
    >>> df = pd.DataFrame(intensity, index=times, columns=['intensity'])
    >>> ts = TimeSeries(df)
    >>> ts.peek()   # doctest: +SKIP

    References
    ----------
    * `Pandas Documentation <https://pandas.pydata.org/pandas-docs/stable/>`_

    """

    # Class attribute used to specify the source class of the TimeSeries.
    _source = None

    def __init__(self, data, meta=None, units=None, **kwargs):
        self.data = data
        tr = TimeRange(self.data.index.min(), self.data.index.max())
        # Check metadata input
        if meta is None:
            # No meta given, so default
            self.meta = TimeSeriesMetaData(MetaDict(), tr, list(self.data.columns.values))
        elif isinstance(meta, (dict, OrderedDict, MetaDict)):
            # Given the values for metadata (dict) and infer timerange and colnames from the data
            self.meta = TimeSeriesMetaData(meta, tr, list(self.data.columns.values))
        elif isinstance(meta, tuple):
            # Given the values all in a tuple
            self.meta = TimeSeriesMetaData(meta, tr, list(self.data.columns.values))
        else:
            # Should have a list of 3-tuples giving a complex metadata list.
            self.meta = meta

        if units is None:
            self.units = {}
        else:
            self.units = units

        # Validate input data
        #self._validate_meta()
        #self._validate_units()

# #### Attribute definitions #### #

    @property
    def source(self):
        """
        A string/object used to specify the source class of the TimeSeries.
        """
        return self._source

    @property
    def columns(self):
        """A list of all the names of the columns in the data."""
        return list(self.data.columns.values)

    @property
    def index(self):
        """The time index of the data."""
        return self.data.index

    @property
    def time_range(self):
        """
        The start and end times of the TimeSeries as a `~sunpy.time.TimeRange`
        object
        """
        return TimeRange(self.data.index.min(), self.data.index.max())

# #### Data Access, Selection and Organisation Methods #### #

    def quantity(self, colname, **kwargs):
        """
        Return a `~astropy.units.quantity.Quantity` for the given column.

        Parameters
        ----------
        colname : `str`
            The heading of the column you want output.

        Returns
        -------
        quantity : `~astropy.units.quantity.Quantity`
        """
        values = self.data[colname].values
        unit   = self.units[colname]
        return u.Quantity(values, unit)

    def add_column(self, colname, quantity, unit=False, overwrite=True, **kwargs):
        """
        Return an new TimeSeries with the given column added or updated.

        Parameters
        ----------
        colname : `str`
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
        # Get the expected units from the quantity if required
        if not unit and isinstance(quantity, astropy.units.quantity.Quantity):
            unit = quantity.unit
        elif not unit:
            unit = u.dimensionless_unscaled

        # Make a copy of all the TimeSeries components.
        data  = copy.copy(self.data)
        meta  = TimeSeriesMetaData(copy.copy(self.meta.metadata))
        units = copy.copy(self.units)

        # Add the unit to the units dictionary if already there.
        if not (colname in self.data.columns):
            units[colname] = unit

        # Convert the given quantity into values for given units if necessary.
        values = quantity
        if isinstance(values, astropy.units.quantity.Quantity) and overwrite:
            values = values.to(units[colname]).value

        # Update or add the data.
        if not (colname in self.data.columns) or overwrite:
            data[colname] = values

        # Return a new TimeSeries with the given updated/added column.
        return self.__class__(data, meta, units)

    def sort_index(self, **kwargs):
        """Returns a sorted version of the TimeSeries object.
        Generally this shouldn't be necessary as most TimeSeries operations sort
        the data anyway to ensure consistent behaviour when truncating.

        Returns
        -------
        newts : `~sunpy.timeseries.TimeSeries`
            A new time series in ascending chronological order.
        """
        return GenericTimeSeries(self.data.sort_index(**kwargs), TimeSeriesMetaData(copy.copy(self.meta.metadata)), copy.copy(self.units))

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
            If specified, the integer indicating the slicing intervals.

        Returns
        -------
        newts : `~sunpy.timeseries.TimeSeries`
            A new time series with only the selected times.
        """
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
        truncated_data = self.data.sort_index()[start:end:int]

        # Truncate the metadata
        # Check there is data still
        truncated_meta = TimeSeriesMetaData([])
        if len(truncated_data) > 0:
            tr = TimeRange(truncated_data.index.min(), truncated_data.index.max())
            truncated_meta = TimeSeriesMetaData(copy.deepcopy(self.meta.metadata))
            truncated_meta._truncate(tr)

        # Build similar TimeSeries object and sanatise metadata and units.
        object = self.__class__(truncated_data.sort_index(), truncated_meta, copy.copy(self.units))
        object._sanitize_metadata()
        object._sanitize_units()
        return object

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
            return GenericTimeSeries(self.data[column_name], TimeSeriesMetaData(self.meta.metadata.copy()))
        """
        # Extract column and remove empty rows
        data = self.data[[column_name]].dropna()

        # Build generic TimeSeries object and sanatise metadata and units.
        object = GenericTimeSeries(data.sort_index(), TimeSeriesMetaData(copy.copy(self.meta.metadata)), copy.copy(self.units))
        object._sanitize_metadata()
        object._sanitize_units()
        return object

    def concatenate(self, otherts, **kwargs):
        """Concatenate with another TimeSeries. This function will check and
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

        # check to see if nothing needs to be done
        if self == otherts:
            return self

        # Check the sources match if specified.
        same_source = kwargs.get('same_source', False)
        if same_source and not (isinstance(otherts, self.__class__)):
            raise TypeError("TimeSeries classes must match if specified.")

        # Concatenate the metadata and data
        meta = self.meta.concatenate(otherts.meta)
        data = pd.concat([self.data.copy(), otherts.data], **kwargs)

        # Add all the new units to the dictionary.
        units = OrderedDict()
        units.update(self.units)
        units.update(otherts.units)

        # If sources match then build similar TimeSeries.
        if self.__class__ == otherts.__class__:
            object = self.__class__(data.sort_index(), meta, units)
        else:
            # Build generic time series if the sources don't match.
            object = GenericTimeSeries(data.sort_index(), meta, units)

        # Sanatise metadata and units
        object._sanitize_metadata()
        object._sanitize_units()
        return object

# #### Plotting Methods #### #

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
        """
        # Check we have a timeseries valid for plotting
        self._validate_data_for_ploting()

        # Now make the plot
        figure = plt.figure()
        self.plot(**kwargs)
        figure.show()

    def _validate_data_for_ploting(self):
        """Raises an exception if the timeseries is invalid for plotting.
        To be added into all the peek methods in all source sup-classes.
        Currently only checks if we have an empty timeseries, where:
        len(self.data) == 0

        """
        # Check we have a valid TS
        if len(self.data) == 0:
            raise ValueError('The timeseries can\'t be plotted as it has no data present. (len(self.data) == 0)')

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

    def _sanitize_metadata(self, **kwargs):
        """
        Sanitises the TimeSeriesMetaData object used to store the metadata.
        Primarily this method will:

        Remove entries outside of the datas TimeRange or truncate TimeRanges
        if the metadata overflows past the data,
        Remove column references in the metadata that don't match to a column
        in the data.
        Remove metadata entries that have no columns matching the data.
        """
        warnings.simplefilter('always', Warning)

        # Truncate the metadata
        self.meta._truncate(self.time_range)

        # Remove non-existant columns
        redundant_cols = list(set(self.meta.columns) - set(self.columns))
        self.meta._remove_columns(redundant_cols)

# #### Export/Output Methods #### #

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

    def __eq__(self, other):
        """
        Check two TimeSeries objects are the same, they have matching type, data,
        metadata and units entries.

        Parameters
        ----------
        other : `~sunpy.timeseries.GenericTimeSeries`
            The second TimeSeries object to compare with.

        Returns
        -------
        result : `bool`
        """
        match = True
        if isinstance(other, type(self)):
            if (not self.data.equals(other.data)) or (self.meta != other.meta) or (self.units != other.units):
                match = False
        else:
            match = False
        return match

    def __ne__(self, other):
        """
        Check two TimeSeries objects are not the same, they don't have matching
        type, data, metadata and/or units entries.

        Parameters
        ----------
        other : `~sunpy.timeseries.GenericTimeSeries`
            The second TimeSeries object to compare with.

        Returns
        -------
        result : `bool`
        """
        return not self == other

    @classmethod
    def _parse_file(cls, filepath):
        """Parses a file - to be implemented in any subclass that may use files"""
        return NotImplemented
