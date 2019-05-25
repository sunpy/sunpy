"""
This module provies `sunpy.timeseries.GenericTimeSeries` which all other
`sunpy.timeseries.TimeSeries` classes inherit from.
"""
import copy
import warnings
from collections import OrderedDict

import matplotlib.pyplot as plt
import pandas as pd

import astropy
import astropy.units as u
from astropy.table import Column, Table

from sunpy import config
from sunpy.time import TimeRange
from sunpy.timeseries import TimeSeriesMetaData
from sunpy.util.exceptions import SunpyUserWarning
from sunpy.util.metadata import MetaDict
from sunpy.visualization import peek_show

# define and register a new unit, needed for RHESSI
det = u.def_unit('detector')
u.add_enabled_units([det])

TIME_FORMAT = config.get("general", "time_format")

__all__ = ["GenericTimeSeries"]


class GenericTimeSeries:
    """
    A generic time series object.

    Parameters
    ----------
    data : `~pandas.DataFrame`
        A `pandas.DataFrame` representing one or more fields as a function of time.
    meta : `~sunpy.timeseries.metadata.TimeSeriesMetaData`, optional
        The metadata giving details about the time series data/instrument.
        Defaults to `None`.
    units : `dict`, optional
        A mapping from column names in ``data`` to the physical units of that column.
        Defaults to `None`.

    Attributes
    ----------
    data : `~pandas.DataFrame`
        A `pandas.DataFrame` representing one or more fields as a function of time.
    meta : `~sunpy.timeseries.metadata.TimeSeriesMetaData`
        The metadata giving details about the time series data/instrument.
    units : `dict`
        A mapping from column names in ``data`` to the physical units ofthat column.

    Examples
    --------
    >>> from sunpy.timeseries import TimeSeries
    >>> from sunpy.time import parse_time
    >>> from astropy.time import TimeDelta
    >>> import numpy as np
    >>> import pandas as pd
    >>> times = parse_time("now") - TimeDelta(np.arange(24 * 60)*u.minute)
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, step=(12 * np.pi) / (24 * 60)))
    >>> df = pd.DataFrame(intensity, index=times, columns=['intensity'])
    >>> ts = TimeSeries(df)
    >>> ts.peek()  # doctest: +SKIP

    References
    ----------
    * `Pandas Documentation <https://pandas.pydata.org/pandas-docs/stable/>`_
    """
    # Class attribute used to specify the source class of the TimeSeries.
    _source = None
    _registry = dict()

    def __init_subclass__(cls, **kwargs):
        """
        An __init_subclass__ hook initializes all of the subclasses of a given
        class.

        So for each subclass, it will call this block of code on import.
        This replicates some metaclass magic without the need to be
        aware of metaclasses. Here we use this to register each subclass
        in a dict that has the `is_datasource_for` attribute. This is
        then passed into the TimeSeries Factory so we can register them.
        """
        super().__init_subclass__(**kwargs)
        if hasattr(cls, 'is_datasource_for'):
            cls._registry[cls] = cls.is_datasource_for

    # kwargs are not used here but are passed in for sources.
    def __init__(self, data, meta=None, units=None, **kwargs):
        self.data = data
        tr = self.time_range
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

        # TODO: Fix this?
        # Validate input data
        # self._validate_meta()
        # self._validate_units()

# #### Attribute definitions #### #

    @property
    def source(self):
        """
        A string/object used to specify the source class of the TimeSeries.
        """
        return self._source

    @property
    def columns(self):
        """
        A list of all the names of the columns in the data.
        """
        return list(self.data.columns.values)

    @property
    def index(self):
        """
        The time index of the data.
        """
        return self.data.index

    @property
    def time_range(self):
        """
        The start and end times of the TimeSeries as a `~sunpy.time.TimeRange`.
        """
        if len(self.data) > 0:
            return TimeRange(self.data.index.min(), self.data.index.max())
        else:
            return None

# #### Data Access, Selection and Organisation Methods #### #

    def quantity(self, colname, **kwargs):
        """
        Return a `~astropy.units.quantity.Quantity` for the given column.

        Parameters
        ----------
        colname : `str`
            The heading of the column you want to output.

        Returns
        -------
        `~astropy.units.quantity.Quantity`
        """
        values = self.data[colname].values
        unit = self.units[colname]
        return u.Quantity(values, unit)

    def add_column(self, colname, quantity, unit=False, overwrite=True, **kwargs):
        """
        Return a new `~sunpy.timeseries.TimeSeries` with the given column added
        or updated.

        Parameters
        ----------
        colname : `str`
            The heading of the column you want output.
        quantity : `~astropy.units.quantity.Quantity` or `~numpy.ndarray`
            The values to be placed within the column.
            If updating values only then a numpy array is permitted.
        overwrite : `bool`, optional
            Defaults to `True`, allowing the method to overwrite a column already present in the `~sunpy.timeseries.TimeSeries`.

        Returns
        -------
        `sunpy.timeseries.TimeSeries`
            A new `~sunpy.timeseries.TimeSeries`.
        """
        # Get the expected units from the quantity if required
        if not unit and isinstance(quantity, astropy.units.quantity.Quantity):
            unit = quantity.unit
        elif not unit:
            unit = u.dimensionless_unscaled

        # Make a copy of all the TimeSeries components.
        data = copy.copy(self.data)
        meta = TimeSeriesMetaData(copy.copy(self.meta.metadata))
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
        """
        Returns a sorted version of a `~sunpy.timeseries.TimeSeries`. Generally
        this shouldn't be necessary as most `~sunpy.timeseries.TimeSeries`
        operations sort the data anyway to ensure consistent behavior when
        truncating.

        Returns
        -------
        `~sunpy.timeseries.TimeSeries`
            A new `~sunpy.timeseries.TimeSeries` in ascending chronological order.
        """
        return GenericTimeSeries(self.data.sort_index(**kwargs),
                                 TimeSeriesMetaData(copy.copy(self.meta.metadata)),
                                 copy.copy(self.units))

    def truncate(self, a, b=None, int=None):
        """
        Returns a truncated version of the TimeSeries object.

        Parameters
        ----------
        a : `sunpy.time.TimeRange`, `str`, `int`
            Either a time range to truncate to, or a start time in some format recognized by pandas, or a index integer.
        b : `str` or `int`, optional
            If specified, the end time of the time range in some format recognized by pandas, or a index integer.
            Defaults to `None`.
        int : `int`, optional
            If specified, the integer indicating the slicing intervals.
            Defaults to `None`.

        Returns
        -------
        `~sunpy.timeseries.TimeSeries`
            A new `~sunpy.timeseries.TimeSeries` with only the selected times.
        """
        # Evaluate inputs
        # If given strings, then use to create a sunpy.time.timerange.TimeRange
        # for the SunPy text date parser.
        if isinstance(a, str) and isinstance(b, str):
            a = TimeRange(a, b)
        if isinstance(a, TimeRange):
            # If we have a TimeRange, extract the values
            start = a.start.datetime
            end = a.end.datetime
        else:
            # Otherwise we already have the values
            start = a
            end = b

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
        """
        Returns a new time series with the chosen column.

        Parameters
        ----------
        column_name : `str`
            A valid column name.

        Returns
        -------
        `~sunpy.timeseries.TimeSeries`
            A new `~sunpy.timeseries.TimeSeries` with only the selected column.
        """
        # TODO: allow the extract function to pick more than one column
        # TODO: Fix this?
        # if isinstance(self, pandas.Series):
        #    return self
        # else:
        #    return GenericTimeSeries(self.data[column_name], TimeSeriesMetaData(self.meta.metadata.copy()))

        # Extract column and remove empty rows
        data = self.data[[column_name]].dropna()

        # Build generic TimeSeries object and sanatise metadata and units.
        object = GenericTimeSeries(data.sort_index(),
                                   TimeSeriesMetaData(copy.copy(self.meta.metadata)),
                                   copy.copy(self.units))
        object._sanitize_metadata()
        object._sanitize_units()
        return object

    def concatenate(self, otherts, same_source=False, **kwargs):
        """
        Concatenate with another `~sunpy.timeseries.TimeSeries`. This function
        will check and remove any duplicate times. It will keep the column
        values from the original timeseries to which the new time series is
        being added.

        Parameters
        ----------
        otherts : `~sunpy.timeseries.TimeSeries`
            Another `~sunpy.timeseries.TimeSeries`.
        same_source : `bool`, optional
            Set to `True` to check if the sources of the time series match. Defaults to `False`.

        Returns
        -------
        `~sunpy.timeseries.TimeSeries`
            A new `~sunpy.timeseries.TimeSeries`.

        Notes
        -----
        Extra keywords are passed to `pandas.concat`.
        """
        # TODO: decide if we want to be able to concatenate multiple time series at once.
        # check to see if nothing needs to be done
        if self == otherts:
            return self

        # Check the sources match if specified.
        if same_source and not (isinstance(otherts, self.__class__)):
            raise TypeError("TimeSeries classes must match if specified.")

        # Concatenate the metadata and data
        kwargs['sort'] = kwargs.pop('sort', False)
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
        """
        Plot a plot of the `~sunpy.timeseries.TimeSeries`.

        Parameters
        ----------
        axes : `~matplotlib.axes.Axes`, optional
            If provided the image will be plotted on the given axes.
            Defaults to `None`, so the current axes will be used.
        **plot_args : `dict`, optional
            Any additional plot arguments that should be used when plotting.

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

    @peek_show
    def peek(self, **kwargs):
        """
        Displays a graphical overview of the data in this object for user evaluation.
        For the creation of plots, users should instead use the
        `~sunpy.timeseries.GenericTimeSeries.plot` method and Matplotlib's pyplot framework.

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

        return figure

    def _validate_data_for_ploting(self):
        """
        Raises an exception if the `~sunpy.timeseries.TimeSeries` is invalid
        for plotting.

        This should be added into all `~sunpy.timeseries.TimeSeries`
        peek methods.
        """
        # Check we have a valid TS
        if len(self.data) == 0:
            raise ValueError("The timeseries can't be plotted as it has no data present. "
                             "(len(self.data) == 0)")

# #### Miscellaneous #### #

    def _validate_meta(self):
        """
        Validates the meta-information associated with a
        `~sunpy.timeseries.TimeSeries`.

        This method includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-
        specific validation should be handled in the relevant file in
        the "sunpy.timeseries.sources".
        """
        warnings.simplefilter('always', Warning)

        for meta_property in ('cunit1', 'cunit2', 'waveunit'):
            if (self.meta.get(meta_property) and
                u.Unit(self.meta.get(meta_property),
                       parse_strict='silent').physical_type == 'unknown'):

                warnings.warn(f"Unknown value for {meta_property.upper()}.", SunpyUserWarning)

    def _validate_units(self, units, **kwargs):
        """
        Validates the astropy unit-information associated with a
        `~sunpy.timeseries.TimeSeries`.

        This method includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-
        specific validation should be handled in the relevant file in
        the "sunpy.timeseries.sources".
        """
        warnings.simplefilter('always', Warning)

        result = True
        for key in units:
            if not isinstance(units[key], astropy.units.UnitBase):
                # If this is not a unit then this can't be a valid units dict.
                result = False
                warnings.warn(f"Invalid unit given for {key}.", SunpyUserWarning)

        return result

    def _sanitize_units(self, **kwargs):
        """
        Sanitizes the `collections.OrderedDict` used to store the units.

        Primarily this method will:

        * Remove entries that don't match up to a column.
        * Add unitless entries for columns with no units defined.
        * Re-arrange the order of the dictionary to match the columns.
        """
        warnings.simplefilter('always', Warning)

        # Populate unspecified units:
        for column in set(self.data.columns.tolist()) - set(self.units.keys()):
            # For all columns not present in the units dictionary.
            self.units[column] = u.dimensionless_unscaled
            warnings.warn(f"Unknown units for {column}.", SunpyUserWarning)

        # Re-arrange so it's in the same order as the columns and removed unused.
        units = OrderedDict()
        for column in self.data.columns.tolist():
            units.update({column: self.units[column]})

        # Now use the amended units Ordered Dictionary
        self.units = units

    def _sanitize_metadata(self, **kwargs):
        """
        Sanitizes the `~sunpy.timeseries.TimeSeriesMetaData`  used to store the
        metadata.

        Primarily this method will:

        * Remove entries outside of the dates or truncate if the metadata overflows past the data.
        * Remove column references in the metadata that don't match to a column in the data.
        * Remove metadata entries that have no columns matching the data.
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
        Return an `astropy.table.Table` of the given
        `~sunpy.timeseries.TimeSeries`.

        Returns
        -------
        `~astropy.table.Table`
            A new `astropy.table.Table` containing the data from the `~sunpy.timeseries.TimeSeries`.
            The table will include units where relevant.
        """
        # TODO: Table.from_pandas(df) doesn't include the index column. Add request?
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
        Return a `~pandas.core.frame.DataFrame` of the given
        `~sunpy.timeseries.TimeSeries`.

        Returns
        -------
        `~pandas.core.frame.DataFrame`
            A `~pandas.core.frame.DataFrame` containing the data.
        """
        return self.data

    def to_array(self, **kwargs):
        """
        Return a `numpy.array` of the given `~sunpy.timeseries.TimeSeries`.

        Parameters
        ----------
        kwargs : `dict`
            All keyword arguments are passed to `pandas.DataFrame.to_numpy`.

        Returns
        -------
        `~numpy.ndarray`
            If the data is heterogeneous and contains booleans or objects, the result will be of ``dtype=object``.
        """
        if hasattr(self.data, "to_numpy"):
            return self.data.to_numpy(**kwargs)
        else:
            return self.data.values

    def __eq__(self, other):
        """
        Check two `~sunpy.timeseries.TimeSeries` are the same, they have
        matching type, data, metadata and units entries.

        Parameters
        ----------
        other : `~sunpy.timeseries.TimeSeries`
            The second `~sunpy.timeseries.TimeSeries` to compare with.

        Returns
        -------
        `bool`
        """
        match = True
        if isinstance(other, type(self)):
            if ((not self.data.equals(other.data)) or
                    (self.meta != other.meta) or
                    (self.units != other.units)):
                match = False
        else:
            match = False
        return match

    def __ne__(self, other):
        """
        Check two `~sunpy.timeseries.TimeSeries` are not the same, they don't
        have matching type, data, metadata and/or units entries.

        Parameters
        ----------
        other : `~sunpy.timeseries.TimeSeries`
            The second `~sunpy.timeseries.TimeSeries` to compare with.

        Returns
        -------
        `bool`
        """
        return not self == other

    @classmethod
    def _parse_file(cls, filepath):
        """
        Parses a file - to be implemented in any subclass that may use files.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        return NotImplemented
