"""
This module provies `sunpy.timeseries.GenericTimeSeries` which all other
`sunpy.timeseries.TimeSeries` classes inherit from.
"""
import copy
import html
import textwrap
from collections import OrderedDict
from collections.abc import Iterable

import pandas as pd
from matplotlib.backend_bases import FigureCanvasBase
from matplotlib.figure import Figure

import astropy
import astropy.units as u
from astropy.table import Column, Table

from sunpy import config
from sunpy.time import TimeRange
from sunpy.timeseries import TimeSeriesMetaData
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.exceptions import warn_user
from sunpy.util.metadata import MetaDict
from sunpy.util.util import _figure_to_base64
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
    data : `~pandas.DataFrame` or `numpy.array`
        A `pandas.DataFrame` or `numpy.array` representing one or more fields as a function of time.
    meta : `~sunpy.timeseries.metadata.TimeSeriesMetaData`, optional
        The metadata giving details about the time series data/instrument.
        Defaults to `None`.
    units : `dict`, optional
        A mapping from column names in ``data`` to the physical units of that column.
        Defaults to `None`.

    Attributes
    ----------
    meta : `~sunpy.timeseries.metadata.TimeSeriesMetaData`
        The metadata giving details about the time series data/instrument.
    units : `dict`
        A mapping from column names in ``data`` to the physical units of that column.

    Examples
    --------
    >>> from sunpy.timeseries import TimeSeries
    >>> from sunpy.time import parse_time
    >>> from astropy.time import TimeDelta
    >>> import astropy.units as u
    >>> import numpy as np
    >>> import pandas as pd
    >>> times = parse_time("now") - TimeDelta(np.arange(24 * 60)*u.minute)
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, step=(12 * np.pi) / (24 * 60)))
    >>> df = pd.DataFrame(intensity, index=times, columns=['intensity'])
    >>> header = {}
    >>> units = {'intensity': u.W/u.m**2}
    >>> ts = TimeSeries(df, header, units)
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
        self._data = data
        tr = self.time_range
        # Check metadata input
        if meta is None:
            # No meta given, so default
            self.meta = TimeSeriesMetaData(MetaDict(), tr, self.columns)
        elif isinstance(meta, (dict, OrderedDict, MetaDict)):
            # Given the values for metadata (dict) and infer timerange and colnames from the data
            self.meta = TimeSeriesMetaData(meta, tr, self.columns)
        elif isinstance(meta, tuple):
            # Given the values all in a tuple
            self.meta = TimeSeriesMetaData(meta, tr, self.columns)
        else:
            # Should have a list of 3-tuples giving a complex metadata list.
            self.meta = meta

        if units is None:
            self.units = {}
        else:
            self.units = units

        for col in self.columns:
            if col not in self.units:
                warn_user(f'Unknown units for {col}')
                self.units[col] = u.dimensionless_unscaled

        # TODO: Fix this?
        # Validate input data
        # self._validate_meta()
        # self._validate_units()

    def _text_summary(self):
        start_time, end_time = self.time_range.start, self.time_range.end
        center_time = self.time_range.center
        duration = str(self.time_range.days.value) + " days"
        columns = ", ".join(self.columns)

        return textwrap.dedent(f"""\
                   Start Time:\t\t {start_time}
                   End Time:\t\t {end_time}
                   Center Time:\t\t {center_time}
                   Duration:\t\t {duration}
                   Columns:\t\t {columns}""")

    def __str__(self):
        return f"{self.meta}\n{self._data.__repr__()}"

    def __repr__(self):
        return f"{object.__repr__(self)}\n{self}"

    def _repr_html_(self):
        """
        Produce an HTML summary with plots for use in Jupyter notebooks.
        """

        # Convert the text repr to an HTML table
        partial_html = self._text_summary().replace('\n', '</td></tr><tr><th>')\
            .replace(':\t', '</th><td>')
        text_to_table = textwrap.dedent(f"""\
            <table style='text-align:left'>
                <tr><th>{partial_html}</td></tr>
            </table>""").replace('\n', '')

        # Plot the image in pixel space
        fig = Figure()
        # Figure instances in matplotlib<3.1 do not create a canvas by default
        if fig.canvas is None:
            FigureCanvasBase(fig)
        ax = fig.subplots()
        self.plot(axes=ax)
        img_src = _figure_to_base64(fig)

        return textwrap.dedent(f"""\
            <pre>{html.escape(object.__repr__(self))}</pre>
            <table>
                <tr>
                    <td>{text_to_table}</td>
                    <td>
                        <div align=center>
                           <img src='data:image/png;base64,{img_src}'
                        />
                        </div>
                    </td>
                </tr>
            </table>""")

# #### Attribute definitions #### #

    @property
    def data(self):
        """
        A `pandas.DataFrame` representing one or more fields as a function of time.
        """
        warn_user("Using .data to access the dataframe is discouraged; "
                  "use .to_dataframe() instead.")
        return self._data

    @data.setter
    def data(self, d):
        self._data = d

    @property
    def source(self):
        """
        A string/object used to specify the source class of the TimeSeries.
        """
        return self._source

    @property
    def observatory(self):
        """
        A string/object used to specify the observatory for the TimeSeries.
        """
        return

    @property
    def columns(self):
        """
        A list of all the names of the columns in the data.
        """
        return list(self._data.columns.values)

    @property
    def index(self):
        """
        The time index of the data.
        """
        return self._data.index

    @property
    def shape(self):
        """
        The shape of the data, a tuple (nrows, ncols).
        """
        return self._data.shape

    @property
    def time_range(self):
        """
        The start and end times of the TimeSeries as a `~sunpy.time.TimeRange`.
        """
        if len(self._data) > 0:
            return TimeRange(self._data.index.min(), self._data.index.max())
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
        values = self._data[colname].values
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
        data = copy.copy(self._data)
        meta = TimeSeriesMetaData(copy.copy(self.meta.metadata))
        units = copy.copy(self.units)

        # Add the unit to the units dictionary if already there.
        if not (colname in self._data.columns):
            units[colname] = unit

        # Convert the given quantity into values for given units if necessary.
        values = quantity
        if isinstance(values, astropy.units.quantity.Quantity) and overwrite:
            values = values.to(units[colname]).value

        # Update or add the data.
        if not (colname in self._data.columns) or overwrite:
            data[colname] = values

        # Return a new TimeSeries with the given updated/added column.
        return self.__class__(data, meta, units)

    def remove_column(self, colname):
        """
        Remove a column.

        Parameters
        ----------
        colname : str
            The heading of the column to remove.

        Returns
        -------
        `sunpy.timeseries.TimeSeries`
            A new `~sunpy.timeseries.TimeSeries`.
        """
        if colname not in self.columns:
            raise ValueError(f'Given column name ({colname}) not in list of columns {self.columns}')
        data = self._data.drop(colname, axis='columns')
        units = self.units.copy()
        units.pop(colname)
        return self.__class__(data, self.meta, units)

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
        return GenericTimeSeries(self._data.sort_index(**kwargs),
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
        truncated_data = self._data.sort_index()[start:end:int]

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
        #    return GenericTimeSeries(self._data[column_name], TimeSeriesMetaData(self.meta.metadata.copy()))

        # Extract column and remove empty rows
        data = self._data[[column_name]].dropna()
        units = {column_name: self.units[column_name]}

        # Build generic TimeSeries object and sanatise metadata and units.
        object = GenericTimeSeries(data.sort_index(),
                                   TimeSeriesMetaData(copy.copy(self.meta.metadata)),
                                   units)
        object._sanitize_metadata()
        return object

    def concatenate(self, others, same_source=False, **kwargs):
        """
        Concatenate with another `~sunpy.timeseries.TimeSeries` or an iterable containing multiple
        `~sunpy.timeseries.TimeSeries`. This function will check and remove any duplicate times.
        It will keep the column values from the original timeseries to which the new time
        series is being added.

        Parameters
        ----------
        others : `~sunpy.timeseries.TimeSeries` or `collections.abc.Iterable`
            Another `~sunpy.timeseries.TimeSeries` or an iterable containing multiple
            `~sunpy.timeseries.TimeSeries`.
        same_source : `bool`, optional
            Set to `True` to check if the sources of the time series match. Defaults to `False`.

        Returns
        -------
        `~sunpy.timeseries.TimeSeries`
            A new `~sunpy.timeseries.TimeSeries`.

        Notes
        -----
        Extra keywords are passed to `pandas.concat`.

        Examples
        --------
        A single `~sunpy.timeseries.TimeSeries` or an `collections.abc.Iterable` containing multiple
        `~sunpy.timeseries.TimeSeries` can be passed to concatenate.

        >>> timeseries_1.concatenate(timeseries_2) # doctest: +SKIP
        >>> timeseries_1.concatenate([timeseries_2, timeseries_3]) # doctest: +SKIP

        Set ``same_source`` to `True` if the sources of the time series are the same.

        >>> timeseries_1.concatenate([timeseries_2, timeseries_3], same_source=True) # doctest: +SKIP
        """
        # Check to see if nothing needs to be done in case the same TimeSeries is provided.
        if self == others:
            return self
        elif isinstance(others, Iterable):
            if len(others) == 1 and self == next(iter(others)):
                return self

        # Check the sources match if specified.
        if (
            same_source
            and isinstance(others, Iterable)
            and not all(isinstance(series, self.__class__) for series in others)
        ):
            raise TypeError("TimeSeries classes must match if 'same_source' is specified.")
        elif (
            same_source
            and not isinstance(others, Iterable)
            and not isinstance(others, self.__class__)
        ):
            raise TypeError("TimeSeries classes must match if 'same_source' is specified.")

        # If an iterable is not provided, it must be a TimeSeries object, so wrap it in a list.
        if not isinstance(others, Iterable):
            others = [others]

        # Concatenate the metadata and data.
        kwargs["sort"] = kwargs.pop("sort", False)
        meta = self.meta.concatenate([series.meta for series in others])
        data = pd.concat(
            [self._data.copy(), *list(series.to_dataframe() for series in others)], **kwargs
        )

        # Add all the new units to the dictionary.
        units = OrderedDict()
        units.update(self.units)
        units.update(
            {k: v for unit in list(series.units for series in others) for k, v in unit.items()}
        )
        units = {k: v for k, v in units.items() if k in data.columns}

        # If sources match then build similar TimeSeries.
        if all(self.__class__ == series.__class__ for series in others):
            object = self.__class__(data.sort_index(), meta, units)
        else:
            # Build generic time series if the sources don't match.
            object = GenericTimeSeries(data.sort_index(), meta, units)

        # Sanatise metadata and units
        object._sanitize_metadata()
        return object

# #### Plotting Methods #### #

    def plot(self, axes=None, columns=None, **plot_args):
        """
        Plot a plot of the `~sunpy.timeseries.TimeSeries`.

        Parameters
        ----------
        axes : `~matplotlib.axes.Axes`, optional
            If provided the image will be plotted on the given axes.
            Defaults to `None`, so the current axes will be used.
        columns : list[str], optional
            If provided, only plot the specified columns.
        **plot_args : `dict`, optional
            Additional plot keyword arguments that are handed to
            :meth:`pandas.DataFrame.plot`.

        Returns
        -------
        `~matplotlib.axes.Axes`
            The plot axes.
        """
        import matplotlib.pyplot as plt

        # Get current axes
        if axes is None:
            axes = plt.gca()

        if columns is None:
            columns = self._data.columns

        axes = self._data[columns].plot(ax=axes, **plot_args)

        units = set([self.units[col] for col in columns])
        if len(units) == 1:
            # If units of all columns being plotted are the same, add a unit
            # label to the y-axis.
            unit = u.Unit(list(units)[0])
            axes.set_ylabel(unit.to_string())
        return axes

    @peek_show
    def peek(self, columns=None, **kwargs):
        """
        Displays a graphical overview of the data in this object for user evaluation.
        For the creation of plots, users should instead use the
        `~sunpy.timeseries.GenericTimeSeries.plot` method and Matplotlib's pyplot framework.

        Parameters
        ----------
        columns : list[str], optional
            If provided, only plot the specified columns.
        **kwargs : `dict`
            Any additional plot arguments that should be used when plotting.
        """
        import matplotlib.pyplot as plt

        # Check we have a timeseries valid for plotting
        self._validate_data_for_plotting()

        # Now make the plot
        figure = plt.figure()
        self.plot(columns=columns, **kwargs)

        return figure

    def _validate_data_for_plotting(self):
        """
        Raises an exception if the `~sunpy.timeseries.TimeSeries` is invalid
        for plotting.

        This should be added into all `~sunpy.timeseries.TimeSeries`
        peek methods.
        """
        # Check we have a valid TS
        if len(self._data) == 0:
            raise ValueError("The timeseries can't be plotted as it has no data present. "
                             "(len(self._data) == 0)")

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
        for meta_property in ('cunit1', 'cunit2', 'waveunit'):
            if (self.meta.get(meta_property) and
                u.Unit(self.meta.get(meta_property),
                       parse_strict='silent').physical_type == 'unknown'):

                warn_user(f"Unknown value for {meta_property.upper()}.")

    def _validate_units(self, units, **kwargs):
        """
        Validates the astropy unit-information associated with a
        `~sunpy.timeseries.TimeSeries`.

        This method includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-
        specific validation should be handled in the relevant file in
        the "sunpy.timeseries.sources".
        """
        result = True
        for key in units:
            if not isinstance(units[key], astropy.units.UnitBase):
                # If this is not a unit then this can't be a valid units dict.
                result = False
                warn_user(f"Invalid unit given for {key}.")

        return result

    def _sanitize_metadata(self, **kwargs):
        """
        Sanitizes the `~sunpy.timeseries.TimeSeriesMetaData`  used to store the
        metadata.

        Primarily this method will:

        * Remove entries outside of the dates or truncate if the metadata overflows past the data.
        * Remove column references in the metadata that don't match to a column in the data.
        * Remove metadata entries that have no columns matching the data.
        """
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
        table = Table.from_pandas(self._data)

        # Get index column and add to table.
        index_col = Column(self._data.index.values, name='date')
        table.add_column(index_col, index=0)

        # Add in units.
        for key in self.units:
            table[key].unit = self.units[key]

        # Output the table
        return table

    def to_dataframe(self, **kwargs):
        """
        Return a `~pandas.DataFrame` of the given
        `~sunpy.timeseries.TimeSeries`.

        Returns
        -------
        `~pandas.DataFrame`
        """
        return self._data

    def to_array(self, **kwargs):
        """
        Return a `numpy.array` of the given `~sunpy.timeseries.TimeSeries`.

        Parameters
        ----------
        **kwargs : `dict`
            All keyword arguments are passed to `pandas.DataFrame.to_numpy`.

        Returns
        -------
        `~numpy.ndarray`
            If the data is heterogeneous and contains booleans or objects, the result will be of ``dtype=object``.
        """
        if hasattr(self._data, "to_numpy"):
            return self._data.to_numpy(**kwargs)
        else:
            return self._data.values

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
            if ((not self._data.equals(other.to_dataframe())) or
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
        raise NoMatchError(f'Could not find any timeseries sources to parse {filepath}')
