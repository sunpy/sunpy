"""
This module provies the `~sunpy.timeseries.TimeSeriesFactory` class.
"""
import os
import copy
import glob
from collections import OrderedDict
from urllib.request import urlopen

import numpy as np
import pandas as pd

import astropy
import astropy.io.fits
import astropy.units as u
from astropy.table import Table
from astropy.time import Time

import sunpy
from sunpy.io.file_tools import UnrecognizedFileTypeError, read_file
from sunpy.io.fits import HDPair
from sunpy.io.header import FileHeader
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util import expand_list
from sunpy.util.config import get_and_create_download_dir
from sunpy.util.datatype_factory_base import (
    BasicRegistrationFactory,
    MultipleMatchError,
    NoMatchError,
    ValidationFunctionError,
)
from sunpy.util.metadata import MetaDict
from sunpy.util.net import download_file

__all__ = ['TimeSeries', 'TimeSeriesFactory', 'NoTimeSeriesFound',
           'InvalidTimeSeriesInput', 'InvalidTimeSeriesType']


class TimeSeriesFactory(BasicRegistrationFactory):
    """
    A factory for generating solar timeseries objects.

    This factory takes a variety of inputs to generate
    `~sunpy.timeseries.GenericTimeSeries` objects.

    Parameters
    ----------
    \\*inputs
        Inputs to parse for timeseries objects. See the example section for a
        detailed list of possible inputs.

    source : `str`, optional
        A string to select the observational source of the data, currently
        necessary to define how files should be read for all instruments.

    concatenate : `bool`, optional
        Defaults to `False`.
        If set, combine any resulting list of TimeSeries objects into a single
        TimeSeries, using successive concatenate methods.

    Returns
    -------
    `sunpy.timeseries.GenericTimeSeries`
        If the input results in a single timeseries object that will be returned, or if ``concatenate=True``.

    `list` of `~sunpy.timeseries.GenericTimeSeries`
        If multiple inputs are parsed, they will be returned in a list, unless
        ``concatenate=True`` is set when they will be combined into a single
        timeseries.

    Examples
    --------
    >>> import sunpy.timeseries
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> my_timeseries = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)  # doctest: +REMOTE_DATA

    The SunPy TimeSeries factory accepts a wide variety of inputs for creating timeseries

    * Preloaded tuples of (data, header) pairs or (data, header, units)

    >>> my_timeseries = sunpy.timeseries.TimeSeries((data, header))   # doctest: +SKIP

    Headers and units must be either a `dict`, `~collections.OrderedDict` or `~sunpy.util.metadata.MetaDict`.

    * data, header pairs, or data, header units triples, not in tuples

    >>> my_timeseries = sunpy.timeseries.TimeSeries(data, header)  # doctest: +SKIP
    >>> my_timeseries = sunpy.timeseries.TimeSeries(data, header, units)  # doctest: +SKIP

    * File names for files understood by `sunpy.io` and those not

    >>> my_timeseries = sunpy.timeseries.TimeSeries('filename.fits')   # doctest: +SKIP
    >>> my_timeseries = sunpy.timeseries.TimeSeries('filename.fits', source='lyra')  # doctest: +SKIP

    * Multiple files can be combined into one TimeSeries, as long as they are the same source

    >>> my_timeseries = sunpy.timeseries.TimeSeries(['goesfile1.fits', 'goesfile2.fits'],
    ...                                             concatenate=True)  # doctest: +SKIP

    * All fits files in a directory by giving a directory

    >>> my_timeseries = sunpy.timeseries.TimeSeries('local_dir/sub_dir')  # doctest: +SKIP

    * Some regex globs

    >>> my_timeseries = sunpy.timeseries.TimeSeries('eit_*.fits')  # doctest: +SKIP

    * URLs

    >>> my_timeseries = sunpy.timeseries.TimeSeries(url)  # doctest: +SKIP

    * Lists of any of the above

    >>> my_timeseries = sunpy.timeseries.TimeSeries(['file1.fits', 'file2.fits',
    ...                                              'file3.fits', 'directory1/'])  # doctest: +SKIP

    * Any mixture of the above not in a list

    >>> my_timeseries = sunpy.timeseries.TimeSeries((data, header), data2, header2,
    ...                                             'file1.fits', url, 'eit_*.fits')  # doctest: +SKIP
    """

    @staticmethod
    def _read_file(fname, **kwargs):
        """
        Reading a file with `sunpy.io` for automatic source detection.

        Parameters
        ----------
        fname : `str`
            The file path to parse.

        Returns
        -------
        parsed : `bool`
            `True` if file has been read.
        pairs : `list` or `str`
            List of ``(data, header)`` pairs if ``parsed`` is `True`,  ``fname`` if ``parsed`` is `False`.
            `False` if the file is not supported or incorrect.
        """
        if 'source' not in kwargs.keys() or not kwargs['source']:
            try:
                pairs = read_file(os.fspath(fname), **kwargs)

                new_pairs = []
                for pair in pairs:
                    filedata, filemeta = pair
                    if isinstance(filemeta, FileHeader):
                        data = filedata
                        meta = MetaDict(filemeta)
                        new_pairs.append(HDPair(data, meta))
                return True, new_pairs
            except UnrecognizedFileTypeError:
                return False, fname
        else:
            return False, fname

    @staticmethod
    def _validate_meta(meta):
        """
        Validate a meta argument for use as metadata.

        Currently only validates by class.
        """
        if isinstance(meta, astropy.io.fits.header.Header):
            return True
        elif isinstance(meta, sunpy.io.header.FileHeader):
            return True
        elif isinstance(meta, dict):
            return True
        elif isinstance(meta, sunpy.timeseries.TimeSeriesMetaData):
            return True
        else:
            return False

    @staticmethod
    def _validate_units(units):
        """
        Validates the astropy unit-information associated with a
        `~sunpy.timeseries.TimeSeries`.

        Should be a dictionary of some form (but not
        `sunpy.util.metadict.MetaDict`) with only `astropy.units` for
        values.
        """
        result = True

        # It must be a dictionary
        if not isinstance(units, dict) or isinstance(units, MetaDict):
            return False

        for key in units:
            if not isinstance(units[key], u.UnitBase):
                # If this is not a unit then this can't be a valid units dict.
                return False

        # Passed all the tests
        return result

    @staticmethod
    def _from_table(t):
        """
        Extract the data, metadata and units from an astropy table for use in
        constructing a `~sunpy.timeseries.TimeSeries`.

        Parameters
        ----------
        t: `~astropy.table.Table`
            The input table. The datetime column must be the first column or the (single) primary key index.
        """
        table = copy.deepcopy(t)
        # Default the time index to the first column
        index_name = table.colnames[0]
        # Check if another column is defined as the index/primary_key
        if table.primary_key:
            # Check there is only one primary_key/index column
            if len(table.primary_key) == 1:
                table.primary_key[0]
            else:
                raise ValueError("Invalid input Table, TimeSeries doesn't support conversion"
                                 " of tables with more then one index column.")

        # Extract, convert and remove the index column from the input table
        index = table[index_name]
        # Convert if the index is given as an astropy Time object
        if isinstance(index, Time):
            index = index.datetime
        index = pd.to_datetime(index)
        table.remove_column(index_name)

        # Extract the column values from the table
        data = {}
        units = {}
        for colname in table.colnames:
            data[colname] = table[colname]
            units[colname] = table[colname].unit

        # Create a dataframe with this and return
        df = pd.DataFrame(data=data, index=index)
        return df, MetaDict(table.meta), units

    def _parse_args(self, *args, **kwargs):
        """
        Parses an `args` list for data-header pairs. `args` can contain any mixture of the following
        entries:

        * tuples of (data, header, unit) (1)
        * data, header not in a tuple (1)
        * filename, which will be read
        * `pathlib.Path`, which is the filepath to a file to be read.
        * directory, from which all files will be read
        * glob, from which all files will be read
        * url, which will be downloaded and read
        * lists containing any of the above.

        (1) header/unit are optional and in either order, but data should be the first entry in each group.

        Examples
        --------
        self._parse_args(data, header,
                         (data, header),
                         ['file1', 'file2', 'file3'],
                         'file4',
                         'directory1',
                         '*.fits')
        """
        data_header_unit_tuples = list()
        data_header_pairs = list()
        already_timeseries = list()
        filepaths = list()

        # Account for nested lists of items. Simply outputs a single list of
        # items, nested lists are expanded to element level.
        args = expand_list(args)

        # For each of the arguments, handle each of the cases
        i = 0
        while i < len(args):
            arg = args[i]

            # Data-header pair in a tuple
            if (isinstance(arg, (np.ndarray, Table, pd.DataFrame))):
                # and self._validate_meta(args[i+1])):
                # Assume a Pandas Dataframe is given
                data = arg
                units = OrderedDict()
                meta = MetaDict()

                # Convert the data argument into a Pandas DataFrame if needed.
                if isinstance(data, Table):
                    # We have an Astropy Table:
                    data, meta, units = self._from_table(data)
                elif isinstance(data, np.ndarray):
                    # We have a numpy ndarray. We assume the first column is a dt index
                    data = pd.DataFrame(data=data[:, 1:], index=Time(data[:, 0]))

                # If there are 1 or 2 more arguments:
                for _ in range(2):
                    if (len(args) > i+1):
                        # If that next argument isn't data but is metaddata or units:
                        if not isinstance(args[i+1], (np.ndarray, Table, pd.DataFrame)):
                            if self._validate_units(args[i+1]):
                                units.update(args[i+1])
                                i += 1  # an extra increment to account for the units
                            elif self._validate_meta(args[i+1]):
                                # if we have an astropy.io FITS header then convert
                                # to preserve multi-line comments
                                if isinstance(args[i+1], astropy.io.fits.header.Header):
                                    args[i+1] = MetaDict(sunpy.io.header.FileHeader(args[i+1]))
                                if isinstance(args[i+1], sunpy.timeseries.TimeSeriesMetaData):
                                    for j in args[i+1].metas:
                                        meta.update(j)
                                else:
                                    meta.update(args[i+1])
                                i += 1  # an extra increment to account for the meta

                # Add a 3-tuple for this TimeSeries.
                data_header_unit_tuples.append((data, meta, units))

            # Filepath
            elif (isinstance(arg, str) and
                  os.path.isfile(os.path.expanduser(arg))):

                path = os.path.expanduser(arg)
                result = self._read_file(path, **kwargs)
                data_header_pairs, filepaths = _apply_result(data_header_pairs, filepaths, result)

            # pathlib.Path
            elif (isinstance(arg, os.PathLike)):
                result = self._read_file(arg.expanduser())
                data_header_pairs, filepaths = _apply_result(data_header_pairs, filepaths, result)

            # Directory
            elif (isinstance(arg, str) and
                  os.path.isdir(os.path.expanduser(arg))):

                path = os.path.expanduser(arg)
                files = [os.path.join(path, elem) for elem in os.listdir(path)]
                for afile in files:
                    # returns a boolean telling us if it were read and either a
                    # tuple or the original filepath for reading by a source
                    result = self._read_file(afile, **kwargs)
                    data_header_pairs, filepaths = _apply_result(data_header_pairs, filepaths,
                                                                 result)

            # Glob
            elif isinstance(arg, str) and '*' in arg:

                files = glob.glob(os.path.expanduser(arg))
                for afile in files:
                    # returns a boolean telling us if it were read and either a
                    # tuple or the original filepath for reading by a source
                    result = self._read_file(afile, **kwargs)
                    data_header_pairs, filepaths = _apply_result(data_header_pairs, filepaths,
                                                                 result)

            # Already a TimeSeries
            elif isinstance(arg, GenericTimeSeries):
                already_timeseries.append(arg)

            # A URL
            elif (isinstance(arg, str) and
                  _is_url(arg)):
                url = arg
                path = download_file(url, get_and_create_download_dir())
                result = self._read_file(path, **kwargs)
                data_header_pairs, filepaths = _apply_result(data_header_pairs, filepaths, result)

            else:
                raise NoMatchError("File not found or invalid input")
            i += 1

        # TODO:
        # In the end, if there are already TimeSeries it should be put in the
        # same order as the input, currently they are not.
        return data_header_unit_tuples, data_header_pairs, already_timeseries, filepaths

    def __call__(self, *args, silence_errors=False, **kwargs):
        """
        Method for running the factory. Takes arbitrary arguments and keyword
        arguments and passes them to a sequence of pre-registered types to
        determine which is the correct `~sunpy.timeseries.TimeSeries` source
        type to build.

        Arguments args and kwargs are passed through to the validation function and to the constructor for the final type.
        For `~sunpy.timeseries.TimeSeries` types, validation function must take a data-header pair as an argument.

        Parameters
        ----------
        silence_errors : `bool`, optional
            If set, ignore data-header pairs which cause an exception.
            Defaults to `False`.

        Notes
        -----
        Extra keyword arguments are passed through to `sunpy.io.read_file` such as `memmap` for FITS files.
        """
        (data_header_unit_tuples, data_header_pairs,
         already_timeseries, filepaths) = self._parse_args(*args, **kwargs)

        new_timeseries = list()

        # The filepaths for unreadable files
        for filepath in filepaths:
            try:
                new_ts = self._check_registered_widgets(filepath=filepath, **kwargs)
                new_timeseries.append(new_ts)
            except (NoMatchError, MultipleMatchError, ValidationFunctionError):
                if not silence_errors:
                    raise
            except Exception:
                raise

        # data_header_pairs is a list of HDUs as read by sunpy.io
        # For each set of HDus find the matching class and read the
        # data_header_unit_tuples by calling the _parse_hdus method
        # of the class.
        for pairs in data_header_pairs:
            # Pairs may be x long where x is the number of HDUs in the file.
            headers = [pair.header for pair in pairs]

            types = []
            for header in headers:
                try:
                    match = self._get_matching_widget(meta=header, **kwargs)
                    if not match == GenericTimeSeries:
                        types.append(match)
                except (MultipleMatchError, NoMatchError):
                    continue

            if not types:
                # If no specific classes have been found we can read the data
                # if we only have one data header pair:
                if len(pairs) == 1:
                    already_timeseries.append(GenericTimeSeries(pairs[0]._data,
                                                                pairs[0].header))
                else:
                    raise NoMatchError("Input read by sunpy.io can not find a "
                                       "matching class for reading multiple HDUs")
            if len(set(types)) > 1:
                raise MultipleMatchError("Multiple HDUs return multiple matching classes.")

            cls = types[0]

            data_header_unit_tuples.append(cls._parse_hdus(pairs))

        # Loop over each registered type and check to see if WidgetType
        # matches the arguments.  If it does, use that type
        for triple in data_header_unit_tuples:
            data, header, units = triple
            # Make a MetaDict from various input types
            meta = header
            if isinstance(meta, astropy.io.fits.header.Header):
                meta = sunpy.io.header.FileHeader(meta)
            meta = MetaDict(meta)

            try:
                new_ts = self._check_registered_widgets(data=data, meta=meta,
                                                        units=units, **kwargs)
                new_timeseries.append(new_ts)
            except (NoMatchError, MultipleMatchError, ValidationFunctionError):
                if not silence_errors:
                    raise
            except Exception:
                raise

        new_timeseries += already_timeseries

        # Concatenate the timeseries into one if specified.
        concatenate = kwargs.get('concatenate', False)
        if concatenate:
            # Merge all these timeseries into one.
            full_timeseries = new_timeseries.pop(0)
            for timeseries in new_timeseries:
                full_timeseries = full_timeseries.concatenate(timeseries)

            new_timeseries = [full_timeseries]

        # Sanitize any units OrderedDict details
        for timeseries in new_timeseries:
            timeseries._sanitize_units()

        # Only return single time series, not in a list if we only have one.
        if len(new_timeseries) == 1:
            return new_timeseries[0]
        return new_timeseries

    def _get_matching_widget(self, **kwargs):
        candidate_widget_types = list()

        for key in self.registry:
            # Call the registered validation function for each registered class
            if self.registry[key](**kwargs):
                candidate_widget_types.append(key)

        n_matches = len(candidate_widget_types)

        if n_matches == 0:
            if self.default_widget_type is None:
                raise NoMatchError("No types match specified arguments and no default is set.")
            else:
                candidate_widget_types = [self.default_widget_type]
        elif n_matches > 1:
            raise MultipleMatchError("Too many candidate types identified ({})."
                                     "Specify enough keywords to guarantee unique type "
                                     "identification.".format(n_matches))

        # Only one suitable source class is found
        return candidate_widget_types[0]

    def _check_registered_widgets(self, **kwargs):
        """
        Checks the (instrument) source(s) that are compatible with this given
        file/data.

        Only if exactly one source is compatible will a
        `~sunpy.timeseries.TimeSeries` be returned.
        """
        WidgetType = self._get_matching_widget(**kwargs)

        # Dealing with the fact that timeseries filetypes are less consistent
        # (then maps), we use a _parse_file() method embedded into each
        # instrument subclass.
        filepath = kwargs.pop('filepath', None)
        data = kwargs.pop('data', None)
        meta = kwargs.pop('meta', None)
        units = kwargs.pop('units', None)
        if filepath:
            data, meta, units = WidgetType._parse_file(filepath)

        # Now return a TimeSeries from the given file.
        return WidgetType(data, meta, units, **kwargs)


def _apply_result(data_header_pairs, filepaths, result):
    read, result = result
    if read:
        data_header_pairs.append(result)
    else:
        filepaths.append(result)

    return data_header_pairs, filepaths


def _is_url(arg):
    try:
        urlopen(arg)
    except Exception:
        return False
    return True


class InvalidTimeSeriesInput(ValueError):
    """
    Exception to raise when input variable is not a
    `~sunpy.timeseries.TimeSeries` instance and does not point to a valid
    TimeSeries input file.
    """


class InvalidTimeSeriesType(ValueError):
    """
    Exception to raise when an invalid type of timeseries is requested with
    `~sunpy.timeseries.TimeSeries`.
    """


class NoTimeSeriesFound(ValueError):
    """
    Exception to raise when input does not point to any valid
    `~sunpy.timeseries.TimeSeries` or files.
    """


TimeSeries = TimeSeriesFactory(registry=GenericTimeSeries._registry,
                               default_widget_type=GenericTimeSeries,
                               additional_validation_functions=['is_datasource_for'])
