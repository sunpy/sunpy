"""
This module provies the `~sunpy.timeseries.TimeSeriesFactory` class.
"""
import os
import copy
import pathlib
from collections import OrderedDict
from urllib.request import Request

import numpy as np
import pandas as pd

import astropy
import astropy.io.fits
import astropy.units as u
from astropy.table import Table
from astropy.time import Time

import sunpy
from sunpy.io.file_tools import UnrecognizedFileTypeError, detect_filetype, read_file
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
from sunpy.util.functools import seconddispatch
from sunpy.util.io import is_url, parse_path, possibly_a_path
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
        pairs : `list` or `str`
            List of ``(data, header)`` pairs or ``fname`` if the file is not supported or incorrect.
        """
        if 'source' not in kwargs.keys() or not kwargs['source']:
            try:
                if detect_filetype(fname) == 'cdf':
                    # Put import here to ensure there is no import dependency
                    # on cdflib for TimeSeries
                    from sunpy.io.cdf import read_cdf
                    return read_cdf(os.fspath(fname), **kwargs)
            except UnrecognizedFileTypeError:
                pass

            try:
                pairs = read_file(os.fspath(fname), **kwargs)

                new_pairs = []
                for pair in pairs:
                    filedata, filemeta = pair
                    if isinstance(filemeta, FileHeader):
                        data = filedata
                        meta = MetaDict(filemeta)
                        new_pairs.append(HDPair(data, meta))
                return [new_pairs]
            except UnrecognizedFileTypeError:
                return [fname]
        else:
            return [fname]

    @staticmethod
    def _is_metadata(meta):
        """
        Return `True` if ``meta`` is an object that could store metadata.
        """
        return isinstance(meta, (astropy.io.fits.header.Header,
                                 sunpy.io.header.FileHeader,
                                 dict,
                                 sunpy.timeseries.TimeSeriesMetaData))

    @staticmethod
    def _is_units(units):
        """
        Return `True` if ``units`` is an object that could store units.

        Should be a dictionary of some form (but not `sunpy.util.metadict.MetaDict`)
        with only `astropy.units` for values.
        """
        # Must be a dict and all items must be a unit
        return (isinstance(units, dict)
                and not isinstance(units, MetaDict)
                and all(isinstance(units[key], u.UnitBase) for key in units))

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

    def _parse_meta(self, meta):
        """
        Parse different metadata objects into a MetaDict.
        """
        if isinstance(meta, astropy.io.fits.header.Header):
            meta = MetaDict(sunpy.io.header.FileHeader(meta))
        if isinstance(meta, sunpy.timeseries.TimeSeriesMetaData):
            new_meta = MetaDict()
            for m in meta.metas:
                new_meta.update(m)
            meta = new_meta
        return meta

    def _sanitise_args(self, args):
        """
        Sanitise a list of args so that a single argument corresponds to either:

        - (data, header, units) tuple.
        - path-like `pathlib.Path` (e.g. a filename, directory, glob etc.).
        - `urllib.request.Request`.
        - `GenericTimeSeries`.
        """
        # Account for nested lists of items. Simply outputs a single list of
        # items, nested lists are expanded to element level.
        args = expand_list(args)

        # Sanitise the input so that each 'type' of input corresponds to a different
        # class, so single dispatch can be used later
        i = 0
        while i < len(args):
            arg = args[i]
            if isinstance(arg, (np.ndarray, Table, pd.DataFrame)):
                # Extract data and metadata
                # The next item is data
                data = args[i]
                meta = MetaDict()
                units = OrderedDict()
                if isinstance(data, Table):
                    # We have an Astropy Table:
                    data, new_meta, new_units = self._from_table(data)
                    units.update(new_units)
                    meta.update(new_meta)
                elif isinstance(data, np.ndarray):
                    # We have a numpy ndarray. We assume the first column is a dt index
                    data = pd.DataFrame(data=data[:, 1:], index=Time(data[:, 0]))

                # The next two could be metadata or units
                for _ in range(2):
                    j = i+1
                    if j < len(args):
                        arg = args[j]
                        if self._is_units(arg):
                            units.update(arg)
                            args.pop(j)
                        elif self._is_metadata(arg):
                            meta.update(self._parse_meta(arg))
                            args.pop(j)

                args[i] = (data, meta, units)

            elif isinstance(arg, str) and is_url(arg):
                args[i] = Request(arg)
            elif possibly_a_path(arg):
                args[i] = pathlib.Path(arg)
            i += 1

        return args

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
        args = self._sanitise_args(args)
        all_ts = []
        for arg in args:
            try:
                all_ts += self._parse_arg(arg, **kwargs)
            except (NoMatchError, MultipleMatchError, ValidationFunctionError):
                if self.silence_errors:
                    continue
                else:
                    raise
            except Exception:
                raise

        return all_ts

    @seconddispatch
    def _parse_arg(self, arg, **kwargs):
        """
        Parse a single arg and return a list of timeseries.
        """
        raise NoMatchError("File not found or invalid input")

    @_parse_arg.register(GenericTimeSeries)
    def _parse_ts(self, ts, **kwargs):
        return [ts]

    @_parse_arg.register(Request)
    def _parse_url(self, request, **kwargs):
        path = download_file(request.full_url, get_and_create_download_dir())
        return self._parse_path(pathlib.Path(path), **kwargs)

    @_parse_arg.register(pathlib.Path)
    def _parse_path(self, path, **kwargs):
        results = parse_path(path, self._read_file)
        all_ts = []

        # r can be either a TimeSeries, path, or a data, header pair
        for r in results:
            if isinstance(r, GenericTimeSeries):
                all_ts += [r]
            elif isinstance(r, pathlib.Path):
                all_ts += [self._check_registered_widgets(filepath=r, **kwargs)]
            else:
                pairs = r
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
                        all_ts += [GenericTimeSeries(pairs[0]._data, pairs[0].header)]
                        continue
                    else:
                        raise NoMatchError("Input read by sunpy.io can not find a "
                                           "matching class for reading multiple HDUs")
                if len(set(types)) > 1:
                    raise MultipleMatchError("Multiple HDUs return multiple matching classes.")

                cls = types[0]

                data_header_unit_tuple = cls._parse_hdus(pairs)
                all_ts += self._parse_arg(data_header_unit_tuple)

        return all_ts

    @_parse_arg.register(tuple)
    def _parse_tuple(self, tup, **kwargs):
        data, header, units = tup
        # Make a MetaDict from various input types
        meta = header
        if isinstance(meta, astropy.io.fits.header.Header):
            meta = sunpy.io.header.FileHeader(meta)
        meta = MetaDict(meta)
        return [self._check_registered_widgets(data=data, meta=meta, units=units, **kwargs)]

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
        self.silence_errors = silence_errors
        new_timeseries = self._parse_args(*args, **kwargs)

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
