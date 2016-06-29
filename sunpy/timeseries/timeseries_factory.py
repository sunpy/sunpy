from __future__ import absolute_import, division, print_function

__authors__ = ["Alex Hamilton, Russell Hewett, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

import warnings
import os
import glob
from collections import OrderedDict

import numpy as np
import pandas as pd
import astropy.io.fits
from astropy.table import Table

import sunpy
from sunpy.timeseries.timeseriesbase import GenericTimeSeries, TIMESERIES_CLASSES
from sunpy.map.header import MapMeta

from sunpy.io.file_tools import read_file, UnrecognizedFileTypeError
from sunpy.io.header import FileHeader

from sunpy.util.net import download_file
from sunpy.util import expand_list

from sunpy.util.datatype_factory_base import BasicRegistrationFactory
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.datatype_factory_base import MultipleMatchError
from sunpy.util.datatype_factory_base import ValidationFunctionError
from sunpy.extern import six

from sunpy.extern.six.moves.urllib.request import urlopen

__all__ = ['TimeSeries', 'TimeSeriesFactory']

class TimeSeriesFactory(BasicRegistrationFactory):
    """
    TimeSeries(*args, **kwargs)

    TimeSeries factory class.  Used to create a variety of TimeSeries objects.
    Valid time series types are specified by registering them with the factory.

    Parameters
    ----------

    source : `str`, optional
        A string to select the observational source of the data, currently
        necessary to define how files should be read for all instruments.

    concatenate : `bool`, optional, default:False
        If set, combine any resulting list of TimeSeries objects into a single
        TimeSeries, using successive concatenate methods.

    Examples
    --------
    >>> import sunpy.timeseries
    >>> sunpy.data.download_sample_data(overwrite=False)   # doctest: +SKIP
    >>> import sunpy.data.sample
    >>> my_timeseries = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_LIGHTCURVE)

    The SunPy TimeSeries factory accepts a wide variety of inputs for creating time series

    * Preloaded tuples of (data, header) pairs or (data, header, units)

    >>> my_timeseries = sunpy.timeseries.TimeSeries((data, header))   # doctest: +SKIP

    headers and units are some base of `dict` or `~collections.OrderedDict`.

    * data, header pairs, or data, header units triples, not in tuples

    >>> my_timeseries = sunpy.timeseries.TimeSeries(data, header)
    >>> my_timeseries = sunpy.timeseries.TimeSeries(data, header, units)

    * File names for files understood by sunpy.io and those not

    >>> my_timeseries = sunpy.timeseries.TimeSeries('filename.fits')   # doctest: +SKIP
    >>> my_timeseries = sunpy.timeseries.TimeSeries('filename.fits', source='lyra')

    * Multiple files can be combined into one TimeSeries, as long as they are
    the same source   # doctest: +SKIP

    >>> my_timeseries = sunpy.timeseries.TimeSeries(['goesfile1.fits', 'goesfile2.fits'], concatenate=True)

    * All fits files in a directory by giving a directory

    >>> my_timeseries = sunpy.timeseries.TimeSeries('local_dir/sub_dir')   # doctest: +SKIP

    * Some regex globs

    >>> my_timeseries = sunpy.timeseries.TimeSeries('eit_*.fits')   # doctest: +SKIP

    * URLs

    >>> my_timeseries = sunpy.timeseries.TimeSeries(url)   # doctest: +SKIP

    * Lists of any of the above

    >>> my_timeseries = sunpy.timeseries.TimeSeries(['file1.fits', 'file2.fits', 'file3.fits', 'directory1/'])   # doctest: +SKIP

    * Any mixture of the above not in a list

    >>> my_timeseries = sunpy.timeseries.TimeSeries((data, header), data2, header2, 'file1.fits', url, 'eit_*.fits')   # doctest: +SKIP
    """

    def _read_file(self, fname, **kwargs):
        """ Read in a file name and return the list of (data, meta, unit) tuples in
            that file. """
        
        # ToDo: implement this for the TimeSeries using either Pandas or AstroPy parser.
        
        # File gets read here.  This needs to be generic enough to seamlessly
        #call a fits file or a jpeg2k file, etc
        try:
            pairs = read_file(fname, **kwargs)

            new_pairs = []
            for pair in pairs:
                filedata, filemeta = pair
                assert isinstance(filemeta, FileHeader)
                # ToDo Validate data before adding it.
                #if len(np.shape(filedata)) > 1:
                if True:
                    data = filedata
                    meta = MapMeta(filemeta)
                    new_pairs.append((data, meta))
            return True, new_pairs
        except UnrecognizedFileTypeError:
            return False, fname


    def _validate_meta(self, meta):
        """
        Validate a meta argument.
        """
        if isinstance(meta, astropy.io.fits.header.Header):
            return True
        elif isinstance(meta, dict):
            return True
        else:
            return False

    def _validate_units(self, units, **kwargs):
        """
        Validates the astropy unit-information associated with a TimeSeries.
        """
    
        warnings.simplefilter('always', Warning)
        
        result = True
        
        # It must be a dictionary
        if not isinstance(units, dict):
            return False
        
        for key in units:
            if not isinstance(units[key], astropy.units.UnitBase):
                # If this is not a unit then this can't be a valid units dict.
                warnings.warn("Invalid unit given for \""+str(key)+"\"", Warning)
                return False
        
        # Passed all the tests
        return result

    def _parse_args(self, *args, **kwargs):
        """
        Parses an args list for data-header pairs.  args can contain any
        mixture of the following entries:
        * tuples of data,header
        * data, header not in a tuple
        * filename, which will be read
        * directory, from which all files will be read
        * glob, from which all files will be read
        * url, which will be downloaded and read
        * lists containing any of the above.

        Example
        -------
        self._parse_args(data, header,
                         (data, header),
                         ['file1', 'file2', 'file3'],
                         'file4',
                         'directory1',
                         '*.fits')

        """

        data_header_unit_tuples = list()
        already_timeseries = list()
        filepaths = list()
        
        # Take source kwarg if defined
        source = kwargs.get('source', None)

        # Account for nested lists of items
        args = expand_list(args)

        # For each of the arguments, handle each of the cases
        i = 0
        while i < len(args):

            arg = args[i]

            # Data-header pair in a tuple
            if ((type(arg) in [tuple, list]) and
                (len(arg) == 2 or len(arg) == 3) and
                isinstance(arg[0], (np.ndarray, Table, pd.DataFrame)) and
                self._validate_meta(arg[1])):
                # Assume a Pandas Dataframe is given.
                data = arg[0]
                units = OrderedDict({})
                # Convert the data argument into a Pandas DataFrame if needed.
                if isinstance(data, Table):
                    # We have an AstroPy Table:
                    data = arg[0].to_pandas()
                    
                    # Extract Units from table:
                    for colname in arg[0].colnames:
                        # Only add the unit if specified.
                        if arg[0][colname].unit:
                            units.update({colname:arg[0][colname].unit})
                    
                elif isinstance(data, np.ndarray):
                    # We have a numpy ndarray:
                    data = pd.DataFrame(data=arg[0])
                    # ToDo: should this include an index? Maybe use the first column?
                
                # The second argument will be the metadata/header.
                meta = OrderedDict(arg[1])
                
                # Check if we're given a third argument for units
                if (len(arg) == 3) and self._validate_meta(arg[2]):
                    # This combines with values gathered from an input Table.
                    units.update(arg[2])
                
                # Add a 3-tuple for this TimeSeries.
                data_header_unit_tuples.append((data, meta, units))

            # Data-header pair not in a tuple
            elif (isinstance(arg, (np.ndarray, Table, pd.DataFrame)) and
                  self._validate_meta(args[i+1])):
                # Assume a Pandas Dataframe is given.
                data = arg
                units = OrderedDict()
                # Convert the data argument into a Pandas DataFrame if needed.
                if isinstance(data, Table):
                    # We have an AstroPy Table:
                    data = arg.to_pandas()
                    # Extract Units from table:
                    for colname in arg.colnames:
                        # Only add the unit if specified.
                        if arg[colname].unit:
                            units.update({colname:arg[colname].unit})
                elif isinstance(data, np.ndarray):
                    # We have a numpy ndarray:
                    data = pd.DataFrame(data=arg)
                    # ToDo: should this include an index? Maybe use the first column?
                    
                # The second argument will be the metadata/header.
                meta = OrderedDict(args[i+1])
                
                # Check if we're given a third argument for units
                if (len(args) > i+2) and self._validate_units(args[i+2]):
                    units.update(args[i+2])
                    i += 1 # an extra increment to account for the units
                
                # Add a 3-tuple for this TimeSeries.
                data_header_unit_tuples.append((data, meta, units))
                i += 1 # an extra increment to account for the header
                
            
            # Filepath
            elif (isinstance(arg,six.string_types) and
                  os.path.isfile(os.path.expanduser(arg))):
                path = os.path.expanduser(arg)
                
                # Only read using generic _read_file() if no source defined.
                if not source:
                    temp = self._read_file(path, **kwargs)
                else:
                    temp = ( False, path )
                
                # If the file was successfully read in:
                if temp[0]:
                    pairs = temp[1]
                    data_header_unit_tuples += pairs
                else:
                    # Unsuccessfully read files are listed to be read by the
                    # source instrument specific file_paser().
                    filepaths.append(temp[1])

            # Directory
            elif (isinstance(arg,six.string_types) and
                  os.path.isdir(os.path.expanduser(arg))):
                path = os.path.expanduser(arg)
                files = [os.path.join(path, elem) for elem in os.listdir(path)]
                for afile in files:
                    data_header_unit_tuples += self._read_file(afile, **kwargs)

            # Glob
            elif (isinstance(arg,six.string_types) and '*' in arg):
                files = glob.glob( os.path.expanduser(arg) )
                for afile in files:
                    data_header_unit_tuples += self._read_file(afile, **kwargs)

            #### Potentially unnecessary, as there won't be a time series cube.
            # Already a TimeSeries
            elif isinstance(arg, GenericTimeSeries):
                already_timeseries.append(arg)
            
            # A URL
            elif (isinstance(arg,six.string_types) and
                  _is_url(arg)):
                default_dir = sunpy.config.get("downloads", "download_dir")
                url = arg
                path = download_file(url, default_dir)
                pairs = self._read_file(path, **kwargs)
                data_header_unit_tuples += pairs

            else:
                raise ValueError("File not found or invalid input")

            i += 1
        #TODO:
        # In the end, if there are already TimeSeries it should be put in the
        # same order as the input, currently they are not.
        return data_header_unit_tuples, already_timeseries, filepaths


    def __call__(self, *args, **kwargs):
        """ Method for running the factory. Takes arbitrary arguments and
        keyword arguments and passes them to a sequence of pre-registered types
        to determine which is the correct TimeSeries source type to build.

        Arguments args and kwargs are passed through to the validation
        function and to the constructor for the final type.  For TimeSeries
        types, validation function must take a data-header pair as an argument.

        Parameters
        ----------

        silence_errors : `bool`, optional
            If set, ignore data-header pairs which cause an exception.

        Notes
        -----
        Extra keyword arguments are passed through to `sunpy.io.read_file` such
        as `memmap` for FITS files.
        """

        # Hack to get around Python 2.x not backporting PEP 3102.
        silence_errors = kwargs.pop('silence_errors', False)

        data_header_unit_tuples, already_timeseries, filepaths = self._parse_args(*args, **kwargs)

        new_timeseries = list()
        
        # The filepaths for unreadable files
        for filepath in filepaths:
            try:
                new_ts = self._check_registered_widgets(filepath=filepath, **kwargs)
            except (NoMatchError, MultipleMatchError, ValidationFunctionError):
                if not silence_errors:
                    raise
            except:
                raise

            new_timeseries.append(new_ts)

        # Loop over each registered type and check to see if WidgetType
        # matches the arguments.  If it does, use that type.
        for pair in data_header_unit_tuples:
            data, header, units = pair
            meta = MapMeta(header)

            try:
                new_ts = self._check_registered_widgets(data=data, meta=meta, units=units, **kwargs)
            except (NoMatchError, MultipleMatchError, ValidationFunctionError):
                if not silence_errors:
                    raise
            except:
                raise

            new_timeseries.append(new_ts)

        new_timeseries += already_timeseries

        # Concatenate the timeseries into one if specified.
        concatenate = kwargs.get('concatenate', False)
        if concatenate:
            # Merge all these timeseries into one.
            # ToDo: consider metadata output carfully.
            full_timeseries = new_timeseries.pop(0)
            for timeseries in new_timeseries:
                full_timeseries = full_timeseries.concatenate(timeseries)
            
            new_timeseries = [ full_timeseries ]

        # Sanitize any units OrderedDict details
        for timeseries in new_timeseries:
            timeseries._sanitize_units()
        
        # Only return single time series, not in a list if we only have one.
        if len(new_timeseries) == 1:
            return new_timeseries[0]
        
        return new_timeseries

    def _check_registered_widgets(self, **kwargs):
        """Checks the (instrument) source/s that are compatible with this given file/data.
        Only if exactly one source is compatible will a TimeSeries be returned."""
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
            raise MultipleMatchError("Too many candidate types identified ({0}).  Specify enough keywords to guarantee unique type identification.".format(n_matches))

        # Only one suitable source class is found
        WidgetType = candidate_widget_types[0]

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

def _is_url(arg):
    try:
        urlopen(arg)
    except:
        return False
    return True

class InvalidTimeSeriesInput(ValueError):
    """Exception to raise when input variable is not a TimeSeries instance and
    does not point to a valid TimeSeries input file."""
    pass

class InvalidTimeSeriesType(ValueError):
    """Exception to raise when an invalid type of time series is requested with
    TimeSeries."""
    pass

class NoTimeSeriesFound(ValueError):
    """Exception to raise when input does not point to any valid time series or
    files."""
    pass

TimeSeries = TimeSeriesFactory(default_widget_type=GenericTimeSeries,
                 additional_validation_functions=['is_datasource_for'])
TimeSeries.registry = TIMESERIES_CLASSES