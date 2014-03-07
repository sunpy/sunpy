import os
import glob
import urllib2
import datetime

import numpy as np

from astropy.utils.misc import isiterable

import sunpy
from sunpy.database.tables import DatabaseEntry
from sunpy.util.net import download_file
from sunpy.time import parse_time
from sunpy.util import expand_list
from sunpy.io.file_tools import read_file, UnrecognizedFileTypeError
from sunpy.io.header import FileHeader

from sunpy.map.header import MapMeta

from sunpy.util.datatype_factory_base import BasicRegistrationFactory
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.datatype_factory_base import MultipleMatchError
from sunpy.util.datatype_factory_base import ValidationFunctionError

from . lightcurve import GenericLightCurve

def _is_url(arg):
    try:
        urllib2.urlopen(arg)
    except ValueError:
        return False
    return True

def _is_time_string(arg):
    try:
        parse_time(arg)
    except (ValueError, TypeError):
        return False
    return True

def _is_time_date(arg):
    """
    Checks if time is timerange time string or time strinf pair or datetime.date
    """
    if isiterable(arg) and len(arg) <= 2 and all([_is_time_string(x) for x in arg]) or \
       isinstance(arg, sunpy.time.TimeRange) or isinstance(arg, datetime.date) \
       or _is_time_string(arg):
           return True

    return False

def _is_file(arg):
    return isinstance(arg,basestring) and os.path.isfile(os.path.expanduser(arg))

class LightCurveFactory(BasicRegistrationFactory):
    """
    LightCurve factory class.
    Used to create a variety of LightCurve objects.
    Valid LightCurve types are specified by registering them with the factory.

    LightCurves can be created from a variety of inputs and then truncated by
    timerange.

    Examples
    --------
    >>> import sunpy.lightcurve
    >>> lc = sunpy.lightcurve.LightCurve('g15_xrs_2s_20120101_20120102.csv',
    ...                                  source='goes')

    * A file that can be read by sunpy.io

    >>> lc = sunpy.lightcurve.LightCurve('lyra_20110810-000000_lev2_std.fits')

    * A file that can't be auto-detected

    >>> lc = sunpy.lightcurve.LightCurve('g15_xrs_2s_20120101_20120102.csv',
    ...                                  source='goes')

    * Generate multiple LightCurves

    >>> lc = sunpy.lightcurve.LightCurve('g15_xrs_2s_20120101_20120102.csv',
    ...                                  'g15_xrs_2s_20031021_20120102.csv',
    ...                                  source='goes')

    * Generate multiple LightCurves from different sources

    >>> lc = sunpy.lightcurve.LightCurve('lyra_20110810-000000_lev2_std.fits',
    ...                                  'g15_xrs_2s_20031021_20120102.csv',
    ...                                  source=['lyra','goes'])

    Note: if you only specify one source, LightCurve will assume that all other
    arguments are of that type, so if not you will have to manually override.

    * Multiple Files to be read into one LightCurve

    >>> lc = sunpy.lightcurve.LightCurve(['lyra_20120101-000000_lev2_std.fits',
    ...                                   'lyra_20120102-000000_lev2_std.fits'])

    * Multiple Files to be read into multiple LightCurves

    >>> lc = sunpy.lightcurve.LightCurve(['lyra_20120101-000000_lev2_std.fits',
    ...                                   'lyra_20120102-000000_lev2_std.fits'],
    ...                                  ['lyra_20130521-000000_lev2_std.fits',
    ...                                   'lyra_20130521-000000_lev2_std.fits'])

    * Create a LightCurve from a file and then truncate to a time range.

    >>> lc = sunpy.lightcurve.LightCurve('lyra_20120101-000000_lev2_std.fits',
    ...                                   timerange=["2012/1/1T00:00:00", "2012/1/1T12:00:00"])

    Note: Only one time range can be specified so it will truncate all files to
    that time range.



    """
    def _read_file(self, fname, source=None, **kwargs):
        try:
            return self._io_read_file(fname, **kwargs)
        except UnrecognizedFileTypeError:
            return [self._check_registered_widgets(source=source)[0]._read_file(fname)]

    def _io_read_file(self, fname, **kwargs):
        """
        Read in a file name and return the list of (data, meta) pairs in that
        file.
        """
        pairs = read_file(fname, fold_blank_hdu=self.fold_blank_hdu, **kwargs)

        new_pairs = []
        for ind, pair in enumerate(pairs):
            filedata, filemeta = pair
            assert isinstance(filemeta, FileHeader)
            data = filedata
            if data is None and self.fold_blank_hdu: #Skip blank
                continue
            meta = MapMeta(filemeta)
            new_pairs.append((data, meta))

        return new_pairs

    def __call__(self, *args, **kwargs):
        """
        This is a tester for the input parsing of the new LightCurveFactory object

        Unlike map one lc can be made up of many files
        """
        # Hack to get around Python 2.x not backporting PEP 3102.
        silence_errors = kwargs.pop('silence_errors', False)
        self.fold_blank_hdu = kwargs.pop('fold_blank_hdu', True)
        self.source = kwargs.pop('source', None)
        self.timerange = kwargs.pop('timerange', None)
        if not isiterable(self.source) or isinstance(self.source, basestring):
            #Make source a list that is the number of args long.
            if len(args): #Don't multiply by 0
                self.source = [self.source]*len(args)
            else:
                self.source = [self.source]

        elif len(self.source) != len(args) and self.timerange is None:
            raise ValueError("""Source argument must be a string,
            or a list the same length as the input unless timerange is given""")

        data_header_pairs = list() #each lc is in a list in here
        already_lcs = list()
        args = list(args)

        #Parse timerange, so that self.timerange is always a TimeRange
        if self.timerange is not None:
            _one_day = datetime.timedelta(hours=23, minutes=59, seconds=59)
            #Is a time string pari
            if isiterable(self.timerange) and len(self.timerange) == 2 and \
               all([_is_time_string(x) for x in self.timerange]):
                self.timerange = sunpy.time.TimeRange(self.timerange[0],self.timerange[1])

            #Is one time string
            elif _is_time_string(self.timerange):
                #If one time is specified, take it as the start time and run
                #till the end of that day.
                end = datetime.datetime.combine(parse_time(self.timerange).date(),
                                                datetime.datetime.min.time()) + _one_day
                self.timerange = sunpy.time.TimeRange(self.timerange, end)

            #Is a datetime.date
            elif isinstance(self.timerange, datetime.date):
                self.timerange = sunpy.time.TimeRange(self.timerange, self.timerange + _one_day)

            #Is timerange
            elif isinstance(self.timerange, sunpy.time.TimeRange):
                pass

            else:
                raise ValueError("""timerange must be one of: a TimeRange object,
                a datetime.date object, one time / date string or a pair of time / date strings.""")

            if not len(args) == len(self.source):
                for asource in self.source[len(args):]:
                    args.append(self._check_registered_widgets(source=asource)[0]._get_url_from_timerange(self.timerange, **kwargs))

        # For each of the arguments, handle each of the cases
        i = 0
        while i < len(args):
            arg = args[i]

            dhp, lcs = self._process_single_lc_args(arg, source=self.source[i])
            data_header_pairs.append(dhp)
            already_lcs.append(lcs)
            i += 1

        #At this point we have finished the input processing and we now need
        # to create LightCurves from the (data,header) pairs and then concat
        # them with any existing LCs in that group.
        new_lc_sets = list()
        for i, lc_dhp in enumerate(data_header_pairs):
            new_lcs = list()
            # Loop over each registered type and check to see if WidgetType
            # matches the arguments.  If it does, use that type.
            for pair in lc_dhp:
                data, header = pair
                meta = header
                meta = MapMeta(header) #TODO: LC this

                try:
                    kwergs = {'source':self.source[i]}
                    kwergs.update(kwargs)
                    new_lc = self._get_registered_widget(data, meta, **kwergs)
                except (NoMatchError, MultipleMatchError, ValidationFunctionError):
                    if not silence_errors:
                        raise
                except:
                    raise

                new_lcs.append(new_lc)

            new_lc_sets.append(new_lcs)

        #This should give a list of lists containg all lc's to be concatenated into
        #one lightcurve
        all_lc_sets = list()
        for old,new in zip(new_lc_sets, already_lcs):
            all_lc_sets.append(old + new)

        lightcurves = self._concatenate_lightcurves(all_lc_sets)

        #If a time range is specified then cut down the lc now:
        if self.timerange:
            lightcurves = [lightcurve.truncate(self.timerange) for lightcurve in lightcurves]

        if len(lightcurves) == 1:
            return lightcurves[0]
        return lightcurves

    def _concatenate_lightcurves(self, lc_sets):
        lightcurves = list()

        for lc_set in lc_sets:
            concatted = lc_set.pop(0)
            for lc in lc_set:
                concatted += lc #TDO: implement __add__ in GenericLightCurve
            lightcurves.append(concatted)

        return lightcurves


    def _process_single_lc_args(self, *args, **kwargs):
        """
        This function processes all the input to make one lc object,
        more like the processing in Map
        """
        source = kwargs.pop('source', None)
        data_header_pairs = list() #each lc is in a list in here
        already_lcs = list()

        args = expand_list(args)
        # For each of the arguments, handle each of the cases
        i = 0
        while i < len(args):
            arg = args[i]

            # Data-header pair in a tuple #TODO: pretty sure this won't work here
            if ((type(arg) in [tuple, list]) and
                 len(arg) == 2 and
                 isinstance(arg[0],np.ndarray) and
                 isinstance(arg[1],dict)):
                data_header_pairs.append(arg)

            # Data-header pair not in a tuple
            elif (isinstance(arg, np.ndarray) and
                  isinstance(args[i+1],dict)):
                pair = (args[i], args[i+1])
                data_header_pairs.append(pair)
                i += 1 # an extra increment to account for the data-header pairing

            # File name
            elif _is_file(arg):
                path = os.path.expanduser(arg)
                pairs = self._read_file(path, source=source, **kwargs)
                data_header_pairs += pairs

            # Directory
            elif (isinstance(arg,basestring) and
                  os.path.isdir(os.path.expanduser(arg))):
                path = os.path.expanduser(arg)
                files = [os.path.join(path, elem) for elem in os.listdir(path)]
                for afile in files:
                    data_header_pairs += self._read_file(afile, source=source, **kwargs)

            # Glob
            elif (isinstance(arg,basestring) and '*' in arg):
                files = glob.glob( os.path.expanduser(arg) )
                for afile in files:
                    data_header_pairs += self._read_file(afile, source=source, **kwargs)

            # Already a Map
            elif isinstance(arg, GenericLightCurve):
                already_lcs.append(arg)

            # A URL
            elif isinstance(arg,basestring) and _is_url(arg):
                default_dir = sunpy.config.get("downloads", "download_dir") #TODO: target filename and dir
                url = arg
                path = download_file(url, default_dir)
                pairs = self._read_file(path, source=source, **kwargs)
                data_header_pairs += pairs

            # A database Entry
            elif isinstance(arg, DatabaseEntry):
                data_header_pairs += self._read_file(arg.path, source=source, **kwargs)

            else:
                raise ValueError("File not found or invalid input")
            i += 1

        return data_header_pairs, already_lcs

    def _check_registered_widgets(self, data=None, meta=None, source=None, **kwargs):
        candidate_widget_types = list()

        for key in self.registry:
            # Call the registered validation function for each registered class
            if self.registry[key](data, meta, source=source, **kwargs):
                candidate_widget_types.append(key)

        n_matches = len(candidate_widget_types)

        if n_matches == 0:
            if self.default_widget_type is None:
                raise NoMatchError("No types match specified arguments and no default is set.")
            else:
                return [self.default_widget_type]
        elif n_matches > 1:
            raise MultipleMatchError("""Too many candidate types idenfitied ({0}).
            Specify enough keywords to guarantee unique type identification.""".format(n_matches))

        return candidate_widget_types

    def _get_registered_widget(self, data, meta, source=None, **kwargs):
        candidate_widget_types = self._check_registered_widgets(data, meta, source=source)
        # Only one is found
        WidgetType = candidate_widget_types[0]

        return WidgetType(data, meta, **kwargs)




LightCurve = LightCurveFactory(default_widget_type=GenericLightCurve,
                 additional_validation_functions=['is_datasource_for'])