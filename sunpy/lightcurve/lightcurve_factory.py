import os
import glob
import urllib2

import numpy as np

from astropy.utils.misc import isiterable

import sunpy
from sunpy.database.tables import DatabaseEntry
from sunpy.util.net import download_file
from sunpy.time import parse_time
from sunpy.util import expand_list
from sunpy.io.file_tools import read_file
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
    except: #ValueEror?
        return False
    return True

def _is_time(arg):
    try:
        parse_time(arg)
    except ValueError:
        return False
    return True

def _is_file(arg):
    return isinstance(arg,basestring) and os.path.isfile(os.path.expanduser(arg))

class LightCurveFactory(BasicRegistrationFactory):
    """
    This is the lightcurve factory
    """

    def _read_file(self, fname, **kwargs):
        """ Read in a file name and return the list of (data, meta) pairs in
            that file. """

        # File gets read here.  This needs to be generic enough to seamlessly
        #call a fits file or a jpeg2k file, etc
        pairs = read_file(fname, **kwargs)

        _no_data = list()

        new_pairs = []
        for ind, pair in enumerate(pairs):
            filedata, filemeta = pair
            assert isinstance(filemeta, FileHeader)
            if filedata is None:
                _no_data.append(ind)
            data = filedata
            meta = MapMeta(filemeta)
            new_pairs.append((data, meta))

        for pair in new_pairs:
            if pair[0] is None:

        return new_pairs

    def __call__(self, *args, **kwargs):
        """
        This is a tester for the input parsing of the new LightCurveFactory object

        Unlike one lc can be made up of many files
        """
        # Hack to get around Python 2.x not backporting PEP 3102.
        silence_errors = kwargs.pop('silence_errors', False)
        concat_hdus = kwargs.pop('concatenate_hdus', False)

        data_header_pairs = list() #each lc is in a list in here
        already_lcs = list()

        # For each of the arguments, handle each of the cases
        i = 0
        while i < len(args):
            arg = args[i]

            #A iterable arg
            if isiterable(arg):
                if all([_is_time(x) for x in arg]):
                    if len(arg) != 2:
                        raise ValueError("Times must be specified in pairs")
                    #Do somthing with time
                    pass
                else:
                    dhp, lcs = self._process_single_lc_args(arg)
                    data_header_pairs.append(dhp)
                    already_lcs.append(lcs)
            else:
                dhp, lcs = self._process_single_lc_args(arg)
                data_header_pairs.append(dhp)
                already_lcs.append(lcs)

            i += 1

        #At this point we have finished the input processing and we now need
        # to create LightCurves from the (data,header) pairs and then concat
        # them with any existing LCs in that group.
#        import pdb; pdb.set_trace()

        new_lc_sets = list()
        for lc_dhp in data_header_pairs:
            assert isinstance(lc_dhp, list) #This dosen't have to stay

            new_lcs = list()
            # Loop over each registered type and check to see if WidgetType
            # matches the arguments.  If it does, use that type.
            for pair in lc_dhp:
                data, header = pair
                meta = header
                meta = MapMeta(header) #TODO: LC this

                try:
                    new_lc = self._check_registered_widgets(data, meta, **kwargs)
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
        pass#TODO: Implement

        return lightcurves

    def _concatenate_lightcurves(self, lc_sets):
        lightcurves = list()

        for lc_set in lc_sets:
            concatted = lc_set.pop(1)
#            for lc in lc_set:
#                concatted += lc #TDO: implement __add__ in GenericLightCurve
            lightcurves.append(concatted)

        if len(lightcurves) == 1:
            return lightcurves[0]

        return lightcurves


    def _process_single_lc_args(self, *args, **kwargs):
        """
        This function processes all the input to make one lc object,
        more like the processing in Map
        """
        data_header_pairs = list() #each lc is in a list in here
        already_lcs = list()

        #TODO: Check that this dosen't futz with things
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
                pairs = self._read_file(path, **kwargs)
                data_header_pairs += pairs

            # Directory
            elif (isinstance(arg,basestring) and
                  os.path.isdir(os.path.expanduser(arg))):
                path = os.path.expanduser(arg)
                files = [os.path.join(path, elem) for elem in os.listdir(path)]
                for afile in files:
                    data_header_pairs += self._read_file(afile, **kwargs)

            # Glob
            elif (isinstance(arg,basestring) and '*' in arg):
                files = glob.glob( os.path.expanduser(arg) )
                for afile in files:
                    data_header_pairs += self._read_file(afile, **kwargs)

            # Already a Map
            elif isinstance(arg, GenericLightCurve):
                already_lcs.append(arg)

            # A URL
            elif (isinstance(arg,basestring) and
                  _is_url(arg)):
                default_dir = sunpy.config.get("downloads", "download_dir")
                url = arg
                path = download_file(url, default_dir)
                pairs = self._read_file(path, **kwargs)
                data_header_pairs += pairs

            # A database Entry
            elif isinstance(arg, DatabaseEntry):
                data_header_pairs += self._read_file(arg.path, **kwargs)

            else:
                raise ValueError("File not found or invalid input")
            i += 1

            return data_header_pairs, already_lcs

    def _check_registered_widgets(self, data, meta, **kwargs):

        candidate_widget_types = list()

        for key in self.registry:

            # Call the registered validation function for each registered class
            if self.registry[key](data, meta, **kwargs):
                candidate_widget_types.append(key)

        n_matches = len(candidate_widget_types)

        if n_matches == 0:
            if self.default_widget_type is None:
                raise NoMatchError("No types match specified arguments and no default is set.")
            else:
                candidate_widget_types = [self.default_widget_type]
        elif n_matches > 1:
            raise MultipleMatchError("Too many candidate types idenfitied ({0}).  Specify enough keywords to guarantee unique type identification.".format(n_matches))

        # Only one is found
        WidgetType = candidate_widget_types[0]

        return WidgetType(data, meta, **kwargs)




LightCurve = LightCurveFactory(default_widget_type=GenericLightCurve,
                 additional_validation_functions=['is_datasource_for'])