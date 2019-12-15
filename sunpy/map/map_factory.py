import os
from collections import OrderedDict
import warnings
import pathlib
import glob

import numpy as np
import astropy.io.fits
from astropy.wcs import WCS

import sunpy
from sunpy import log
from sunpy.map.mapbase import GenericMap, MapMetaValidationError
from sunpy.map.compositemap import CompositeMap
from sunpy.map.mapsequence import MapSequence

from sunpy.data import cache

from sunpy.io.file_tools import read_file
from sunpy.io.header import FileHeader

from sunpy.util import expand_list
from sunpy.util.metadata import MetaDict
from sunpy.util.exceptions import SunpyDeprecationWarning
from sunpy.util.types import DatabaseEntryType

from sunpy.util.datatype_factory_base import BasicRegistrationFactory
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.datatype_factory_base import MultipleMatchError
from sunpy.util.datatype_factory_base import ValidationFunctionError
from urllib.request import urlopen

SUPPORTED_ARRAY_TYPES = (np.ndarray,)
try:
    import dask.array
    SUPPORTED_ARRAY_TYPES += (dask.array.Array,)
except ImportError:
    pass

__authors__ = ["Russell Hewett, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

__all__ = ['Map', 'MapFactory']


class MapFactory(BasicRegistrationFactory):
    """
    Map(\\*args, \\*\\*kwargs)

    Map factory class.  Used to create a variety of Map objects.  Valid map types
    are specified by registering them with the factory.


    Examples
    --------
    >>> import sunpy.map
    >>> from astropy.io import fits
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> mymap = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA

    The SunPy Map factory accepts a wide variety of inputs for creating maps

    * Preloaded tuples of (data, header) pairs

    >>> mymap = sunpy.map.Map((data, header))   # doctest: +SKIP

    headers are some base of `dict` or `collections.OrderedDict`, including
    `sunpy.io.header.FileHeader` or `sunpy.util.metadata.MetaDict` classes.

    * data, header pairs, not in tuples

    >>> mymap = sunpy.map.Map(data, header)   # doctest: +SKIP

    * data, wcs object, in tuple

    >>> from astropy.wcs import WCS
    >>> wcs = WCS(sunpy.data.sample.AIA_171_ROLL_IMAGE)     # doctest: +REMOTE_DATA
    >>> data = fits.getdata(sunpy.data.sample.AIA_171_ROLL_IMAGE)    # doctest: +REMOTE_DATA
    >>> mymap = sunpy.map.Map((data, wcs))    # doctest: +REMOTE_DATA

    * data, wcs object, not in tuple

    >>> from astropy.wcs import WCS
    >>> wcs = WCS(sunpy.data.sample.AIA_171_ROLL_IMAGE)     # doctest: +REMOTE_DATA
    >>> data = fits.getdata(sunpy.data.sample.AIA_171_ROLL_IMAGE)    # doctest: +REMOTE_DATA
    >>> mymap = sunpy.map.Map(data, wcs)   # doctest: +REMOTE_DATA

    * File names

    >>> mymap = sunpy.map.Map('file1.fits')   # doctest: +SKIP

    * All fits files in a directory by giving a directory

    >>> mymap = sunpy.map.Map('local_dir/sub_dir')   # doctest: +SKIP

    * A filesystem path expressed as a `pathlib.Path`

    >>> import pathlib
    >>> mymap = sunpy.map.Map(pathlib.Path('file1.fits'))   # doctest: +SKIP
    >>> sub_dir = pathlib.Path('local_dir/sub_dir')
    >>> mymap = sunpy.map.Map(sub_dir)   # doctest: +SKIP
    >>> mymap = sunpy.map.Map(sub_dir / 'file3.fits')   # doctest: +SKIP

    * Some regex globs

    >>> mymap = sunpy.map.Map('eit_*.fits')   # doctest: +SKIP

    * URLs

    >>> mymap = sunpy.map.Map(url_str)   # doctest: +SKIP

    * DatabaseEntry

    >>> mymap = sunpy.map.Map(db_result)   # doctest: +SKIP

    * Lists of any of the above

    >>> mymap = sunpy.map.Map(['file1.fits', 'file2.fits', 'file3.fits', 'directory1/'])  # doctest: +SKIP

    * Any mixture of the above not in a list

    >>> mymap = sunpy.map.Map(((data, header), data2, header2, 'file1.fits', url_str, 'eit_*.fits'))  # doctest: +SKIP
    """  # noqa

    def _read_file(self, fname, **kwargs):
        """ Read in a file name and return the list of (data, meta) pairs in
            that file. """

        # File gets read here.  This needs to be generic enough to seamlessly
        # call a fits file or a jpeg2k file, etc
        # NOTE: use os.fspath so that fname can be either a str or pathlib.Path
        # This can be removed once read_file supports pathlib.Path
        log.debug(f'Reading {fname}')
        pairs = read_file(os.fspath(fname), **kwargs)

        new_pairs = []
        for pair in pairs:
            filedata, filemeta = pair
            assert isinstance(filemeta, FileHeader)
            # This tests that the data is more than 1D
            if len(np.shape(filedata)) > 1:
                data = filedata
                meta = MetaDict(filemeta)
                new_pairs.append((data, meta))
        return new_pairs

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

    def _parse_args(self, *args, **kwargs):
        """
        Parses an args list for data-header pairs.  args can contain any
        mixture of the following entries:
        * tuples of data,header
        * data, header not in a tuple
        * data, wcs object in a tuple
        * data, wcs object not in a tuple
        * filename, as a str or pathlib.Path, which will be read
        * directory, as a str or pathlib.Path, from which all files will be read
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

        data_header_pairs = list()
        already_maps = list()

        # Account for nested lists of items
        args = expand_list(args)

        # For each of the arguments, handle each of the cases
        i = 0
        while i < len(args):

            arg = args[i]

            # Data-header or data-WCS pair
            if isinstance(arg, SUPPORTED_ARRAY_TYPES):
                arg_header = args[i+1]
                if isinstance(arg_header, WCS):
                    arg_header = args[i+1].to_header()

                if self._validate_meta(arg_header):
                    pair = (args[i], OrderedDict(arg_header))
                    data_header_pairs.append(pair)
                    i += 1    # an extra increment to account for the data-header pairing

            # A database Entry
            elif isinstance(arg, DatabaseEntryType):
                data_header_pairs += self._read_file(arg.path, **kwargs)

            # Already a Map
            elif isinstance(arg, GenericMap):
                already_maps.append(arg)

            # URL
            elif isinstance(arg, str) and _is_url(arg):
                url = arg
                path = str(cache.download(url).absolute())
                pairs = self._read_file(path, **kwargs)
                data_header_pairs += pairs

            # File system path (file or directory or glob)
            elif _possibly_a_path(arg):
                path = pathlib.Path(arg).expanduser()
                if _is_file(path):
                    pairs = self._read_file(path, **kwargs)
                    data_header_pairs += pairs
                elif _is_dir(path):
                    for afile in sorted(path.glob('*')):
                        data_header_pairs += self._read_file(afile, **kwargs)
                elif glob.glob(os.path.expanduser(arg)):
                    for afile in sorted(glob.glob(os.path.expanduser(arg))):
                        data_header_pairs += self._read_file(afile, **kwargs)

                else:
                    raise ValueError(f'Did not find any files at {arg}')

            else:
                raise ValueError(f"Invalid input: {arg}")

            i += 1

        # TODO:
        # In the end, if there are already maps it should be put in the same
        # order as the input, currently they are not.
        return data_header_pairs, already_maps

    def __call__(self, *args, composite=False, sequence=False, silence_errors=False, **kwargs):
        """ Method for running the factory. Takes arbitrary arguments and
        keyword arguments and passes them to a sequence of pre-registered types
        to determine which is the correct Map-type to build.

        Arguments args and kwargs are passed through to the validation
        function and to the constructor for the final type.  For Map types,
        validation function must take a data-header pair as an argument.

        Parameters
        ----------
        composite : `bool`, optional
            Indicates if collection of maps should be returned as a `~sunpy.map.CompositeMap`.
            Default is `False`.

        sequence : `bool`, optional
            Indicates if collection of maps should be returned as a `sunpy.map.MapSequence`.
            Default is `False`.

        silence_errors : `bool`, optional
            If set, ignore data-header pairs which cause an exception.
            Default is ``False``.

        Notes
        -----
        Extra keyword arguments are passed through to `sunpy.io.read_file` such
        as `memmap` for FITS files.
        """

        data_header_pairs, already_maps = self._parse_args(*args, **kwargs)

        new_maps = list()

        # Loop over each registered type and check to see if WidgetType
        # matches the arguments.  If it does, use that type.
        for pair in data_header_pairs:
            data, header = pair
            meta = MetaDict(header)

            try:
                new_map = self._check_registered_widgets(data, meta, **kwargs)
                new_maps.append(new_map)
            except (NoMatchError, MultipleMatchError, ValidationFunctionError):
                if not silence_errors:
                    raise
            except MapMetaValidationError as e:
                warnings.warn(f"One of the data, header pairs failed to validate with: {e}")
            except Exception:
                raise

        new_maps += already_maps

        if not len(new_maps):
            raise RuntimeError('No maps loaded')

        # If the list is meant to be a sequence, instantiate a map sequence
        if sequence:
            return MapSequence(new_maps, **kwargs)

        # If the list is meant to be a composite map, instantiate one
        if composite:
            return CompositeMap(new_maps, **kwargs)

        if len(new_maps) == 1:
            return new_maps[0]

        return new_maps

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
            raise MultipleMatchError("Too many candidate types identified ({})."
                                     "Specify enough keywords to guarantee unique type"
                                     "identification.".format(n_matches))

        # Only one is found
        WidgetType = candidate_widget_types[0]

        return WidgetType(data, meta, **kwargs)


def _is_url(arg):
    try:
        urlopen(arg)
    except Exception:
        return False
    return True


def _possibly_a_path(arg):
    """
    Check if arg can be coerced into a Path object.
    Does *not* check if the path exists.
    """
    try:
        is_path = pathlib.Path(arg)
        return True
    except Exception:
        return False


# In python<3.8 paths with un-representable chars (ie. '*' on windows)
# raise an error, so make our own version that returns False instead of
# erroring. These can be removed when we support python >= 3.8
# https://docs.python.org/3/library/pathlib.html#methods
def _is_file(path):
    try:
        return path.is_file()
    except Exception:
        return False


def _is_dir(path):
    try:
        return path.is_dir()
    except Exception:
        return False


class InvalidMapInput(ValueError):
    """Exception to raise when input variable is not a Map instance and does
    not point to a valid Map input file."""
    pass


class InvalidMapType(ValueError):
    """Exception to raise when an invalid type of map is requested with Map
    """
    pass


class NoMapsFound(ValueError):
    """Exception to raise when input does not point to any valid maps or files
    """
    pass


Map = MapFactory(registry=GenericMap._registry, default_widget_type=GenericMap,
                 additional_validation_functions=['is_datasource_for'])
