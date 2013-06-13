from __future__ import absolute_import

import os
import glob
import urllib2

import numpy as np

from sunpy.map import GenericMap, MapBase
from sunpy.map.header import MapMeta
#from sunpy.map.mapcube import MapCube
#from sunpy.map.compositemap import CompositeMap
#from sunpy.map.sources import *

from sunpy.io.file_tools import read_file
from sunpy.io.header import FileHeader

from sunpy.util.net import download_file

from sunpy.util.datatype_factory_base import RegisteredFactoryBase

__all__ = ['Map']

class Map(RegisteredFactoryBase):
	
    GenericWidgetType = GenericMap
    

    @classmethod
    def _read_file(cls, fname):
        # File gets read here.  This needs to be generic enough to seamlessly
        #call a fits file or a jpeg2k file, etc
        
        filedata, filemeta  = read_file(fname)
        
        assert isinstance(filemeta, FileHeader)
        
        data = filedata
        meta = MapMeta(filemeta)
        
        return (data, meta)

    @classmethod
    def _parse_args(cls, *args, **kwargs):
        
        data_header_pairs = list()
        already_maps = list()
        
        # For each of the arguments, handle each of the cases
        i = 0
        while i < len(args):
            
            arg = args[i]
            
            # Data-header pair in a tuple
            if ((type(arg) in [tuple, list]) and 
                 isinstance(arg[0],np.ndarray) and # or NDData or something else?
                 isinstance(arg[1],SunpyMetaBase)): # FITSHeader, JP2kHeader, OrderedDict, dict?
                data_header_pairs.append(arg)
            
            # Data-header pair not in a tuple
            elif (isinstance(arg, np.ndarray) and # or NDData or something else?
                  isinstance(args[i+1],SunpyMetaBase)): # FITSHeader, JP2kHeader, OrderedDict, dict? 
                pair = (args[i], args[i+1])
                data_header_pairs.append(pair)
                i += 1 # an extra increment to account for the data-header pairing
            
            # File name
            elif (isinstance(arg,basestring) and 
                  os.path.isfile(os.path.expanduser(arg))):
                path = os.path.expanduser(arg)
                pair = cls._read_file(path)
                data_header_pairs.append(pair)
            
            # Directory
            elif (isinstance(arg,basestring) and 
                  os.path.isdir(os.path.expanduser(arg))):
                path = os.path.expanduser(arg)
                files = [os.path.join(path, elem) for elem in os.listdir(path)]
                data_header_pairs += map(cls._read_files, files)
            
            # Glob
            elif (isinstance(arg,basestring) and '*' in arg):
                files = glob.glob( os.path.expanduser(arg) )
                data_header_pairs += map(cls._read_file, files)
            
            # Already a Map
            elif isinstance(arg, MapBase):
                already_maps.append(arg)
                
            # A URL
            elif (isinstance(arg,basestring) and 
                  urllib2.urlopen(arg)):
                default_dir = sunpy.config.get("downloads", "download_dir")
                url = arg
                path = download_file(url, default_dir)
                pair = cls._read_file(path)
                data_header_pairs.append(pair)
                
            i += 1
        #TODO:
        # In the end, if there are aleady maps it should be put in the same
        # order as the input, currently they are not.
        
        return data_header_pairs, already_maps
    
    
    def __new__(cls, *args, **kwargs):

        # Hack to get around Python 2.x not backporting PEP 3102.
        composite = kwargs.pop('composite', False)
        cube = kwargs.pop('cube', False)

        if cls is Map:
            
            # Get list of data-header pairs, e.g., [(d1, h1), (d2, h2), ...]
            data_header_pairs, already_maps = cls._parse_args(*args, **kwargs)
            
            # If the list is meant to be a cube, instantiate a map cube
            if cube:
                return MapCube(*data_header_pairs, **kwargs)

            # If the list is meant to be a composite mape, instantiate one
            if composite:
                return CompositeMap(*data_header_pairs, **kwargs)
                
            # Otherwise, each pair in the list gets built on its own
            new_maps = list()
            
            for data, header in data_header_pairs:
                # Test to see which type of Map this pair is.  If none of the
                # registered Map types match, use a generic map.
                WidgetType = None
                for key in cls.registry:
                    
                    if cls.registry[key](data, header, **kwargs):
                        WidgetType = key
                        break
                else:
                    WidgetType = cls.GenericWidgetType
                    
                # Instantiate the new map.
                new_maps.append(WidgetType(data, header, **kwargs))
            
            new_maps += already_maps
            
            # If there was only one map instantiated, return that, otherwise
            # return the list of them.
            return new_maps[0] if len(new_maps) == 1 else new_maps
        else:
            return super(Map, cls).__new__(cls, *args, **kwargs)
 


#def make_map(*args, **kwargs):
#    """Processes one or more inputs and returns a Map, MapCube, or CompositeMap
#    instance.
#
#    Parameters
#    ----------
#    args : filepath(s), data array
#        The data source used to create the map object. This can be either a
#        filepath to an image, a 2d list, or an ndarray.
#    type : {'composite' | 'cube'}
#        Type of multimap to construct when passed more than one input. The
#        default choice is a CompositeMap which is more lenient with respect
#        to how similar the input data is.
#
#    Returns
#    -------
#    out : Map, MapCube, CompositeMap
#        Returns a  subclass instance
#
#    Examples
#    --------
#    >>> import sunpy
#    >>> sunpy.make_map("file.fts")
#    >>> sunpy.make_map("file1.fts", "file2.fts",..)
#    >>> sunpy.make_map(["file1.fts", "file2.fts",..])
#    >>> sunpy.make_map("path/to/files/*.fts")
#    >>> sunpy.make_map(Map)
#    >>> sunpy.make_map(Map1, Map2,..)
#    >>> sunpy.make_map([[0, 1],[2, 3]], {'telescop': 'sunpy',..})
#
#    """
#    if len(args) is 0:
#        raise TypeError("Invalid input.")
#
#    # First check to see if data/header were passed in    
#    if isinstance(args[0], list) or isinstance(args[0], np.ndarray):
#        data = None
#
#        # n-dimensional list
#        if isinstance(args[0][0], list) or isinstance(args[0], np.ndarray):
#            data = args[0]
#        else:
#            try:
#                float(args[0][0])
#            except (ValueError, TypeError):
#                pass
#            else:
#                # 1-dimensional data
#                data = args[0]
#
#        # if either of the above cases hold, then create a new Map
#        if data is not None:
#            if len(args) > 1:
#                return Map(args[0], args[1])
#            else:
#                return Map(args[0], {})
#
#
#    # If not, check for one or more maps or filepaths
#    if len(args) == 1:
#        # String
#        if isinstance(args[0], basestring):
#            filepath = os.path.expanduser(args[0])
#
#            # Wildcard string
#            if filepath.find("*") != -1:
#                maps = glob.glob(filepath)
#            # Directory (use all files)
#            elif os.path.isdir(filepath):
#                maps = [os.path.join(filepath, x) for x in os.listdir(filepath)]
#
#            # Filepath
#            else:
#                return Map.read(filepath)
#
#        # Map/MapCube/CompositeMap
#        elif (isinstance(args[0], Map) or 
#              isinstance(args[0], CompositeMap) or 
#              isinstance(args[0], MapCube)):
#            return args[0]
#
#        # List of filepaths or Maps
#        elif isinstance(args[0], list):
#            # list of maps or filepaths
#            maps = args[0]
#
#        # Unrecognized input
#        else:
#            raise InvalidMapInput("Invalid input for make_map. Please specify "
#                                  "one or more filepaths, Maps, directories, "
#                                  "or wildcard expressions.")
#    else:
#        maps = args
#
#    # Make sure we found some data
#    if len(maps) is 0:
#        raise NoMapsFound("Specified path contains no valid files.")
#
#    mtype = kwargs.get("type", "composite")
#
#    # MapCube
#    if mtype == "cube":
#        return MapCube(*maps)
#    # CompositeMap (default)
#    elif mtype == "composite":
#        return CompositeMap(*maps)
#    else:
#        raise InvalidMapType("Invalid multi-map type specified. Please choose "
#                             "between 'composite' or 'cube'.")

#def read_header(filepath):
#    """Parses a file header and return some important parameters"""
#    return Map.read_header(filepath)

class InvalidMapInput(ValueError):
    """Exception to raise when input variable is not a Map instance and does
    not point to a valid Map input file."""
    pass

class InvalidMapType(ValueError):
    """Exception to raise when an invalid type of map is requested with make_map
    """
    pass

class NoMapsFound(ValueError):
    """Exception to raise when input does not point to any valid maps or files 
    """
    pass