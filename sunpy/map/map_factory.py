from __future__ import absolute_import

__authors__ = ["Russell Hewett, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

import os
import glob
import urllib2

from collections import OrderedDict

import numpy as np

from sunpy.map import GenericMap
from sunpy.map.header import MapMeta
from sunpy.map.compositemap import CompositeMap
from sunpy.map.mapcube import MapCube

from sunpy.io.file_tools import read_file
from sunpy.io.header import FileHeader

from sunpy.util.net import download_file
from sunpy.util import expand_list
from sunpy.util.datatype_factory_base import RegisteredFactoryBase
from sunpy.util import Deprecated

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
        
        # Account for nested lists of items
        args = expand_list(args)
        
        # For each of the arguments, handle each of the cases
        i = 0
        while i < len(args):
            
            arg = args[i]
            
            # Data-header pair in a tuple
            if ((type(arg) in [tuple, list]) and 
                 isinstance(arg[0],np.ndarray) and # or NDData or something else?
                 isinstance(arg[1],OrderedDict)): # FITSHeader, JP2kHeader, OrderedDict, dict?
                data_header_pairs.append(arg)
            
            # Data-header pair not in a tuple
            elif (isinstance(arg, np.ndarray) and # or NDData or something else?
                  isinstance(args[i+1],OrderedDict)): # FITSHeader, JP2kHeader, OrderedDict, dict? 
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
                data_header_pairs += map(cls._read_file, files)
            
            # Glob
            elif (isinstance(arg,basestring) and '*' in arg):
                files = glob.glob( os.path.expanduser(arg) )
                data_header_pairs += map(cls._read_file, files)
            
            # Already a Map
            elif isinstance(arg, GenericMap):
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
            
            # If the list is meant to be a cube, instantiate a map cube
            if cube:
                return MapCube(new_maps, **kwargs)

            # If the list is meant to be a composite mape, instantiate one
            if composite:
                return CompositeMap(new_maps, **kwargs)
            
            # If there was only one map instantiated, return that, otherwise
            # return the list of them.
            return new_maps[0] if len(new_maps) == 1 else new_maps
        else:
            return super(Map, cls).__new__(cls, *args, **kwargs)
 

@Deprecated("Please use the new factory sunpy.Map")
def make_map(*args, **kwargs):
    __doc__ = Map.__doc__
    Map(*args, **kwargs)
    


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