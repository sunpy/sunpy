from __future__ import absolute_import

__authors__ = ["Russell Hewett, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

import os
import glob
import urllib2

from sunpy.util.odict import OrderedDict

import numpy as np

import sunpy
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
    """
    Map factory class.  Used to create a variety of Map objects.  Valid map types
    are specified by registering them with the factory.
    
    Example
    -------
    >>> import sunpy
    >>> mymap = sunpy.Map(sunpy.AIA_171_IMAGE)
    
    """

    
    DefaultWidgetType = GenericMap

    @classmethod
    def _read_file(cls, fname, **kwargs):
        """ Read in a file name and return the list of (data, meta) pairs in
            that file. """
        
        # File gets read here.  This needs to be generic enough to seamlessly
        #call a fits file or a jpeg2k file, etc
        pairs = read_file(fname, **kwargs)
        
        new_pairs = []
        for pair in pairs:
            filedata, filemeta = pair
            assert isinstance(filemeta, FileHeader)
            #This tests that the data is more than 1D
            if len(np.shape(filedata)) > 1:
                data = filedata
                meta = MapMeta(filemeta)
                new_pairs.append((data, meta))
        return new_pairs

    @classmethod
    def _parse_args(cls, *args, **kwargs):
        """
        Parses an args list for data-header pairs.  args can contain any mixture
        of the following entries:
        * tuples of data,header
        * data, header not in a tuple
        * filename, which will be read
        * directory, from which all files will be read
        * glob, from which all files will be read
        * url, which will be downloaded and read
        * lists containing any of the above.
        
        Example
        -------
        cls._parse_args(data, header, (data, header), ['file1', 'file2', 'file3'], 'file4', 'directory1', '*.fits')
        
        """ 
        
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
            elif (isinstance(arg,basestring) and 
                  os.path.isfile(os.path.expanduser(arg))):
                path = os.path.expanduser(arg)
                pairs = cls._read_file(path, **kwargs)
                data_header_pairs += pairs
            
            # Directory
            elif (isinstance(arg,basestring) and 
                  os.path.isdir(os.path.expanduser(arg))):
                path = os.path.expanduser(arg)
                files = [os.path.join(path, elem) for elem in os.listdir(path)]
                for afile in files:
                    data_header_pairs += cls._read_file(afile, **kwargs)
            
            # Glob
            elif (isinstance(arg,basestring) and '*' in arg):
                files = glob.glob( os.path.expanduser(arg) )
                for afile in files:
                    data_header_pairs += cls._read_file(afile, **kwargs)
            
            # Already a Map
            elif isinstance(arg, GenericMap):
                already_maps.append(arg)
                
            # A URL
            elif (isinstance(arg,basestring) and 
                  _is_url(arg)):
                default_dir = sunpy.config.get("downloads", "download_dir")
                url = arg
                path = download_file(url, default_dir)
                pairs = cls._read_file(path, **kwargs)
                data_header_pairs += pairs
            
            else:
                raise ValueError("File not found or invalid input")
            
            i += 1
        #TODO:
        # In the end, if there are aleady maps it should be put in the same
        # order as the input, currently they are not.
        return data_header_pairs, already_maps
    
    
    def __new__(cls, *args, **kwargs):
        """
        Factory-construction method.  Takes arbitrary arguments and keyword 
        arguments and passes them to a sequence of pre-registered types to
        determine which is the correct Map-type to build.
        
        Map-type registration must take a data-header pair as an argument.
        
        Parameters
        ----------
        composite: boolean, optional
            Indicates if collection of maps should be returned as a CompositeMap
        cube: boolean, optional
            Indicates if collection of maps should be returned as a MapCube
        
        """

        # Hack to get around Python 2.x not backporting PEP 3102.
        composite = kwargs.pop('composite', False)
        cube = kwargs.pop('cube', False)
        
        if cls is Map:
            
            # Get list of data-header pairs, e.g., [(d1, h1), (d2, h2), ...]
            data_header_pairs, already_maps = cls._parse_args(*args, **kwargs)
                            
            new_maps = list()
            
            # Loop over each registered type and check to see if WidgetType
            # matches the arguments.  If it does, use that type.
            for pair in data_header_pairs:
                data, header = pair
                # Test to see which type of Map this pair is.  If none of the
                # registered Map types match, use a generic map.
                WidgetType = None
                for key in cls.registry:
                    
                    if cls.registry[key](data, header, **kwargs):
                        WidgetType = key
                        break
                else:
                    WidgetType = cls.DefaultWidgetType
                
                # Make sure the map header is a MapMeta, useful if the user
                # passed in a dict
                meta = MapMeta(header)
                # Instantiate the new map.
                new_maps.append(WidgetType(data, meta, **kwargs))
            
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

def _is_url(arg):
    try:
        urllib2.urlopen(arg)
    except:
        return False
    return True

@Deprecated("Please use the new factory sunpy.Map")
def make_map(*args, **kwargs):
    __doc__ = Map.__doc__
    return Map(*args, **kwargs)

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
