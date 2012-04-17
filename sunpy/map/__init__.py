"""SunPy Maps"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import os
from sunpy.map.header import MapHeader
from sunpy.map.basemap import BaseMap
from sunpy.map.mapcube import MapCube
from sunpy.map.compositemap import CompositeMap
from sunpy.map.sources import *

def make_map(*args, **kwargs):
    """Processes one or more inputs and returns a Map, MapCube, or CompositeMap
    instance.
    
    Parameters
    ----------
    args : filepath(s), data array
        The data source used to create the map object. This can be either a
        filepath to an image, a 2d list, or an ndarray.
    type : {'composite' | 'cube'}
        Type of multimap to construct when passed more than one input. The
        default choice is a CompositeMap which is more lenient with respect
        to how similar the input data is.
        
    Returns
    -------
    out : Map, MapCube, CompositeMap
        Returns a  subclass instance
        
    Examples
    --------
    >>> import sunpy
    >>> sunpy.make_map("file.fts")
    >>> sunpy.make_map("file1.fts", "file2.fts",..)
    >>> sunpy.make_map(["file1.fts", "file2.fts",..])
    >>> sunpy.make_map("path/to/files/*.fts")
    >>> sunpy.make_map(Map)
    >>> sunpy.make_map(Map1, Map2,..)

    """
    # Single Map or wildcard string
    if len(args) == 1:
        # String
        if isinstance(args[0], basestring):
            filepath = os.path.expanduser(args[0])
            
            # Wildcard string
            if filepath.find("*") != -1:
                import glob
                maps = glob.glob(filepath)
            # Directory (use all files)
            elif os.path.isdir(filepath):
                maps = [os.path.join(filepath, x) for x in os.listdir(filepath)]
                
            # Filepath
            else:
                return BaseMap.read(filepath)

        # Map/MapCube/CompositeMap
        elif (isinstance(args[0], BaseMap) or 
              isinstance(args[0], CompositeMap) or 
              isinstance(args[0], MapCube)):
            return args[0]
        
        # List of filepaths or Maps
        elif isinstance(args[0], list):
            maps = args[0]

        # Unrecognized input
        else:
            raise InvalidMapInput("Invalid input for make_map. Please specify "
                                  "one or more filepaths, Maps, directories, "
                                  "or wildcard expressions.")
    else:
        maps = args
        
    # Make sure we found some data
    if len(maps) is 0:
        raise NoMapsFound("Specified path contains no valid files.")
        
    mtype = kwargs.get("type", "composite")
        
    # MapCube
    if mtype == "cube":
        return MapCube(*maps)
    # CompositeMap (default)
    elif mtype == "composite":
        return CompositeMap(*maps)
    else:
        raise InvalidMapType("Invalid multi-map type specified. Please choose "
                             "between 'composite' or 'cube'.")
        
def read_header(filepath):
    """Parses a file header and return some important parameters"""
    return BaseMap.read_header(filepath)
    
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
