"""SunPy Maps"""
from __future__ import absolute_import

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.map.mapcube import MapCube
from sunpy.map.compositemap import CompositeMap

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
            # Wildcard string
            if args[0].find("*") != -1:
                import glob
                maps = glob.glob(args[0])
            else:
                # Filepath
                return BaseMap.map_from_filepath(args[0])

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
            raise InvalidMapInput
    else:
        maps = args
        
    mtype = kwargs.get("type", "composite")
        
    # MapCube
    if mtype == "cube":
        return MapCube(*maps)
    # CompositeMap (default)
    elif mtype == "composite":
        return CompositeMap(*maps)
    else:
        raise InvalidMapType
    
class InvalidMapInput(ValueError):
    """Exception to raise when input variable is not a Map instance and does
    not point to a valid Map input file. """
    pass

class InvalidMapType(ValueError):
    """Exception to raise when an invalid type of map is requested with make_map
    """
    pass

if __name__ == "__main__":
    import sunpy
    sunpy.make_map(sunpy.AIA_171_IMAGE, sunpy.RHESSI_IMAGE).show()
