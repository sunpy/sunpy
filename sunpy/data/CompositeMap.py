"""A Composite Map class

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.data.MapCube import MapCube

class CompositeMap(MapCube):
    """
    A subclass of MapCube for dealing with composite or stacked map data.
    
    Parameters
    ----------
    coalign : [ None | 'diff' ] 
    """
    def __init__(self, input_, coalign='diff', **kwargs):
        MapCube.__init__(self, input_, coalign=coalign, **kwargs)
        
    def plot(self):
        """TODO: implement plot method for composite maps"""
        pass