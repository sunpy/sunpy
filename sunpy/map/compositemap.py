"""A Composite Map class

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
from __future__ import absolute_import

import sunpy

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

class CompositeMap:
    """Class representing a stack of several Maps"""    
    def __init__(self, *args):
        self._maps = []
        
        # Parse input Maps/filepaths
        zindex = 0
        
        for input_ in args:
            # Parse map
            m = sunpy.Map(input_)
            
            # Give the map a z-index and opacity
            m.zindex = zindex
            m.opacity = 1

            self._maps.append(m)
            
            zindex += 10

    def add_map(self, input_, zindex=None, opacity=1):
        """Adds a map to the CompositeMap"""
        if zindex is None:
            zindex = max([m.zindex for m in self._maps]) + 10
        
        m = sunpy.Map(input_)
        m.zindex = zindex
        m.opacity = opacity
        
        self._maps.append(m)
        
    def remove_map(self, index):
        """Removes and returns the map with the given index"""
        return self._maps.pop(index)
    
    def list_maps(self):
        """Prints a list of the currently included maps"""
        print self._maps

        
