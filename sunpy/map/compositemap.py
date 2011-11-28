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
        for input_ in args:
            self._maps.append(sunpy.Map(input_))
