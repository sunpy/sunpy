# -*- coding: utf-8 -*-
"""
Provides processing routines for data captured with the AIA instrument on SDO.
"""
from sunpy.map.sources.sdo import AIAMap

def aiaprep(aiamap):
    """    
    Processes a level 1 AIAMap into a level 1.5 AIAMap

    Parameters
    ----------
    aiamap: AIAMap instance
        A sunpy map from AIA

    Returns
    -------
    A level 1.5 copy of aiamap
    """
    assert isinstance(aiamap, AIAMap)

    scale_factor = aiamap.scale['x'] / 0.6

    newmap = aiamap.rotate(recenter=True, scale=scale_factor,
                           missing=aiamap.min())

    # Update header values as needed
    newmap.meta['lvl_num'] = 1.5

    return newmap
