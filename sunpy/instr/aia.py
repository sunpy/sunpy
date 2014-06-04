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

    scale_ref = 0.6
    scale_factor = aiamap.scale['x'] / scale_ref

    newmap = aiamap.rotate(recenter=True, scale=scale_factor,
                           missing=aiamap.min())

    # Update header values as needed
    newmap.meta['cdelt1'] = scale_ref
    newmap.meta['cdelt2'] = scale_ref
    newmap.meta['lvl_num'] = 1.5

    return newmap
