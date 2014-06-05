# -*- coding: utf-8 -*-
"""
Provides processing routines for data captured with the AIA instrument on SDO.
"""
from sunpy.map.sources.sdo import AIAMap

def aiaprep(aiamap):
    """    
    Processes a level 1 AIAMap into a level 1.5 AIAMap. Rotates, scales and
    translates the image so that solar North is aligned with the y axis, each
    pixel is 0.6 arcsec across, and the centre of the sun is at the centre of
    the image. The actual transformation is done by :func:`Map.rotate`.
    
    This function is equivalent in functionality to aia_prep() in SSWIDL, but
    it does not use the same transformation to rotate the image and should
    therefore not be expected to produce the same results.

    Parameters
    ----------
    aiamap : AIAMap instance
        A sunpy map from AIA

    Returns
    -------
    newmap : A level 1.5 copy of aiamap
    """
    assert isinstance(aiamap, AIAMap)
    scale_factor = aiamap.scale['x'] / 0.6
    newmap = aiamap.rotate(recenter=True, scale=scale_factor, missing=aiamap.min())
    newmap.meta['lvl_num'] = 1.5

    return newmap
