# -*- coding: utf-8 -*-
"""
Provides processing routines for data captured with the AIA instrument on SDO.
"""
from sunpy.map.sources.sdo import AIAMap

def aiaprep(aiamap):
    """    
    Processes a level 1 AIAMap into a level 1.5 AIAMap. Rotates, scales and
    translates the image so that solar North is aligned with the y axis, each
    pixel is 0.6 arcsec across, and the center of the sun is at the center of
    the image. The actual transformation is done by Map's
    :meth:`~sunpy.map.mapbase.GenericMap.rotate` method.
    
    This function is similar in functionality to aia_prep() in SSWIDL, but
    it does not use the same transformation to rotate the image and it handles
    the meta data differently. It should therefore not be expected to produce
    the same results.

    Parameters
    ----------
    aiamap : AIAMap instance
        A sunpy map from AIA

    Returns
    -------
    newmap : A level 1.5 copy of aiamap
    
    Notes
    -----
    This routine makes use of Map's :meth:`~sunpy.map.mapbase.GenericMap.rotate`
    method, which modifes the header information to the standard PCi_j WCS
    formalism.
    The FITS header resulting in saving a file after this procedure will
    therefore differ from the original file.
    """
    
    if not isinstance(aiamap, AIAMap):
        raise ValueError("Input must be an AIAMap")

    # Taget scale is 0.6 arcsec/pixel, but this needs to be adjusted if the map
    # has already been rescaled.
    if round(aiamap.scale['x']/0.6) != 1.0 and aiamap.shape != (4096, 4096):
        scale = round(aiamap.scale['x']/0.6) * 0.6
    else:
        scale = 0.6 # pragma: no cover # can't test this because it needs a full res image
    scale_factor = aiamap.scale['x'] / scale
    
    newmap = aiamap.rotate(recenter=True, scale=scale_factor, missing=aiamap.min())
    newmap.meta['lvl_num'] = 1.5

    return newmap
