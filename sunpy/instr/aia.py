# -*- coding: utf-8 -*-
"""
Provides processing routines for data captured with the AIA instrument on SDO.
"""
import numpy as np
from sunpy.image.rotate import affine_transform
from sunpy.map.sources.sdo import AIAMap
from copy import deepcopy

def aiaprep(aiamap):
    """Processes a level 1 AIAMap into a level 1.5 AIAMap

    Parameters
    ----------
    aiamap: AIAMap instance
        A sunpy map from AIA

    Returns
    -------
    A level 1.5 copy of aiamap
    """
    assert isinstance(aiamap, AIAMap)

    angle = np.radians(-aiamap.rotation_angle['y'])

    c = np.cos(angle)
    s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])

    rotation_center = aiamap.reference_pixel['y'], aiamap.reference_pixel['x']
    center = (np.array(aiamap.data.shape)-1)/2.0
    shift = np.array(rotation_center) - np.array(center)
    print shift

    newmap = deepcopy(aiamap)
    newmap.data = affine_transform(aiamap.data.copy(), rmatrix=rmatrix, shift=shift)

    # Update header values as needed
    newmap.meta['crpix1'] = newmap.shape[1]/2.0 + 0.5
    newmap.meta['crpix2'] = newmap.shape[0]/2.0 + 0.5
    #
    #
    newmap.meta['crota1'] = 0.0
    newmap.meta['crota2'] = 0.0
    newmap.meta['lvl_num'] = 1.5

    return newmap
