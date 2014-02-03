"""
Provides processing routines for data captured with the AIA instrument on SDO.
"""
import numpy as np

from sunpy.map.sources.sdo import AIAMap

__all__ = ['aiaprep']

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
    assert isinstance(aiamap,AIAMap)
    
    #I believe this is the target pixel scale for SDO data
    scale_ref = 0.6
    #Extract useful parameters from the map
    scale_factor = aiamap.scale['x'] / scale_ref
    x0 = aiamap.reference_pixel['x']
    y0 = aiamap.reference_pixel['y']
    instrot = aiamap.rotation_angle['y'] #This needs to be pulled out into map.
    
    fovx = aiamap.shape[1] * aiamap.scale['x']
    fovy = aiamap.shape[0] * aiamap.scale['y']
    
    #Missing data is repesented as a very negative number
    missing_val = aiamap.min()
    
    #if image has been cropped or resized then do not re-centre image
    if ((fovx < 2456.5) or (fovx > 2458.5) or
        (fovy < 2456.5) or (fovy > 2458.5)):
        
        recenter = False
    else:
        recenter = True

    aiamap = aiamap.rotate(np.radians(-instrot), scale=scale_factor,
                           rotation_center=(y0,x0), recenter=recenter,
                           missing=missing_val, interpolation='bicubic',
                           interp_param=-0.5)
    
        
    
    #; Update header tag values as needed:
    aiamap.meta['crpix1'] = aiamap.shape[1]/2. + 0.5
    aiamap.meta['crpix2'] = aiamap.shape[0]/2. + 0.5
    aiamap.meta['cunit1'] = scale_ref
    aiamap.meta['cunit2'] = scale_ref
    aiamap.meta['crota1'] = 0.0
    aiamap.meta['crota2'] = 0.0
    aiamap.meta['lvl_num']= 1.5
    
    return aiamap
