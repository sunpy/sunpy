from __future__ import absolute_import

"""
Solar WCS provides functions to parse a World Coordinate System (WCS) fits header.
"""

__authors__ = ["Steven Christe"]
__email__ = "steven.d.christe@nasa.gov"

import numpy as np
from datetime import datetime
from sunpy.sun import constants

def solar_limb(header):
    """Return the angular size of the Sun viewed from Earth (in arcsec)"""
    
    # check to see if the answer is provided in the fits header itself
    rsun_obs = header.get('RSUN_OBS')
    
    if rsun_obs is None: rsun_obs = 960.0
    
    return rsun_obs

def observer_position(header):
    """Return the observer distance from the Sun."""
    
    # check to see if the answer is provided in the fits header itself
    dsun_obs = header.get('DSUN_OBS')
    
    return dsun_obs

def get_center(header, axis = None):
    
    x = header.get('cdelt1')*header.get('naxis1')/2 + header.get('crval1') - header.get('crpix1')*header.get('cdelt1')
    y = header.get('cdelt2')*header.get('naxis2')/2 + header.get('crval2') - header.get('crpix2')*header.get('cdelt2')
    
    if axis is 'x': return x
    if axis is 'y': return y
    if axis is None: return [x,y]
    
def get_units(header, axis = None):
    
    xunits = header.get('cunit1')
    yunits = header.get('cunit2')
    
    if xunits is None: xunits = header.get('ctype1')
    if yunits is None: yunits = header.get('ctype2')
    
    if axis is 'x': return xunits
    if axis is 'y': return yunits
    if axis is None: return [xunits,yunits]