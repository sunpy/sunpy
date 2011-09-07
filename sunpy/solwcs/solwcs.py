from __future__ import absolute_import

"""
Solar WCS provides functions to parse a World Coordinate System (WCS) fits header.

Note that SOLAR_B0, HGLT_OBS, and CRLT_OBS are all synonyms.
"""

__authors__ = ["Steven Christe"]
__email__ = "steven.d.christe@nasa.gov"

import numpy as np
from datetime import datetime
from sunpy.sun import constants

import pyfits
import sunpy

def solar_limb(header):
    """Return the angular size of the Sun viewed from Earth (in arcsec)"""
    
    # check to see if the answer is provided in the fits header itself
    rsun_obs = header.get('RSUN_OBS')
    
    # As a last resort just set it to a default value
    # TODO: don't hardcode this number
    if rsun_obs is None: rsun_obs = 960.0
    
    return rsun_obs

def observer_position(header):
    """Return the observer distance from the Sun."""
    
    # check to see if the answer is provided in the fits header itself
    dsun_obs = header.get('DSUN_OBS')
    
    return dsun_obs

def get_center(header, axis = None):
    """Return the center of the map."""
    x = header.get('cdelt1')*header.get('naxis1')/2 + header.get('crval1') - header.get('crpix1')*header.get('cdelt1')
    y = header.get('cdelt2')*header.get('naxis2')/2 + header.get('crval2') - header.get('crpix2')*header.get('cdelt2')
    
    if axis is 'x': return x
    if axis is 'y': return y
    if axis is None: return [x,y]
    
def get_units(header, axis = None):
    """Return the units used for crpix, crdelt in the header."""
    xunits = header.get('cunit1')
    yunits = header.get('cunit2')
    
    if xunits is None: xunits = header.get('ctype1')
    if yunits is None: yunits = header.get('ctype2')
    
    if axis is 'x': return xunits
    if axis is 'y': return yunits
    if axis is None: return [xunits,yunits]
    
def get_platescale(header, axis = None):
    """Return the plate scale of the image, i.e. the size of the pixels in unit."""
    xscale = header.get('cdelt1')
    yscale = header.get('cdelt2')
     
    if axis is 'x': return xscale
    if axis is 'y': return yscale
    if axis is None: return [xscale,yscale]
    
def get_solar_b0(header):
    """Return the solar B0 angle which is simply the heliographic latitude of the observer."""
    
    # check to see if the answer is provided by the header
    solar_b0 = header.get('HGLT_OBS')
    
    # check in another place
    if solar_b0 is None: solar_b0 = header.get('CRLT_OBS')
    
    # check in another place
    if solar_b0 is None: solar_b0 = header.get('SOLAR_B0')
    
    # As a last resort just set it to a default value
    if solar_b0 is None: solar_b0 = 0
    
    return solar_b0

def get_carrington_hg_longitude(header):
    """Return the Carrington heliographic longitude of the observer."""
    
    # check to see if the answer is provided by the header
    crln_obs = header.get('CRLN_OBS')

    # TODO: as a last resort, assume the observer is at Earth
    crln_obs = 0
    
    return crln_obs
    
def convert_ang_units(type = 'hpc', unit = 'arcsec'):
    '''Determine the conversion factor between the data and radians.'''
    
    if unit == 'arcmin': return np.deg2rad(1)/60.0
    if unit == 'arcsec': return np.deg2rad(1)/(60*60.0)
    if unit == 'mas': return np.deg2rad(1)/(60*60*1000.0)

def get_projection(header, axis = 'x'):
    """Return the projection that the data was taken in."""
    xtype = header.get('ctype1')
    ytype = header.get('ctype2')
     
    if axis is 'x': return xtype
    if axis is 'y': return ytype
    # Assume that the projection is the same in both axis
    # TODO: Remove assumption of same projection in both axis
    if axis is None: return xtype

def get_shape(header):
    """Return the shape of the data array."""
    return [header.get('naxis1'), header.get('naxis2')]

def convert_data_to_coord(header, pixel_index = None):
    """This procedure takes a WCS-compliant header, and calculates the data coordinates at each pixel position."""
    if pixel_index is not None: pixel_index = np.array(pixel_index)
    naxis = get_shape(header)
    cdelt = get_platescale(header)
    crpix = [header.get('crpix1'), header.get('crpix2')]
    
    # first assume that coord is just [x,y]
    coord = np.zeros(pixel_index.shape)
    
    coord[0] = pixel_index[0] - (crpix[0] - 1)
    coord[1] = pixel_index[1] - (crpix[1] - 1)

    coord[0] = coord[0] * cdelt[0]
    coord[1] = coord[1] * cdelt[1]
        
    if pixel_index is None:
        
        xcoord = np.zeros(naxis)
        ycoord = np.zeros(naxis)
            
    # check to see what projection is being used
    projection = get_projection(header)
    if  projection.count('TAN'):    
        coord = proj_tan(header, coord)
        
    return coord

def test(x, y, hpln, hplt):
    
    fits = pyfits.open(sunpy.AIA_171_IMAGE)
    header = fits[0].header
    
    r = convert_data_to_coord(header,pixel_index = [5,0])
    print(r)
    
    r = convert_hpc_hcc(header, x,y)
    print(r)
    r = convert_hcc_hpc(header, hpln, hplt)
    print(r)
    
def convert_hpc_hcc(header, hpln, hplt, distance = None):
    """This routine converts Helioprojective-Cartesian (HPC) coordinates into Heliocentric-Cartesian (HCC) coordinates, using equations 15 in Thompson (2006), A&A, 449, 791-803."""
    
    cx = convert_ang_units(unit = get_units(header, axis = 'x'))
    cy = convert_ang_units(unit = get_units(header, axis = 'y'))
    
    lon = cx * hpln
    lat = cy * hplt

    # Calculate the sines and cosines.
    cosx = np.cos(lon)
    sinx = np.sin(lon)
    cosy = np.cos(lat)
    siny = np.sin(lat)
    
    dsun = header.get('dsun_obs')
    # Should we use the rsun_ref defined in the fits file or our local (possibly different/more correct) definition
    rsun = header.get('rsun_ref')
    
    if distance is None: 
        q = dsun * cosy * cosx
        distance = q ** 2 - dsun ** 2 + rsun ** 2
        # Need to check if there are values which are negative and get rid of them
        # IDL code
        # w = where(distance lt 0, count)
        # if count gt 0 then flag_missing, distance, 
        distance = q - np.sqrt(distance) 
    
    x = distance * cosy * sinx
    y = distance * siny
    z = dsun - distance * cosy * cosx
    return [x,y,z]

def convert_hcc_hpc(header, x, y, distance = None):
    """Convert Heliocentric-Cartesian (HCC) to angular Helioprojective-Cartesian (HPC) coordinates (in degrees)."""
    
    dsun = header.get('dsun_obs')
    # Should we use the rsun_ref defined in the fits file or our local (possibly different/more correct) definition
    rsun = header.get('rsun_ref')
    
    # Calculate the z coordinate by assuming that it is on the surface of the Sun
    z = rsun ** 2 - x ** 2 - y ** 2
    z = np.sqrt( z )
    
    zeta = dsun - z
    distance = np.sqrt(x ** 2 + y ** 2 + zeta ** 2)
    hpln = np.arctan(x / zeta)
    hplt = np.arcsin(y / distance)
    
    # convert the results to degrees
    result = np.rad2deg([hpln, hplt])
    return result

def proj_tan(header, coord, force = False):
    """Applies the gnomonic (TAN) projection to intermediate relative coordinates."""
    
    # if pixels are within 3 degrees of the Sun than skip the calculatin unless force is True
    # This applies to all sdo images so this function is just here as a place holder
    return coord
    