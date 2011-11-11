from __future__ import absolute_import

"""
The WCS package provides functions to parse a World Coordinate System (WCS) fits 
header for solar images as well as convert between various solar coordinate
systems. The solar coordinates supported are

* Helioprojective-Cartesian (HPC): The most often used solar coordinate system.
    Describes positions on the Sun as angles measured from the center of the 
    solar disk (usually in arcseconds) using cartesian coordinates (X, Y)
* Helioprojective-Radial (HPR): Describes positions on the Sun using angles, 
    similar to HPC, but uses a radial coordinate (rho, psi) system centered on
    solar disk where psi is measured in the counter clock wise direction.
* Heliocentric-Cartesian (HCC): The same as HPC but with positions expressed in
    true (deprojected) physical distances instead of angles on the celestial 
    sphere.
* Heliocentric-Radial (HCR): The same as HPR but with rho expressed in
    true (deprojected) physical distances instead of angles on the celestial 
    sphere.
* Stonyhurst-Heliographic (HG): Expressed positions on the Sun using longitude 
    and latitude on the solar sphere but with the origin which is at the 
    intersection of the solar equator and the central meridian as seen from Earth. 
    This means that the coordinate system remains fixed with respect to Earth 
    while the Sun rotates underneath it.
* Carrington-Heliographic (HG, /CARRINGTON): Carrington longitude is offset 
        from Stonyhurst longitude by a time-dependent scalar value, L0. At the 
        start of each Carrington rotation, L0 = 360, and steadily decreases 
        until it reaches L0 = 0, at which point the next Carrington 
        rotation starts. 

References
----------
* `Thompson (2006), A&A, 449, 791 <http://www.aanda.org/index.php?option=com_article&access=doi&doi=10.1051/0004-6361:20054262&Itemid=129>` 
    `PDF <http://fits.gsfc.nasa.gov/wcs/coordinates.pdf>`

Note that SOLAR_B0, HGLT_OBS, and CRLT_OBS are all synonyms.
"""

__authors__ = ["Steven Christe"]
__email__ = "steven.d.christe@nasa.gov"

import numpy as np
from sunpy.sun import constants as con
from decimal import *

def get_solar_limb(header):
    """Return the angular size of the Sun viewed from Earth (in arcsec)"""
    # khughitt: Perhaps rsun should be handled in source-specific logic, and
    #           passed in?
    return (header.get('RSUN_OBS') or 
            header.get('SOLAR_R') or 
            header.get('RADIUS', con.average_angular_size))

def get_observer_position(header):
    """Return the observer distance from the Sun."""
    return header.get('DSUN_OBS')

def get_center(header, axis=None):
    """Return the center of the map."""
    x = (header.get('cdelt1') * header.get('naxis1') / 2 + 
         header.get('crval1') - header.get('crpix1') * header.get('cdelt1'))
    
    y = (header.get('cdelt2') * header.get('naxis2') / 2 + 
         header.get('crval2') - header.get('crpix2') * header.get('cdelt2'))
    
    if axis is 'x':
        return x
    elif axis is 'y':
        return y
    else:
        return [x,y]
    
def get_units(header, axis=None):
    """Return the units used for crpix, crdelt in the header."""
    xunits = header.get('cunit1', header.get('ctype1'))
    yunits = header.get('cunit2', header.get('ctype2'))
    
    if axis is 'x':
        return xunits
    elif axis is 'y':
        return yunits
    else:
        return [xunits,yunits]
    
def get_platescale(header, axis=None):
    """Return the plate scale of the image, 
    i.e. the size of the pixels in unit."""
    xscale = header.get('cdelt1')
    yscale = header.get('cdelt2')
     
    if axis is 'x':
        return xscale
    elif axis is 'y':
        return yscale
    else:
        return [xscale,yscale]
    
def get_solar_b0(header):
    """Return the solar B0 angle which is simply the heliographic latitude of 
    the observer."""
    
    return (header.get('HGLT_OBS') or
            header.get('CRLT_OBS') or
            header.get('SOLAR_B0', 0))

def get_solar_l0(header, carrington=False):
    """Return the Carrington heliographic longitude of the observer."""
    if carrington is False:
        return header.get('hgln_obs', 0)    
    if carrington is True:
        return header.get('CRLN_OBS', 0)
    
def convert_ang_units(type='hpc', unit='arcsec'):
    """Determine the conversion factor between the data and radians."""
    
    if unit == 'arcmin':
        return np.deg2rad(1) / 60.0
    elif unit == 'arcsec':
        return np.deg2rad(1) / (60 * 60.0)
    elif unit == 'mas':
        return np.deg2rad(1) / (60 * 60 * 1000.0)

def get_projection(header, axis='x'):
    """Return the projection that the data was taken in."""
    xtype = header.get('ctype1')
    ytype = header.get('ctype2')

    # Assume that the projection is the same in both axis
    # TODO: Remove assumption of same projection in both axis     
    if axis is 'x':
        return xtype
    elif axis is 'y':
        return ytype
    else:
        return xtype

def get_shape(header):
    """Return the shape of the data array."""
    return [header.get('naxis1'), header.get('naxis2')]

def convert_pixel_to_data(header, pixel_index=None):
    """This procedure takes a WCS-compliant header, and calculates the data coordinates at each pixel position."""

    naxis = np.array(get_shape(header))
    cdelt = np.array(get_platescale(header))
    crpix = np.array([header.get('crpix1'), header.get('crpix2')])
    
    # first assume that coord is just [x,y]
    if pixel_index is not None:
        pixel_index = np.array(pixel_index)
        coord = (pixel_index - (crpix - 1)) * cdelt
    else:
        tempa = (np.arange(naxis[0]*naxis[1]) %  naxis[0])
        tempb = tempa.reshape(naxis[0],naxis[1]).transpose().reshape(naxis[0]*naxis[1])
        pixel = np.array(zip(tempa,tempb))
        coord = (pixel - (crpix - 1) )* cdelt
            
    # check to see what projection is being used
    projection = get_projection(header)
    if  projection.count('TAN'):    
        coord = proj_tan(header, coord)
        
    return coord

def convert_data_to_pixel(header, coord):
    """This procedure takes a WCS-compliant header, and calculates the pixel coordinates for given data coordinates."""
    naxis = np.array(get_shape(header))
    cdelt = np.array(get_platescale(header))
    crpix = np.array([header.get('crpix1'), header.get('crpix2')])
    # De-apply any tabular projections.
    # coord = inv_proj_tan(header,coord)
    
    pixel = coord/cdelt + crpix - 1

    return pixel

def convert_hpc_hcc(header, hp_longitude, hp_latitude, distance=None):
    """This routine converts Helioprojective-Cartesian (HPC) coordinates into 
    Heliocentric-Cartesian (HCC) coordinates, using equations 15 in 
    Thompson (2006), A&A, 449, 791-803.
    """
    
    c = np.array([convert_ang_units(unit=get_units(header, axis='x')), 
                  convert_ang_units(unit=get_units(header, axis='y'))])

    longitude = hp_longitude * c[0]
    latitude = hp_latitude * c[1]

    cosx = np.cos(longitude)
    sinx = np.sin(longitude)
    cosy = np.cos(latitude)
    siny = np.sin(latitude)

    dsun = header.get('dsun_obs')
    rsun = header.get('rsun_ref')

    if distance is None: 
        q = dsun * cosy * cosx
        distance = q ** 2 - dsun ** 2 + rsun ** 2
        # distance[np.where(distance < 0)] = np.sqrt(-1)
        distance = q - np.sqrt(distance)

    x = distance * cosy * sinx
    y = distance * siny
    z = dsun - distance * cosy * cosx

    return np.array([x,y]).T

def convert_hcc_hpc(header, x, y, distance=None):
    """Convert Heliocentric-Cartesian (HCC) to angular 
    Helioprojective-Cartesian (HPC) coordinates (in degrees)."""

    #Distance to the Sun but should we use our own?
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
    result = np.rad2deg([hpln, hplt]).T
    return result

def convert_hcc_hg(header, x, y):
    """Convert Heliocentric-Cartesian (HCC) to Heliographic coordinates HG (given in degrees)."""

    rsun = header.get('rsun_ref')

    z = np.sqrt(rsun ** 2 - x ** 2 - y ** 2)
    # z[np.where(z < 0)] = Decimal('Nan')

    b0 = get_solar_b0(header)
    l0 = get_solar_l0(header)
    cosb = np.cos(np.deg2rad(b0))
    sinb = np.sin(np.deg2rad(b0))

    hecr = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    hgln = np.arctan(x / (z*cosb - y*sinb)) + l0
    hglt = np.arcsin( (y*cosb + z*sinb) / hecr )
    
    return np.rad2deg(np.array([hgln,hglt]).T)

def convert_hg_hcc(header, hgln, hglt):
    """Convert Heliographic coordinates (given in degrees) to Heliocentric-Cartesian."""
    # using equations 11 in Thompson (2006), A&A, 449, 791-803

    cx = np.deg2rad(1)
    cy = np.deg2rad(1)
    
    hecr = header.get('rsun_ref')
    
    lon = cx * hgln
    lat = cy * hglt
    
    b0 = get_solar_b0(header)
    l0 = get_solar_l0(header)

    cosb = np.cos(np.deg2rad(b0))
    sinb = np.sin(np.deg2rad(b0))

    lon = lon - l0

    cosx = np.cos(lon)
    sinx = np.sin(lon)
    cosy = np.cos(lat)
    siny = np.sin(lat)
    
    # Perform the conversion.
    x = hecr * cosy * sinx
    y = hecr * (siny*cosb - cosy*cosx*sinb)
    z = hecr * (siny*sinb + cosy*cosx*cosb)

    return np.array([x,y]).T

def convert_hg_hpc(header, hglon, hglat):
    """Convert Helioprojective-Cartesian (HPC) to Heliographic coordinates 
    (HG)"""
    temp_result = convert_hg_hcc(header, hglon, hglat)
    result = convert_hcc_hpc(header, temp_result[0], temp_result[1])

    return result

def proj_tan(header, coord, force=False):
    """Applies the gnomonic (TAN) projection to intermediate relative 
    coordinates."""
    # if pixels are within 3 degrees of the Sun then skip the calculatin unless 
    # force is True. This applies to all sdo images so this function is just 
    # here as a place holder
    return coord
    
def convert_to_coord( header, x, y, fromcoord = 'hpc', tocoord = 'hg'):
    """Apply a coordinate transform to coordinate in header 
    to coordinate coord. Right now can only do hpc to hcc to hg"""
    
    #coord = np.array(convert_pixel_to_data(header))
    print(coord.shape)
    #temp = np.array(convert_hpc_hcc(header, coord[:,:,0], coord[:,:,1]))
    print(temp.shape)
    new_coord = np.array(convert_hcc_hg(header, temp[0,:,:], temp[1,:,:]))
                
    return new_coord