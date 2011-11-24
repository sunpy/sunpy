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
from __future__ import absolute_import

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
    
def convert_angle_units(type='hpc', unit='arcsec'):
    """Determine the conversion factor between the data and radians."""
    
    if unit == 'deg':
        return np.deg2rad(1)
    elif unit == 'arcmin':
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

def convert_pixel_to_data(header, x = None, y = None):
    """This procedure takes a WCS-compliant header, and calculates the 
        data coordinates at each x and y pixels. If no x and y are given
        then return the entire detector."""

    naxis = np.array(get_shape(header))
    cdelt = np.array(get_platescale(header))
    crpix = np.array([header.get('crpix1'), header.get('crpix2')])
    crval = np.array([header.get('crval1'), header.get('crval2')])
    
    # first assume that coord is just [x,y]
    if (x is None) and (y is None):
        x, y = np.meshgrid(np.arange(get_shape(header)[0]), np.arange(get_shape(header)[1]))

    coordx = (x - (crpix[0] - 1) ) * cdelt[0] + crval[0]
    coordy = (y - (crpix[1] - 1) ) * cdelt[1] + crval[1]
            
    # check to see what projection is being used
    projection = get_projection(header)
    if  projection.count('TAN'):    
        coordx, coordy = proj_tan(header, coordx, coordy)
        
    return coordx, coordy

def convert_data_to_pixel(header, x, y):
    """This procedure takes a WCS-compliant header, and calculates the pixel 
    coordinates for given data coordinates."""
    # TODO: Needs to check what coordinate system the data is given in
    naxis = np.array(get_shape(header))
    cdelt = np.array(get_platescale(header))
    crpix = np.array([header.get('crpix1'), header.get('crpix2')])
    # De-apply any tabular projections.
    # coord = inv_proj_tan(header,coord)
    
    pixelx = x/cdelt[0] + crpix[1] - 1
    pixely = y/cdelt[1] + crpix[1] - 1

    return pixelx, pixely

def convert_hpc_hcc(header, hp_longitude, hp_latitude, distance=None):
    """This routine converts Helioprojective-Cartesian (HPC) coordinates into 
    Heliocentric-Cartesian (HCC) coordinates, using equations 15 in 
    Thompson (2006), A&A, 449, 791-803.
    """
    
    c = np.array([convert_angle_units(unit=get_units(header, axis='x')), 
                  convert_angle_units(unit=get_units(header, axis='y'))])

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

    return x, y

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
    return np.rad2deg(hpln), np.rad2deg(hplt)

def convert_hcc_hg(header, x, y):
    """Convert Heliocentric-Cartesian (HCC) to Heliographic coordinates HG 
    (given in degrees)."""

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
    
    return np.rad2deg(hgln), np.rad2deg(hglt)

def convert_hg_hcc(header, hgln, hglt, occultation = False):
    """Convert Heliographic coordinates (given in degrees) to Heliocentric-Cartesian."""
    x, y, z = convert_hg_hcc_xyz(header, hgln, hglt)
    
    if occultation:
        index = (z < 0)
        x[index] = np.nan
        y[index] = np.nan
    
    return x, y


def convert_hg_hcc_xyz(header, hgln, hglt):
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
    
    return x, y, z

def test(): 
    # number of points in the line
    grid_spacing = 10
    hg_longitude_deg = np.arange(-90,90, grid_spacing, dtype = 'float')
    hg_latitude_deg = np.arange(-90, 90, grid_spacing, dtype = 'float')

    xx, yy = np.meshgrid(hg_longitude_deg, hg_latitude_deg)

    return xx, yy

def convert_hg_hpc(header, hglon, hglat, units = None, occultation = False):
    """Convert Heliographic coordinates (HG) to Helioprojective-Cartesian (HPC)"""
    tempx, tempy = convert_hg_hcc(header, hglon, hglat, occultation)
    x, y = convert_hcc_hpc(header, tempx, tempy)

    if units == 'arcsec':
        x = 60*60*x
        y = 60*60*y

    return x, y

def convert_hpc_hg(header, x, y):
    """Convert Helioprojective-Cartesian (HPC) to Heliographic coordinates (HG)"""
    tempx, tempy = convert_hpc_hcc(header, x, y)
    lon, lat = convert_hcc_hg(header, tempx, tempy)
    return lon, lat

def proj_tan(header, x, y, force=False):
    """Applies the gnomonic (TAN) projection to intermediate relative 
    coordinates."""
    # if pixels are within 3 degrees of the Sun then skip the calculatin unless 
    # force is True. This applies to all sdo images so this function is just 
    # here as a place holder
    return x, y
    
def convert_to_coord( header, x, y, fromcoord = 'hpc', tocoord = 'hg'):
    """Apply a coordinate transform to coordinates in header 
    to coordinate coord. Right now can only do hpc to hcc to hg"""
     
    #coord = np.array(convert_pixel_to_data(header))
    #temp = np.array(convert_hpc_hcc(header, coord[:,:,0], coord[:,:,1]))
    x, y = convert_hcc_hg(header, x, y)
                
    return x, y