"""
    The WCS package provides functions to parse a World Coordinate System (WCS) 
    coordinates for solar images as well as convert between various solar 
    coordinate systems. The solar coordinates supported are

    * Helioprojective-Cartesian (HPC): The most often used solar coordinate 
        system. Describes positions on the Sun as angles measured from the 
        center of the solar disk (usually in arcseconds) using cartesian 
        coordinates (X, Y)
    * Helioprojective-Radial (HPR): Describes positions on the Sun using angles, 
        similar to HPC, but uses a radial coordinate (rho, psi) system centered 
        on solar disk where psi is measured in the counter clock wise direction.
    * Heliocentric-Cartesian (HCC): The same as HPC but with positions expressed
        in true (deprojected) physical distances instead of angles on the 
        celestial sphere.
    * Heliocentric-Radial (HCR): The same as HPR but with rho expressed in
        true (deprojected) physical distances instead of angles on the celestial 
        sphere.
    * Stonyhurst-Heliographic (HG): Expressed positions on the Sun using 
        longitude and latitude on the solar sphere but with the origin which is 
        at the intersection of the solar equator and the central meridian as 
        seen from Earth. This means that the coordinate system remains fixed 
        with respect to Earth while the Sun rotates underneath it.
    * Carrington-Heliographic (HG, /CARRINGTON): Carrington longitude is offset 
        from Stonyhurst longitude by a time-dependent scalar value, L0. At the 
        start of each Carrington rotation, L0 = 360, and steadily decreases 
        until it reaches L0 = 0, at which point the next Carrington rotation 
        starts. 

    Note that SOLAR_B0, HGLT_OBS, and CRLT_OBS are all synonyms.

    References
    ----------
    | Thompson (2006), A&A, 449, 791 <http://dx.doi.org/10.1051/0004-6361:20054262>
    | PDF <http://fits.gsfc.nasa.gov/wcs/coordinates.pdf>
"""
from __future__ import absolute_import

__all__ = ["convert_angle_units", "convert_pixel_to_data", "convert_hpc_hg",
           "convert_data_to_pixel", "convert_hpc_hcc", "convert_hcc_hpc",
           "convert_hcc_hg", "convert_hg_hcc", "convert_hg_hcc_xyz", "proj_tan",
           "convert_hg_hpc",  "convert_to_coord", "convert_hpc_hcc_xyz", 
           "get_center"]

import numpy as np

def convert_angle_units(unit='arcsec'):
    """Determine the conversion factor between the data and radians."""
    
    if unit == 'deg':
        return np.deg2rad(1)
    elif unit == 'arcmin':
        return np.deg2rad(1) / 60.0
    elif unit == 'arcsec':
        return np.deg2rad(1) / (60 * 60.0)
    elif unit == 'mas':
        return np.deg2rad(1) / (60 * 60 * 1000.0)

def convert_pixel_to_data(width, height, scale_x, scale_y, crpix1, crpix2, 
                          crval1, crval2, ctype, x=None, y=None):
    """This procedure takes WCS-compliant header tags, and calculates the 
        data coordinates at each x and y pixel centers. If no x and y are given
        then return the entire detector."""
    cdelt = np.array([scale_x, scale_y])
    crpix = np.array([crpix1, crpix2])
    crval = np.array([crval1, crval2])
    
    # first assume that coord is just [x,y]
    if (x is None) and (y is None):
        x, y = np.meshgrid(np.arange(width), np.arange(height))

    # note that crpix[] counts pixels starting at 1
    coordx = (x - (crpix[0] - 1)) * cdelt[0] + crval[0]
    coordy = (y - (crpix[1] - 1)) * cdelt[1] + crval[1]
    
    # check to see what projection is being used
    if ctype.count('TAN'):    
        coordx, coordy = proj_tan(coordx, coordy)

    return coordx, coordy

def convert_data_to_pixel(scale_x, scale_y, crpix1, crpix2, 
                          crval1, crval2, x, y):
    """This procedure takes WCS-compliant header tags, and calculates the pixel 
    coordinates for given data coordinates."""
    # TODO: Needs to check what coordinate system the data is given in
    cdelt = np.array([scale_x, scale_y])
    crpix = np.array([crpix1, crpix2])
    crval = np.array([crval1, crval2])
    # De-apply any tabular projections.
    # coord = inv_proj_tan(coord)
    
    # note that crpix[] counts pixels starting at 1
    pixelx = (x - crval[0]) / cdelt[0] + (crpix[1] - 1)
    pixely = (y - crval[1]) / cdelt[1] + (crpix[1] - 1)

    return pixelx, pixely

def convert_hpc_hcc(rsun, dsun, units_x, units_y, hpx, hpy, 
                    distance=None):
    """This routine converts Helioprojective-Cartesian (HPC) coordinates into 
    Heliocentric-Cartesian (HCC) coordinates, using equations 15 in 
    Thompson (2006), A&A, 449, 791-803. Returns only x and y.
    """
    x, y, z = convert_hpc_hcc_xyz(rsun, dsun, units_x, units_y, hpx, hpy)
    return x, y

def convert_hpc_hcc_xyz(rsun, dsun, units_x, units_y, hpx, hpy, distance=None):
    """This routine converts Helioprojective-Cartesian (HPC) coordinates into 
    Heliocentric-Cartesian (HCC) coordinates, using equations 15 in 
    Thompson (2006), A&A, 449, 791-803.
    """
    c = np.array([convert_angle_units(unit=units_x), 
                  convert_angle_units(unit=units_y)])

    cosx = np.cos(hpx * c[0])
    sinx = np.sin(hpx * c[0])
    cosy = np.cos(hpy * c[1])
    siny = np.sin(hpy * c[1])

    if distance is None: 
        q = dsun * cosy * cosx
        distance = q ** 2 - dsun ** 2 + rsun ** 2
        # distance[np.where(distance < 0)] = np.sqrt(-1)
        distance = q - np.sqrt(distance)

    x = distance * cosy * sinx
    y = distance * siny
    z = dsun - distance * cosy * cosx

    return x, y, z

def convert_hcc_hpc(rsun, dsun, x, y, units=None, distance=None):
    """Convert Heliocentric-Cartesian (HCC) to angular 
    Helioprojective-Cartesian (HPC) coordinates (in degrees)."""

    # Should we use the rsun_ref defined in the fits file or our 
    # local (possibly different/more correct) definition
    
    # Calculate the z coordinate by assuming that it is on the surface of the 
    # Sun
    z = np.sqrt(rsun ** 2 - x ** 2 - y ** 2)

    zeta = dsun - z
    distance = np.sqrt(x**2 + y**2 + zeta**2)
    hpcx = np.rad2deg(np.arctan2(x, zeta))
    hpcy = np.rad2deg(np.arcsin(y / distance))
    
    if units == 'arcsec':
        hpcx = 60 * 60 * hpcx
        hpcy = 60 * 60 * hpcy
    
    return hpcx, hpcy

def convert_hcc_hg(rsun, b0, l0, x, y, z=None):
    """Convert Heliocentric-Cartesian (HCC) to Heliographic coordinates (HG) 
    (given in degrees)."""
    if z is None:
        z = np.sqrt(rsun**2 - x**2 - y**2)
    # z[z < 0] = np.NAN

    b0 = np.deg2rad(b0)
    l0 = np.deg2rad(l0)
    cosb = np.cos(b0)
    sinb = np.sin(b0)

    hecr = np.sqrt(x**2 + y**2 + z**2)
    hgln = np.arctan2(x, z * cosb - y * sinb) + l0
    hglt = np.arcsin((y * cosb + z * sinb) / hecr)
    
    return np.rad2deg(hgln), np.rad2deg(hglt)

def convert_hg_hcc(rsun, b0, l0, hgln, hglt, occultation=False):
    """Convert Heliographic coordinates (given in degrees) to 
    Heliocentric-Cartesian."""
    x, y, z = convert_hg_hcc_xyz(rsun, b0, l0, hgln, hglt)
    
    if occultation:
        index = (z < 0)
        x[index] = np.nan
        y[index] = np.nan
    
    return x, y

def convert_hg_hcc_xyz(rsun, b0, l0, hgln, hglt):
    """Convert Heliographic coordinates (given in degrees) to 
    Heliocentric-Cartesian."""
    # using equations 11 in Thompson (2006), A&A, 449, 791-803

    cx = np.deg2rad(1)
    cy = np.deg2rad(1)
    
    lon = cx * hgln
    lat = cy * hglt
    
    b0 = np.deg2rad(b0)
    l0 = np.deg2rad(l0)

    cosb = np.cos(b0)
    sinb = np.sin(b0)

    lon = lon - l0

    cosx = np.cos(lon)
    sinx = np.sin(lon)
    cosy = np.cos(lat)
    siny = np.sin(lat)
    
    # Perform the conversion.
    x = rsun * cosy * sinx
    y = rsun * (siny * cosb - cosy * cosx * sinb)
    z = rsun * (siny * sinb + cosy * cosx * cosb)
    
    return x, y, z

def convert_hg_hpc(rsun, dsun, b0, l0, hglon, hglat, units=None, 
                   occultation=False):
    """Convert Heliographic coordinates (HG) to Helioprojective-Cartesian 
    (HPC)"""
    tempx, tempy = convert_hg_hcc(rsun, b0, l0, hglon, hglat, occultation)
    x, y = convert_hcc_hpc(rsun, dsun, tempx, tempy, units=units)
    return x, y

def convert_hpc_hg(rsun, dsun, units_x, units_y, b0, l0, x, y):
    """Convert Helioprojective-Cartesian (HPC) to Heliographic coordinates 
    (HG)"""
    tempx, tempy = convert_hpc_hcc(rsun, dsun, units_x, units_y, x, y)
    lon, lat = convert_hcc_hg(rsun, b0, l0, tempx, tempy)
    return lon, lat

def get_center(size, scale, crpix, crval):
    """Returns the center of the map."""
    return scale * (size - 1) / 2. + crval - (crpix - 1) * scale

def proj_tan(x, y, force=False):
    """Applies the gnomonic (TAN) projection to intermediate relative 
    coordinates."""
    # if pixels are within 3 degrees of the Sun then skip the calculatin unless 
    # force is True. This applies to all sdo images so this function is just 
    # here as a place holder for the future
    # TODO: write proj_tan function
    return x, y
    
def convert_to_coord(x, y, rsun, b0, l0, fromto):
    """Apply a coordinate transform to coordinates. Right now can only do hpc 
    to hcc to hg"""
    
    #coord = np.array(convert_pixel_to_data(...))
    #temp = np.array(convert_hpc_hcc(rsun, coord[:,:,0], coord[:,:,1]))
    x, y = convert_hcc_hg(rsun, b0, l0, x, y)
                
    return x, y
