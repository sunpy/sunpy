"""
    The WCS package provides functions to parse World Coordinate System (WCS) 
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
    * Carrington-Heliographic (HG): Carrington longitude is offset 
        from Stonyhurst longitude by a time-dependent scalar value, L0. At the 
        start of each Carrington rotation, L0 = 360, and steadily decreases 
        until it reaches L0 = 0, at which point the next Carrington rotation 
        starts. 
        
    Some definition
    
    * b0: Tilt of the solar North rotational axis toward the observer 
        (helio- graphic latitude of the observer). Note that SOLAR_B0, 
        HGLT_OBS, and CRLT_OBS are all synonyms.
    * l0: Carrington longitude of central meridian as seen from Earth.
    * dsun_meters: Distance between observer and the Sun. Default is 1 AU.
    * rsun: Radius of the Sun in meters. Default is 6.955e8 meters. This valued is stored
      locally in this module and can be modified if necessary.
    
    References
    ----------
    | Thompson (2006), A&A, 449, 791 <http://dx.doi.org/10.1051/0004-6361:20054262>
    | PDF <http://fits.gsfc.nasa.gov/wcs/coordinates.pdf>
"""
from __future__ import absolute_import

import numpy as np
import sunpy.sun as sun

rsun_meters = sun.constants.radius

__all__ = ['_convert_angle_units', 'convert_pixel_to_data', 'convert_hpc_hg',
           'convert_data_to_pixel', 'convert_hpc_hcc', 'convert_hcc_hpc',
           'convert_hcc_hg', 'convert_hg_hcc', 'convert_hg_hcc_xyz', 'proj_tan',
           'convert_hg_hpc',  'convert_to_coord', 'convert_hpc_hcc_xyz', 
           'get_center']

def _convert_angle_units(unit='arcsec'):
    """Determine the conversion factor between the data units and radians."""
    if unit == 'deg':
        return np.deg2rad(1)
    elif unit == 'arcmin':
        return np.deg2rad(1) / 60.0
    elif unit == 'arcsec':
        return np.deg2rad(1) / (60 * 60.0)
    elif unit == 'mas':
        return np.deg2rad(1) / (60 * 60 * 1000.0)

def convert_pixel_to_data(size, scale, reference_pixel, 
                          reference_coordinate, x=None, y=None):
    """Calculate the data coordinate for particular pixel indices.
    
    Parameters
    ----------
    size : 2d ndarray
        Number of pixels in width and height.
   scale : 2d ndarray
        The size of a pixel (dx,dy) in data coordinates (equivalent to WCS/CDELT)
    reference_pixel : 2d ndarray
        The reference pixel (x,y) at which the reference coordinate is given (equivalent to WCS/CRPIX)
    reference_coordinate: 2d ndarray
        The data coordinate (x, y) as measured at the reference pixel (equivalent to WCS/CRVAL)
    x,y : int or ndarray
        The pixel values at which data coordinates are requested. If none are given, 
        returns coordinates for every pixel.
        
    Returns
    -------
    out : ndarray
        The data coordinates at pixel (x,y).
    
    Notes
    -----
    This function assumes a gnomic projection which is correct for a detector at the focus
    of an optic observing the Sun.
        
    Examples
    --------
    
    """
    cdelt = np.array(scale)
    crpix = np.array(reference_pixel)
    crval = np.array(reference_coordinate)
    
    # first assume that coord is just [x,y]
    if (x is None) and (y is None):
        x, y = np.meshgrid(np.arange(size[0]), np.arange(size[1]))

    # note that crpix[] counts pixels starting at 1
    coordx = (x - (crpix[0] - 1)) * cdelt[0] + crval[0]
    coordy = (y - (crpix[1] - 1)) * cdelt[1] + crval[1]
    
    # Correct for Gnomic projection
    coordx, coordy = proj_tan(coordx, coordy)
    
    return coordx, coordy

def get_center(size, scale, reference_pixel, reference_coordinate):
    """Returns the center of the image in data coordinates.
      
    Parameters
    ----------
    size : 2d ndarray
        Number of pixels in width and height.
    scale : 2d ndarray
        The size of a pixel (dx,dy) in data coordinates (equivalent to WCS/CDELT)
    reference_pixel : 2d ndarray
        The reference pixel (x,y) at which the reference coordinate is given (equivalent to WCS/CRPIX)
    reference_coordinate: 2d ndarray
        The data coordinate (x, y) as measured at the reference pixel (equivalent to WCS/CRVAL)
          
    Returns
    -------
    out : ndarray
        The data coordinates
            
    Examples
    --------
    
    """
    return scale * (size - 1) / 2. + reference_coordinate - (reference_pixel - 1) * scale

def convert_data_to_pixel(x, y, scale, reference_pixel, reference_coordinate):
    """Calculate the pixel indices for a given data coordinate.

    Parameters
    ----------
    x, y : float
        Data coordinate in same units as reference coordinate
    scale : 2d ndarray
        The size of a pixel (dx,dy) in data coordinates (equivalent to WCS/CDELT)
    reference_pixel : 2d ndarray
        The reference pixel (x,y) at which the reference coordinate is given (equivalent to WCS/CRPIX)
    reference_coordinate: 2d ndarray
        The data coordinate (x, y) as measured at the reference pixel (equivalent to WCS/CRVAL)
        
    Returns
    -------
    out : ndarray
        The  pixel coordinates (x,y) at that data coordinate.
            
    Examples
    --------
    
    """
    
    # TODO: Needs to check what coordinate system the data is given in
    cdelt = np.array([scale_x, scale_y])
    crpix = np.array([reference_pixel_x, reference_pixel_y])
    crval = np.array([reference_coordinate_x, reference_coordinate_y])
    # De-apply any tabular projections.
    # coord = inv_proj_tan(coord)
    
    # note that crpix[] counts pixels starting at 1
    pixelx = (x - crval[0]) / cdelt[0] + (crpix[1] - 1)
    pixely = (y - crval[1]) / cdelt[1] + (crpix[1] - 1)

    return pixelx, pixely

def convert_hpc_hcc(x, y, dsun_meters=None, angle_units='arcsec', distance=None):
    """Converts between Helioprojective-Cartesian (HPC) coordinates (hpx, hpy) into 
    Heliocentric-Cartesian (HCC) coordinates, using equations 15 in 
    Thompson (2006), A&A, 449, 791-803. Returns two dimensional coordinates in meters.

    Parameters
    ----------
    x, y : float
        Data coordinate in angle units (default is arcsec)
    dsun_meters : float
        Distance from the observer to the Sun in meters
    angle_units : str
        Units of the data coordinates (e.g. arcsec, arcmin, degrees). Default is arcsec.
    distance : float
          
    Returns
    -------
    out : ndarray
        The  data coordinates (x,y) in heliocentric cartesian coordinates in meters.

    Notes
    -----

    Examples
    --------
    
    """
    x, y, z = convert_hpc_hcc_xyz(x, y, dsun_meters, angle_units, distance)
    return x, y

def convert_hpc_hcc_xyz(x, y, dsun_meters=None, angle_units='arcsec', distance=None):
    """Converts from Helioprojective-Cartesian (HPC) coordinates into 
    Heliocentric-Cartesian (HCC) coordinates. Returns all three dimensions, x, y, z in
    meters.
    
    Parameters
    ----------
    x, y : float
        Data coordinate in angle units (default is arcsec)
    dsun_meters : float
        Distance from the observer to the Sun in meters. Default is 1 AU.
    angle_units : str
        Units of the data coordinates (e.g. arcsec, arcmin, degrees). Default is arcsec.
    distance : float
        
    Returns
    -------
    out : ndarray
        The  data coordinates (x,y,z) in heliocentric cartesian coordinates in meters.
    
    Notes
    -----
            
    Examples
    --------
    
    """
    
    c = np.array([_convert_angle_units(unit=angle_units), 
                  _convert_angle_units(unit=angle_units)])

    cosx = np.cos(x * c[0])
    sinx = np.sin(x * c[0])
    cosy = np.cos(y * c[1])
    siny = np.sin(y * c[1])

    if dsun_meters is None:
        dsun_meters = sun.constants.au

    if distance is None: 
        q = dsun_meters * cosy * cosx
        distance = q ** 2 - dsun_meters ** 2 + rsun_meters ** 2
        # distance[np.where(distance < 0)] = np.sqrt(-1)
        distance = q - np.sqrt(distance)

    rx = distance * cosy * sinx
    ry = distance * siny
    rz = dsun_meters - distance * cosy * cosx

    return rx, ry, rz

def convert_hcc_hpc(x, y, dsun_meters=None, angle_units='arcsec'):
    """Convert Heliocentric-Cartesian (HCC) to angular 
    Helioprojective-Cartesian (HPC) coordinates (in degrees).

    Parameters
    ----------
    x, y : float
        Data coordinate in meters.
    dsun_meters : float
        Distance from the observer to the Sun in meters. Default is 1 AU.
    angle_units : str
        Units of the data coordinates (e.g. arcsec, arcmin, degrees). Default is arcsec.
    distance : float
        
    Returns
    -------
    out : ndarray
        The  data coordinates (x,y) in helioprojective cartesian coordinates in arcsec.
    
    Notes
    -----
            
    Examples
    --------
    
    """
    
    # Calculate the z coordinate by assuming that it is on the surface of the 
    # Sun
    z = np.sqrt(rsun_meters ** 2 - x ** 2 - y ** 2)
    
    if dsun_meters is None:
        dsun_meters = sun.constants.au
    zeta = dsun_meters - z
    distance = np.sqrt(x**2 + y**2 + zeta**2)
    hpcx = np.rad2deg(np.arctan2(x, zeta))
    hpcy = np.rad2deg(np.arcsin(y / distance))
    
    if angle_units == 'arcsec':
        hpcx = 60 * 60 * hpcx
        hpcy = 60 * 60 * hpcy
    
    return hpcx, hpcy

def convert_hcc_hg(x, y, b0=0, l0=0):
    """Convert from Heliocentric-Cartesian (HCC) (given in arcsec) to Heliographic 
    coordinates (HG) given in degrees.
    
    Parameters
    ----------
    x, y : float (meters)
        Data coordinate in meters. Z unit is assumed to be on the Sun.
    b0 : float (degrees)
        Tilt of the solar North rotational axis toward the observer 
        (heliographic latitude of the observer). Usually given as SOLAR_B0, 
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0 : float (degrees)
        Carrington longitude of central meridian as seen from Earth. Default is 0.
           
    Returns
    -------
    out : ndarray
        The  data coordinates (x,y) in heliographic coordinates in degrees.
    
    Notes
    -----
            
    Examples
    --------
    
    """
    z = np.sqrt(rsun_meters**2 - x**2 - y**2)

    b0 = np.deg2rad(b0)
    l0 = np.deg2rad(l0)
    cosb = np.cos(b0)
    sinb = np.sin(b0)

    hecr = np.sqrt(x**2 + y**2 + z**2)
    hgln = np.arctan2(x, z * cosb - y * sinb) + l0
    hglt = np.arcsin((y * cosb + z * sinb) / hecr)
    
    return np.rad2deg(hgln), np.rad2deg(hglt)

def convert_hg_hcc(hglongitude, hglatitude, b0=0, l0=0, occultation=False):
    """Convert from Heliographic coordinates (given in degrees) to 
    Heliocentric-Cartesian coordinates.
    
    Parameters
    ----------
    hglongitude, hglatitude : float (degrees)
        Data coordinate in meters. Z unit is assumed to be on the Sun.
    b0 : float (degrees)
        Tilt of the solar North rotational axis toward the observer 
        (heliographic latitude of the observer). Usually given as SOLAR_B0, 
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0 : float (degrees)
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    occultation : Bool
        If true set all points behind the Sun (e.g. not visible) to Nan.
    
    Returns
    -------
    out : ndarray (meters)
        The  data coordinates (x,y) in Heliocentric-Cartesian coordinates.
    
    Notes
    -----
            
    Examples
    --------
    
    """
    x, y, z = convert_hg_hcc_xyz(hglongitude, hglatitude, b0, l0)
    
    if occultation:
        index = (z < 0)
        x[index] = np.nan
        y[index] = np.nan
    
    return x, y

def convert_hg_hcc_xyz(hgln, hglt, b0=0, l0=0):
    """Convert from Heliographic coordinates (given in degrees) to 
    Heliocentric-Cartesian coordinates.
    
    Parameters
    ----------
    hglongitude, hglatitude : float (degrees)
        Data coordinate in meters. Z unit is assumed to be on the Sun.
    b0 : float (degrees)
        Tilt of the solar North rotational axis toward the observer 
        (heliographic latitude of the observer). Usually given as SOLAR_B0, 
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0 : float (degrees)
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    occultation : Bool
        If true set all points behind the Sun (e.g. not visible) to Nan.
    
    Returns
    -------
    out : ndarray (meters)
        The  data coordinates (x,y) in Heliocentric-Cartesian coordinates.
    
    Notes
    -----
    
                
    Examples
    --------
    
    """
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
    x = rsun_meters * cosy * sinx
    y = rsun_meters * (siny * cosb - cosy * cosx * sinb)
    z = rsun_meters * (siny * sinb + cosy * cosx * cosb)
    
    return x, y, z

def convert_hg_hpc(hglon, hglat, b0=0, l0=0, dsun_meters=None, angle_units='arcsec', 
                   occultation=False):
    """Convert from Heliographic coordinates (HG) to Helioprojective-Cartesian 
    (HPC).
    
    Parameters
    ----------
    hglongitude, hglatitude : float (degrees)
        Data coordinate in meters. Z unit is assumed to be on the Sun.
    b0 : float (degrees)
        Tilt of the solar North rotational axis toward the observer 
        (heliographic latitude of the observer). Usually given as SOLAR_B0, 
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0 : float (degrees)
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    occultation : Bool
        If true set all points behind the Sun (e.g. not visible) to Nan.
    
    Returns
    -------
    out : ndarray (arcsec)
        The  data coordinates (x,y) in Helioprojective-Cartesian coordinates.
    
    Notes
    -----
    Uses equations 15 in Thompson (2006), A&A, 449, 791-803. 
            
    Examples
    --------
    
    """
    
    tempx, tempy = convert_hg_hcc(hglon, hglat, b0=b0, l0=b0, occultation=occultation)
    x, y = convert_hcc_hpc(tempx, tempy, dsun_meters=dsun_meters, angle_units=angle_units)
    return x, y

def convert_hpc_hg(x, y, b0=0, l0=0, dsun_meters=None, angle_units='arcsec'):
    """Convert from Helioprojective-Cartesian (HPC) to Heliographic coordinates 
    (HG) in degrees.
    
    Parameters
    ----------
    x, y : float ()
        Data coordinate in angle units.
    b0 : float (degrees)
        Tilt of the solar North rotational axis toward the observer 
        (heliographic latitude of the observer). Usually given as SOLAR_B0, 
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0 : float (degrees)
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    dsun_meters : float (meters)
        If true set all points behind the Sun (e.g. not visible) to Nan.
    angle_units : str
        
    
    Returns
    -------
    out : ndarray (degrees)
        The  data coordinates (hglongitude, hglatitude) in Heliographic coordinates.
    
    Notes
    -----
    Uses equations 15 in Thompson (2006), A&A, 449, 791-803. 
            
    Examples
    --------
    
    """
    
    tempx, tempy = convert_hpc_hcc(x, y, dsun_meters, angle_units)
    lon, lat = convert_hcc_hg(tempx, tempy, b0, l0)
    return lon, lat

def proj_tan(x, y, force=False):
    """Applies the gnomonic (TAN) projection to intermediate relative 
    coordinates. This function is not currently implemented!"""
    # if pixels are within 3 degrees of the Sun then skip the calculatin unless 
    # force is True. This applies to all sdo images so this function is just 
    # here as a place holder for the future
    # TODO: write proj_tan function
    return x, y
    
def convert_to_coord(x, y, from_coord, to_coord, b0=0, l0=0):
    """Apply a coordinate transform to coordinates. Right now can only do hpc 
    to hcc to hg"""
    
    if (from_coord == 'hcc') and (to_coord == 'hg'):
         rx, ry = convert_hcc_hg(x, y, b0, l0)
    elif (from_coord == 'hpc') and (to_coord == 'hg'):
         #rx, ry = convert_hpc_hg(x, y, b0, l0)
         convert_hpc_hg(x, y, b0=0, l0=0, dsun_meters=None, angle_units='arcsec')
    elif (from_coord == 'hg') and (to_coord == 'hcc'):
         rx, ry = convert_hcc_hg(x, y, b0, l0)
    elif (from_coord == 'hcc') and (to_coord == 'hpc'):
         rx, ry = convert_hcc_hg(x, y, b0, l0)
    elif (from_coord == 'hg') and (to_coord == 'hpc'):
         rx, ry = convert_hcc_hg(x, y, b0, l0)
    elif (from_coord == 'hpc') and (to_coord == 'hcc'):
         rx, ry = convert_hcc_hg(x, y, b0, l0)
    
    return rx, ry