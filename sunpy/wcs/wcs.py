from __future__ import absolute_import

import numpy as np
import sunpy.sun as sun

import astropy.units as u

rsun_meters = sun.constants.radius.si

__all__ = ['convert_pixel_to_data', 'convert_hpc_hg',
           'convert_data_to_pixel', 'convert_hpc_hcc', 'convert_hcc_hpc',
           'convert_hcc_hg', 'convert_hg_hcc', 'proj_tan',
           'convert_hg_hpc',  'convert_to_coord',
           'get_center']


def convert_pixel_to_data(size, scale, reference_pixel,
                          reference_coordinate, x=None, y=None):
    """Calculate the data coordinate for particular pixel indices.

    Parameters
    ----------
    size : ~'astropy.units.instance', 2d array
        Number of pixels in width and height.
    scale : ~'astropy.units.instance', 2d array
        The size of a pixel (dx,dy) in data coordinates (equivalent to WCS/CDELT)
    reference_pixel : ~'astropy.units.instance', 2d array
        The reference pixel (x,y) at which the reference coordinate is given (equivalent to WCS/CRPIX)
    reference_coordinate : ~'astropy.units.instance', 2d array
        The data coordinate (x, y) as measured at the reference pixel (equivalent to WCS/CRVAL)
    x,y : ~'astropy.units.instance'
        The pixel values at which data coordinates are requested. If none are given,
        returns coordinates for every pixel.

    Returns
    -------
    out : ~'astropy.units.instance'
        The data coordinates at pixel (x,y).

    Notes
    -----
    This function assumes a gnomic projection which is correct for a detector at the focus
    of an optic observing the Sun.

    Examples
    --------

    """
    if not isinstance(scale or reference_pixel or reference_coordinate
                      or size, u.Quantity):
        raise ValueError("Must be astropy.units instance")
    
    if (x != None) and (y != None):
        if not isinstance(x or y, u.Quantity):
            raise ValueError("Must be astropy.units instance")
    
    if (x is None) and (y is None):
        # first assume that coord is just [x,y]
        x, y = np.meshgrid(np.arange(size[0]), np.arange(size[1]))
        x *= u.pix
        y *= u.pix

    cdelt = np.array(scale) * u.Unit(scale.unit)
    crpix = np.array(reference_pixel) * u.Unit(reference_pixel.unit)
    crval = np.array(reference_coordinate) * u.Unit(reference_coordinate.unit)


    # note that crpix[] counts pixels starting at 1

    coordx = ((x - (crpix[0] - 1*u.pix))).value * cdelt[0] + crval[0]
    coordy = ((y - (crpix[1] - 1*u.pix))).value * cdelt[1] + crval[1]

    # Correct for Gnomic projection
    coordx, coordy = proj_tan(coordx, coordy)

    return np.array([coordx.value, coordy.value]) * u.arcsec

def get_center(size, scale, reference_pixel, reference_coordinate):
    """Returns the center of the image in data coordinates.

    Parameters
    ----------
    size : ~'astropy.units.instance'
        Number of pixels in width and height.
    scale : ~'astropy.units.instance'
        The size of a pixel (dx,dy) in data coordinates (equivalent to WCS/CDELT)
    reference_pixel : ~'astropy.units.instance', 2d array
        The reference pixel (x,y) at which the reference coordinate is given (equivalent to WCS/CRPIX)
    reference_coordinate : ~'astropy.units.instance', array
        The data coordinate (x, y) as measured at the reference pixel (equivalent to WCS/CRVAL)

    Returns
    -------
    out : ~'astropy.units.instance'
        The data coordinates

    Examples
    --------

    """
    if not isinstance(size or scale or reference_pixel or reference_coordinate,
                       u.Quantity):
        raise ValueError("Must be astropy.quantity instance")
    return (scale.value * (size.value - 1) / 2. + reference_coordinate.value 
           - (reference_pixel.value - 1) * scale.value) * scale.unit

def convert_data_to_pixel(x, y, scale, reference_pixel, reference_coordinate):
    """Calculate the pixel indices for a given data coordinate.

    Parameters
    ----------
    x, y : ~'astropy.units.instance'
        Data coordinate in same units as reference coordinate
    scale : ~'astropy.units.instance'
        The size of a pixel (dx,dy) in data coordinates (equivalent to WCS/CDELT)
    reference_pixel : ~'astropy.units.instance'
        The reference pixel (x,y) at which the reference coordinate is given (equivalent to WCS/CRPIX)
    reference_coordinate : ~'astropy.units.instance'
        The data coordinate (x, y) as measured at the reference pixel (equivalent to WCS/CRVAL)

    Returns
    -------
    out : ~'astropy.units.instance'
        The  pixel coordinates (x,y) at that data coordinate.

    Examples
    --------

    """
    
    if not isinstance(scale or reference_pixel or reference_coordinate,
                      u.Quantity):
        raise ValueError("Must be astropy.units instance")
    
    if not isinstance(x or y, u.Quantity):
        raise ValueError("Must be astropy.units instance")
    # TODO: Needs to check what coordinate system the data is given in
    cdelt = np.array(scale) * u.Unit(scale.unit)
    crpix = np.array(reference_pixel) * u.Unit(reference_pixel.unit)
    crval = np.array(reference_coordinate) * u.Unit(reference_coordinate.unit)
    # De-apply any tabular projections.l
    # coord = inv_proj_tan(coord)

    # note that crpix[] counts pixels starting at 1
    pixelx = ((x - crval[0]) / cdelt[0]) * u.pix + (crpix[1] - 1 * u.pix)
    pixely = ((y - crval[1]) / cdelt[1]) * u.pix + (crpix[1] - 1 * u.pix)

    return np.array([pixelx.value, pixely.value]) * u.pix

def convert_hpc_hcc(x, y, dsun_meters=None, z=False):
    """Converts from Helioprojective-Cartesian (HPC) coordinates into
    Heliocentric-Cartesian (HCC) coordinates. Returns all three dimensions, x, y, z in
    meters.

    Parameters
    ----------
    x, y : ~'astropy.units.instance'
        Data coordinate in angle units
    dsun_meters : ~'astropy.units.instance'
        Distance from the observer to the Sun in meters. Default is 1 AU.
    z : Bool
        If true return the z coordinate as well.

    Returns
    -------
    out : ~'astropy.units.instance'
        The data coordinates (x,y,z) in heliocentric cartesian coordinates in meters.

    Notes
    -----
    Implements Eq. (15) of Thompson (2006), A&A, 449, 791.

    Examples
    --------
    >>> sunpy.wcs.convert_hpc_hcc(40.0, 32.0, z=True)
    (28876152.176423457, 23100922.071266972, 694524220.8157959)

    """
    rsun_meters = sun.constants.radius.si
    if not isinstance(x or y, u.Quantity):
        raise ValueError("Must be astropy.Quantity instance")
    cosx = np.cos(x)
    sinx = np.sin(x)
    cosy = np.cos(y)
    siny = np.sin(y)

    if dsun_meters is None:
        dsun_meters = sun.constants.au.si
    if not isinstance(dsun_meters, u.Quantity):
        raise ValueError("Must be astropy.Quantity instance")

    q = dsun_meters * cosy * cosx
    distance = q ** 2 - dsun_meters ** 2 + rsun_meters ** 2
    # distance[np.where(distance < 0)] = np.sqrt(-1)
    distance = q - np.sqrt(distance)

    rx = distance * cosy * sinx
    ry = distance * siny
    rz = dsun_meters - distance * cosy * cosx

    if np.all(z == True):
        return rx, ry, rz
    else:
        return rx, ry

def convert_hcc_hpc(x, y, dsun_meters=None):
    """Convert Heliocentric-Cartesian (HCC) to angular
    Helioprojective-Cartesian (HPC) coordinates (in degrees).

    Parameters
    ----------
    x, y : ~'astropy.units.instance'
        Data coordinate in meters.
    dsun_meters : ~'astropy.units.instance'
        Distance from the observer to the Sun in meters. Default is 1 AU.

    Returns
    -------
    out : ~'astropy.units.instance'
        The  data coordinates (x,y) in helioprojective cartesian coordinates in degrees.

    Notes
    -----
    Implements Eq. (16) of Thompson (2006), A&A, 449, 791.

    Examples
    --------
    >>> sunpy.wcs.convert_hcc_hpc(28748691, 22998953)
    (39.823439773829705, 31.858751644835717)

    """

    rsun_meters = sun.constants.radius.si
    if not isinstance(x and y, u.Quantity):
       raise ValueError("Must be astropy.Quantity instance")

    if x.unit != 'meter':
       x = x.to('meter')
    if y.unit != 'meter':
       y = y.to('meter')

    # Calculate the z coordinate by assuming that it is on the surface of the Sun
    z = np.sqrt(rsun_meters ** 2 - x ** 2 - y ** 2)

    if dsun_meters is None:
        dsun_meters = sun.constants.au.si
    if not isinstance(dsun_meters, u.Quantity):
        raise ValueError("Must be astropy.units.instance")
    dsun_meters = dsun_meters.to('meter')

    zeta = dsun_meters - z
    distance = np.sqrt(x**2 + y**2 + zeta**2)
    hpcx = np.arctan2(x, zeta).to(u.deg)
    hpcy = np.arcsin(y / distance).to(u.deg)

    return hpcx, hpcy

def convert_hcc_hg(x, y, z=None, b0_deg=0 * u.deg, l0_deg=0 * u.deg, radius=False):
    """Convert from Heliocentric-Cartesian (HCC) (given in meters) to
    Stonyhurst Heliographic coordinates (HG) given in degrees, with
    radial output in meters.

    Parameters
    ----------
    x, y : ~'astropy.units.instance'
        Data coordinate in meters.
    z : ~'astropy.units.instance'
        Data coordinate in meters.  If None, then the z-coordinate is assumed
        to be on the Sun.
    b0_deg : ~'astropy.units.instance'
        Tilt of the solar North rotational axis toward the observer
        (heliographic latitude of the observer). Usually given as SOLAR_B0,
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0_deg : ~'astropy.units.instance'
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    radius : Bool
        If true, forces the output to return a triple of (lon, lat, r). If
        false, return (lon, lat) only.

    Returns
    -------
    out : ndarray (degrees, meters)
        if radius is false, return the data coordinates (lon, lat).  If
        radius=True, return the data cordinates (lon, lat, r).  The quantities
        (lon, lat) are the heliographic coordinates in degrees.  The quantity
        'r' is the heliographic radius in meters.

    Notes
    -----
    Implements Eq. (12) of Thompson (2006), A&A, 449, 791.

    Examples
    --------
    >>> sunpy.wcs.convert_hcc_hg(230000.0,45000000.0,
    z=695508000.0 + 8000000.0, radius=True)
    (0.01873188196651189, 3.6599471896203317, 704945784.41465974)
    """

    rsun_meters = sun.constants.radius.si
    if not isinstance(x and y, u.Quantity):
        raise ValueError("Must be astropy.units.instance")
    x = x.to(u.meter)
    y = y.to(u.meter)
    if z is None:
        z = np.sqrt(rsun_meters**2 - x**2 - y**2)
    if not isinstance(z, u.Quantity):
        raise ValueError("Must be astropy.units.instance")
    if not isinstance(b0_deg and l0_deg, u.Quantity):
        raise ValueError("Must be astropy.units.instance")

    cosb = np.cos(b0_deg)
    sinb = np.sin(b0_deg)

    hecr = np.sqrt(x**2 + y**2 + z**2)
    hgln = np.arctan2(x, z * cosb - y * sinb) + l0_deg.to(u.rad)
    hglt = np.arcsin((y * cosb + z * sinb) / hecr)

    if radius:
        return hgln.to(u.deg), hglt.to(u.deg), hecr
    else:
        return hgln.to(u.deg), hglt.to(u.deg)

def convert_hg_hcc(hglon_deg, hglat_deg, b0_deg=0 * u.deg, l0_deg=0 * u.deg, occultation=False,
                   z=False, r=sun.constants.radius.si):
    """Convert from Stonyhurst Heliographic coordinates (given in degrees) to
    Heliocentric-Cartesian coordinates (given in meters).

    Parameters
    ----------
    hglon_deg, hglat_deg : ~'astropy.units.instance'
        Heliographic longitude and latitude in degrees.
    b0_deg : ~'astropy.units.instance'
        Tilt of the solar North rotational axis toward the observer
        (heliographic latitude of the observer). Usually given as SOLAR_B0,
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0_deg : ~'astropy.units.instance'
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    occultation : Bool
        If true set all points behind the Sun (e.g. not visible) to Nan.
    z : Bool
        If true return the z coordinate as well.
    r : ~'astropy.units.instance'
        Heliographic radius

    Returns
    -------
    out : ndarray (meters)
        The data coordinates in Heliocentric-Cartesian coordinates.

    Notes
    -----
    Implements Eq. (11) of Thompson (2006), A&A, 449, 791, with the default
    assumption that the value 'r' in Eq. (11) is identical to the radius of the
    Sun.

    Examples
    --------
    >>> sunpy.wcs.convert_hg_hcc(0.01873188196651189, 3.6599471896203317,
    r=704945784.41465974, z=True)
    (230000.0, 45000000.0, 703508000.0)
    """
    if not isinstance(hglon_deg or hglat_deg, u.Quantity):
        raise ValueError("Must be astropy.units.Quantity instance")
    if not isinstance(b0_deg or l0_deg, u.Quantity):
        raise ValueError("Must be astropy.units.Quantity instance")
    lon = hglon_deg
    lat = hglat_deg

    cosb = np.cos(b0_deg)
    sinb = np.sin(b0_deg)

    lon = lon - l0_deg

    cosx = np.cos(lon)
    sinx = np.sin(lon)
    cosy = np.cos(lat)
    siny = np.sin(lat)

    # Perform the conversion.
    x = np.array([(r * cosy * sinx).value])
    y = np.array([(r * (siny * cosb - cosy * cosx * sinb)).value])
    zz = np.array([(r * (siny * sinb + cosy * cosx * cosb)).value])

    if occultation:
        x[zz < 0] = np.nan
        y[zz < 0] = np.nan

    if np.all(z == True):
        return [x, y, zz] * u.meter
    else:
        return [x, y] * u.meter

def convert_hg_hpc(hglon_deg, hglat_deg, b0_deg=0 * u.deg, l0_deg=0 * u.deg, dsun_meters=None,
                   occultation=False):
    """Convert from Heliographic coordinates (HG) to Helioprojective-Cartesian
    (HPC).

    Parameters
    ----------
    hglon_deg, hglat_deg : ~'astropy.units.instance'
        Heliographic longitude and latitude in degrees.
    b0_deg : ~'astropy.units.instance'
        Tilt of the solar North rotational axis toward the observer
        (heliographic latitude of the observer). Usually given as SOLAR_B0,
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0_deg : ~'astropy.units.instance'
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    occultation : Bool
        If true set all points behind the Sun (e.g. not visible) to Nan.
    dsun_meters : ~'astropy.units.instance'
        Distance between the observer and the Sun.


    Returns
    -------
    out : ~'astropy.units.Quantity instance'
        The data coordinates (x,y) in Helioprojective-Cartesian coordinates.

    Notes
    -----
    Uses equations 11 and 16 in Thompson (2006), A&A, 449, 791-803.

    Examples
    --------
    >>> sunpy.wcs.convert_hg_hpc(34.0, 45.0, b0_deg=-7.064078, l0_deg=0.0)
    (380.05656560308898, 743.78281283290016)
    """
    
    if not isinstance(hglon_deg or hglat_deg, u.Quantity):
        raise ValueError("Must be astropy.units.Quantity instance")
    if not isinstance(l0_deg or dsun_meters, u.Quantity):
        raise ValueError("Must be astropy.units.Quanitity instance")
    tempx, tempy = convert_hg_hcc(hglon_deg, hglat_deg, b0_deg=b0_deg, l0_deg=l0_deg, occultation=occultation)
    x, y = convert_hcc_hpc(tempx, tempy, dsun_meters=dsun_meters)
    return x, y

def convert_hpc_hg(x, y, b0_deg=0 * u.deg, l0_deg=0 * u.deg, dsun_meters=None):
    """Convert from Helioprojective-Cartesian (HPC) to Heliographic coordinates
    (HG) in degrees.

    Parameters
    ----------
    x, y : ~'astropy.units.instance'
        Data coordinate in angle units.
    b0_deg : ~'astropy.units.instance'
        Tilt of the solar North rotational axis toward the observer
        (heliographic latitude of the observer). Usually given as SOLAR_B0,
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0_deg : ~'astropy.units.instance'
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    dsun_meters : ~'astropy.units.instance'
        Distance between the observer and the Sun.

    Returns
    -------
    out : ~'astropy.units.instance'
        The  data coordinates (hglongitude, hglatitude) in Heliographic coordinates.

    Notes
    -----
    Uses equations 15 and 12 in Thompson (2006), A&A, 449, 791-803.

    Examples
    --------
    >>> sunpy.wcs.convert_hg_hpc(382, 748, b0_deg=-7.064078, l0_deg=0.0)
    (34.504653439914669, 45.443143275518182)
    """

    if not isinstance(x or  y, u.Quantity):
        raise ValueError("Must be astropy.units.Quantity instance")
    if not isinstance(b0_deg or l0_deg, u.Quantity):
        raise ValueError("Must be astropy.units.Quantity instance")   
    tempx, tempy = convert_hpc_hcc(x, y, dsun_meters=dsun_meters)
    lon, lat = convert_hcc_hg(tempx, tempy, b0_deg=b0_deg, l0_deg=l0_deg)
    return lon, lat

def proj_tan(x, y, force=False):
    """Applies the gnomonic (TAN) projection to intermediate relative
    coordinates. This function is not currently implemented!"""
    # if pixels are within 3 degrees of the Sun then skip the calculatin unless
    # force is True. This applies to all sdo images so this function is just
    # here as a place holder for the future
    # TODO: write proj_tan function
    return x, y

def convert_to_coord(x, y, from_coord, to_coord, b0_deg=0 * u.deg, l0_deg=0 * u.deg, dsun_meters=None):
    """Apply a coordinate transform to coordinates. Right now can only do hpc
    to hcc to hg"""

    if not isinstance(x or y, u.Quantity):
        raise ValueError("Must be astropy.units.Quanitity instance")
    if not isinstance(b0_deg or l0_deg, u.Quantity):
        raise ValueError("Must be astropy.units.Quantity instance")
    if (from_coord == 'hcc') and (to_coord == 'hg'):
        rx, ry = convert_hcc_hg(x, y, b0_deg=b0_deg, l0_deg=l0_deg)
    elif (from_coord == 'hpc') and (to_coord == 'hg'):
        rx, ry = convert_hpc_hg(x, y, b0_deg=b0_deg, l0_deg=l0_deg, dsun_meters=dsun_meters)
    elif (from_coord == 'hg') and (to_coord == 'hcc'):
        rx, ry = convert_hg_hcc(x, y, b0_deg=b0_deg, l0_deg=l0_deg)
    elif (from_coord == 'hcc') and (to_coord == 'hpc'):
        rx, ry = convert_hcc_hpc(x, y, dsun_meters=dsun_meters)
    elif (from_coord == 'hg') and (to_coord == 'hpc'):
        rx, ry = convert_hg_hpc(x, y, b0_deg=b0_deg, l0_deg=l0_deg, dsun_meters=dsun_meters)
    elif (from_coord == 'hpc') and (to_coord == 'hcc'):
        rx, ry = convert_hpc_hcc(x, y, dsun_meters=dsun_meters)

    return rx, ry
