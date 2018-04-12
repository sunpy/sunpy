from __future__ import absolute_import

import numpy as np

from sunpy.util.decorators import deprecated

import sunpy.sun as sun

import astropy.units as u

rsun_meters = sun.constants.radius.si.value

__all__ = ['_convert_angle_units', 'convert_pixel_to_data', 'convert_hpc_hg',
           'convert_data_to_pixel', 'convert_hpc_hcc', 'convert_hcc_hpc',
           'convert_hcc_hg', 'convert_hg_hcc', 'proj_tan',
           'convert_hg_hpc',  'convert_to_coord',
           'get_center']

def _convert_angle_units(unit='arcsec'):
    """Determine the conversion factor between the data units and radians."""
    if unit == 'degrees':
        return np.deg2rad(1)
    elif unit == 'arcmin':
        return np.deg2rad(1) / 60.0
    elif unit == 'arcsec':
        return np.deg2rad(1) / (60 * 60.0)
    elif unit == 'mas':
        return np.deg2rad(1) / (60 * 60 * 1000.0)
    else:
        raise ValueError("The units specified are either invalid or is not supported at this time.")

@deprecated("0.8.0", alternative="sunpy.map.GenericMap.pixel_to_world")
def convert_pixel_to_data(size, scale, reference_pixel,
                          reference_coordinate, x=None, y=None):
    """
    Calculate the data coordinate for particular pixel indices.

    Parameters
    ----------
    size : 2d ndarray
        Number of pixels in width and height.
    scale : 2d ndarray
        The size of a pixel (dx,dy) in data coordinates (equivalent to WCS/CDELT)
    reference_pixel : 2d ndarray
        The reference pixel (x,y) at which the reference coordinate is given (equivalent to WCS/CRPIX)
    reference_coordinate : 2d ndarray
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

@deprecated("0.8.0")
def get_center(size, scale, reference_pixel, reference_coordinate):
    """
    Returns the center of the image in data coordinates.

    Parameters
    ----------
    size : 2d ndarray
        Number of pixels in width and height.
    scale : 2d ndarray
        The size of a pixel (dx,dy) in data coordinates (equivalent to WCS/CDELT)
    reference_pixel : 2d ndarray
        The reference pixel (x,y) at which the reference coordinate is given (equivalent to WCS/CRPIX)
    reference_coordinate : 2d ndarray
        The data coordinate (x, y) as measured at the reference pixel (equivalent to WCS/CRVAL)

    Returns
    -------
    out : ndarray
        The data coordinates
    """
    return scale * (size - 1 * u.pix) / 2. + reference_coordinate - (reference_pixel - 1 * u.pix) * scale


@deprecated("0.8.0", alternative="sunpy.map.GenericMap.world_to_pixel")
def convert_data_to_pixel(x, y, scale, reference_pixel, reference_coordinate):
    """
    Calculate the pixel indices for a given data coordinate.

    Parameters
    ----------
    x, y : float
        Data coordinate in same units as reference coordinate
    scale : 2d ndarray
        The size of a pixel (dx,dy) in data coordinates (equivalent to WCS/CDELT)
    reference_pixel : 2d ndarray
        The reference pixel (x,y) at which the reference coordinate is given (equivalent to WCS/CRPIX)
    reference_coordinate : 2d ndarray
        The data coordinate (x, y) as measured at the reference pixel (equivalent to WCS/CRVAL)

    Returns
    -------
    out : ndarray
        The  pixel coordinates (x,y) at that data coordinate.
    """

    # TODO: Needs to check what coordinate system the data is given in
    cdelt = np.array(scale)
    crpix = np.array(reference_pixel)
    crval = np.array(reference_coordinate)
    # De-apply any tabular projections.
    # coord = inv_proj_tan(coord)

    # note that crpix[] counts pixels starting at 1
    pixelx = (x - crval[0]) / cdelt[0] + (crpix[0] - 1)
    pixely = (y - crval[1]) / cdelt[1] + (crpix[1] - 1)

    return pixelx, pixely


@deprecated("0.8.0", alternative="sunpy.coordinates")
def convert_hpc_hcc(x, y, dsun_meters=None, angle_units='arcsec', z=False):
    """
    Converts from Helioprojective-Cartesian (HPC) coordinates into
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
    z : Bool
        If true return the z coordinate as well.

    Returns
    -------
    out : ndarray
        The data coordinates (x,y,z) in heliocentric cartesian coordinates in meters.

    Notes
    -----
    Implements Eq. (15) of Thompson (2006), A&A, 449, 791.

    Examples
    --------
    >>> import sunpy.wcs
    >>> sunpy.wcs.convert_hpc_hcc(40.0, 32.0, z=True)
    (28876152.176423457, 23100922.07126697, 694524220.8157959)

    """
    c = np.array([_convert_angle_units(unit=angle_units),
                  _convert_angle_units(unit=angle_units)])

    cosx = np.cos(x * c[0])
    sinx = np.sin(x * c[0])
    cosy = np.cos(y * c[1])
    siny = np.sin(y * c[1])

    if dsun_meters is None:
        dsun_meters = sun.constants.au.si.value
    elif isinstance(dsun_meters, u.Quantity):
        dsun_meters = dsun_meters.si.value

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


@deprecated("0.8.0", alternative="sunpy.coordinates")
def convert_hcc_hpc(x, y, dsun_meters=None, angle_units='arcsec'):
    """
    Convert Heliocentric-Cartesian (HCC) to angular
    Helioprojective-Cartesian (HPC) coordinates (in degrees).

    Parameters
    ----------
    x, y : float (meters)
        Data coordinate in meters.
    dsun_meters : float
        Distance from the observer to the Sun in meters. Default is 1 AU.
    angle_units : str
        Units of the data coordinates (e.g. arcsec, arcmin, degrees). Default is arcsec.

    Returns
    -------
    out : ndarray
        The  data coordinates (x,y) in helioprojective cartesian coordinates in arcsec.

    Notes
    -----
    Implements Eq. (16) of Thompson (2006), A&A, 449, 791.

    Examples
    --------
    >>> import sunpy.wcs
    >>> sunpy.wcs.convert_hcc_hpc(28748691, 22998953)
    (39.823439773829705, 31.858751644835717)

    """

    # Calculate the z coordinate by assuming that it is on the surface of the Sun
    z = np.sqrt(rsun_meters ** 2 - x ** 2 - y ** 2)

    if dsun_meters is None:
        dsun_meters = sun.constants.au.si.value
    elif isinstance(dsun_meters, u.Quantity):
        dsun_meters = dsun_meters.si.value

    zeta = dsun_meters - z
    distance = np.sqrt(x**2 + y**2 + zeta**2)
    hpcx = np.rad2deg(np.arctan2(x, zeta))
    hpcy = np.rad2deg(np.arcsin(y / distance))

    if angle_units == 'arcsec':
        hpcx = 60 * 60 * hpcx
        hpcy = 60 * 60 * hpcy
    elif angle_units == 'arcmin':
        hpcx = 60 * hpcx
        hpcy = 60 * hpcy

    return hpcx, hpcy


@deprecated("0.8.0", alternative="sunpy.coordinates")
def convert_hcc_hg(x, y, z=None, b0_deg=0, l0_deg=0, radius=False):
    """
    Convert from Heliocentric-Cartesian (HCC) (given in meters) to
    Stonyhurst Heliographic coordinates (HG) given in degrees, with
    radial output in meters.

    Parameters
    ----------
    x, y : float (meters)
        Data coordinate in meters.
    z : float (meters)
        Data coordinate in meters.  If None, then the z-coordinate is assumed
        to be on the Sun.
    b0_deg : float (degrees)
        Tilt of the solar North rotational axis toward the observer
        (heliographic latitude of the observer). Usually given as SOLAR_B0,
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0_deg : float (degrees)
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    radius : Bool
        If true, forces the output to return a triple of (lon, lat, r). If
        false, return (lon, lat) only.

    Returns
    -------
    out : ndarray (degrees, meters)
        if radius is false, return the data coordinates (lon, lat).  If
        radius=True, return the data coordinates (lon, lat, r).  The quantities
        (lon, lat) are the heliographic coordinates in degrees.  The quantity
        'r' is the heliographic radius in meters.

    Notes
    -----
    Implements Eq. (12) of Thompson (2006), A&A, 449, 791.

    Examples
    --------
    >>> import sunpy.wcs
    >>> sunpy.wcs.convert_hcc_hg(230000.0,45000000.0,
    ...                          z=695508000.0 + 8000000.0, radius=True)
    (0.01873188196651189, 3.6599471896203317, 704945784.4146597)
    """
    if z is None:
        z = np.sqrt(rsun_meters**2 - x**2 - y**2)

    cosb = np.cos(np.deg2rad(b0_deg))
    sinb = np.sin(np.deg2rad(b0_deg))

    hecr = np.sqrt(x**2 + y**2 + z**2)
    hgln = np.arctan2(x, z * cosb - y * sinb) + np.deg2rad(l0_deg)
    hglt = np.arcsin((y * cosb + z * sinb) / hecr)

    if radius:
        return np.rad2deg(hgln), np.rad2deg(hglt), hecr
    else:
        return np.rad2deg(hgln), np.rad2deg(hglt)


@deprecated("0.8.0", alternative="sunpy.coordinates")
def convert_hg_hcc(hglon_deg, hglat_deg, b0_deg=0, l0_deg=0, occultation=False,
                   z=False, r=rsun_meters):
    """
    Convert from Stonyhurst Heliographic coordinates (given in degrees) to
    Heliocentric-Cartesian coordinates (given in meters).

    Parameters
    ----------
    hglon_deg, hglat_deg : float (degrees)
        Heliographic longitude and latitude in degrees.
    b0_deg : float (degrees)
        Tilt of the solar North rotational axis toward the observer
        (heliographic latitude of the observer). Usually given as SOLAR_B0,
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0_deg : float (degrees)
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    occultation : Bool
        If true set all points behind the Sun (e.g. not visible) to Nan.
    z : Bool
        If true return the z coordinate as well.
    r : float (meters)
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
    >>> import sunpy.wcs
    >>> sunpy.wcs.convert_hg_hcc(0.01873188196651189, 3.6599471896203317,
    ...                          r=704945784.41465974, z=True)
    (230000.0, 45000000.0, 703508000.0)
    """
    lon = np.deg2rad(hglon_deg)
    lat = np.deg2rad(hglat_deg)

    cosb = np.cos(np.deg2rad(b0_deg))
    sinb = np.sin(np.deg2rad(b0_deg))

    lon = lon - np.deg2rad(l0_deg)

    cosx = np.cos(lon)
    sinx = np.sin(lon)
    cosy = np.cos(lat)
    siny = np.sin(lat)

    # Perform the conversion.
    x = r * cosy * sinx
    y = r * (siny * cosb - cosy * cosx * sinb)
    zz = r * (siny * sinb + cosy * cosx * cosb)

    if occultation:
        x[zz < 0] = np.nan
        y[zz < 0] = np.nan

    if np.all(z == True):
        return x, y, zz
    else:
        return x, y


@deprecated("0.8.0", alternative="sunpy.coordinates")
def convert_hg_hpc(hglon_deg, hglat_deg, b0_deg=0, l0_deg=0, dsun_meters=None, angle_units='arcsec',
                   occultation=False):
    """
    Convert from Heliographic coordinates (HG) to Helioprojective-Cartesian
    (HPC).

    Parameters
    ----------
    hglon_deg, hglat_deg : float (degrees)
        Heliographic longitude and latitude in degrees.
    b0_deg : float (degrees)
        Tilt of the solar North rotational axis toward the observer
        (heliographic latitude of the observer). Usually given as SOLAR_B0,
        HGLT_OBS, or CRLT_OBS. Default is 0.
    l0_deg : float (degrees)
        Carrington longitude of central meridian as seen from Earth. Default is 0.
    occultation : Bool
        If true set all points behind the Sun (e.g. not visible) to Nan.
    dsun_meters : float (meters)
        Distance between the observer and the Sun.
    angle_units : str


    Returns
    -------
    out : ndarray (arcsec)
        The data coordinates (x,y) in Helioprojective-Cartesian coordinates.

    Notes
    -----
    Uses equations 11 and 16 in Thompson (2006), A&A, 449, 791-803.

    Examples
    --------
    >>> import sunpy.wcs
    >>> sunpy.wcs.convert_hg_hpc(34.0, 45.0, b0_deg=-7.064078, l0_deg=0.0)
    (380.056565603089, 743.7828128329002)
    """

    tempx, tempy = convert_hg_hcc(hglon_deg, hglat_deg, b0_deg=b0_deg, l0_deg=l0_deg, occultation=occultation)
    x, y = convert_hcc_hpc(tempx, tempy, dsun_meters=dsun_meters, angle_units=angle_units)
    return x, y

@deprecated("0.8.0", alternative="sunpy.coordinates")
def convert_hpc_hg(x, y, b0_deg=0, l0_deg=0, dsun_meters=None, angle_units='arcsec'):
    """
    Convert from Helioprojective-Cartesian (HPC) to Heliographic coordinates
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
        Distance between the observer and the Sun.
    angle_units : str
        Units used for input x and y. Default is arcsec.

    Returns
    -------
    out : ndarray (degrees)
        The  data coordinates (hglongitude, hglatitude) in Heliographic coordinates.

    Notes
    -----
    Uses equations 15 and 12 in Thompson (2006), A&A, 449, 791-803.

    Examples
    --------
    >>> import sunpy.wcs
    >>> sunpy.wcs.convert_hpc_hg(382, 748, b0_deg=-7.064078, l0_deg=0.0)
    (34.50465343991467, 45.44314327551818)
    """
    tempx, tempy = convert_hpc_hcc(x, y, dsun_meters=dsun_meters, angle_units=angle_units)
    lon, lat = convert_hcc_hg(tempx, tempy, b0_deg=b0_deg, l0_deg=l0_deg)
    return lon, lat


@deprecated("0.8.0", alternative="sunpy.map")
def proj_tan(x, y, force=False):
    """
    Applies the gnomonic (TAN) projection to intermediate relative
    coordinates. This function is not currently implemented!"""
    # if pixels are within 3 degrees of the Sun then skip the calculation unless
    # force is True. This applies to all sdo images so this function is just
    # here as a place holder for the future
    # TODO: write proj_tan function
    return x, y

@deprecated("0.8.0", alternative="sunpy.coordinates")
def convert_to_coord(x, y, from_coord, to_coord, b0_deg=0, l0_deg=0, dsun_meters=None, angle_units='arcsec'):
    """
    Apply a coordinate transform to coordinates. Right now can only do hpc
    to hcc to hg
    """

    if (from_coord == 'hcc') and (to_coord == 'hg'):
        rx, ry = convert_hcc_hg(x, y, b0_deg=b0_deg, l0_deg=l0_deg)
    elif (from_coord == 'hpc') and (to_coord == 'hg'):
        rx, ry = convert_hpc_hg(x, y, b0_deg=b0_deg, l0_deg=l0_deg, dsun_meters=dsun_meters, angle_units=angle_units)
    elif (from_coord == 'hg') and (to_coord == 'hcc'):
        rx, ry = convert_hg_hcc(x, y, b0_deg=b0_deg, l0_deg=l0_deg)
    elif (from_coord == 'hcc') and (to_coord == 'hpc'):
        rx, ry = convert_hcc_hpc(x, y, dsun_meters=dsun_meters, angle_units=angle_units)
    elif (from_coord == 'hg') and (to_coord == 'hpc'):
        rx, ry = convert_hg_hpc(x, y, b0_deg=b0_deg, l0_deg=l0_deg, dsun_meters=dsun_meters, angle_units=angle_units)
    elif (from_coord == 'hpc') and (to_coord == 'hcc'):
        rx, ry = convert_hpc_hcc(x, y, dsun_meters=dsun_meters, angle_units=angle_units)

    return rx, ry
