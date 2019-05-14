from sunpy import map
from sunpy.util import MetaDict
from astropy.coordinates import SkyCoord
from astropy import units as u
from sunpy.sun import constants as con
from sunpy.time import parse_time
import numpy as np
import random
import matplotlib.pyplot as plt
import astropy
from sunpy.coordinates import frames
import warnings
from sunpy.util.exceptions import SunpyUserWarning

__all__ = ['meta_keywords', 'make_header', 'make_header_from_skycoord']


def meta_keywords(*key):
    """
    Prints the accepted metadata keywords that are used when
    creating a `sunpy.map.Map`.
    """
    if key:
        for k in key:
            if k in acceptable_meta.keys():
                print(k, " : ", acceptable_meta[k])
            else:
                print(k, " is not in the accepted meta data required for sunpy.map.Map \n The accepted keywords are : \n", list(acceptable_meta.keys()))
        return

    for k in acceptable_meta:
        print(k, " : ", acceptable_meta[k])


def make_header(data):
	return MetaDict({'cdelt1':1, 'cdelt2':1})


@u.quantity_input
def make_header_from_skycoord(data, coordinate, crpix=None, cdelt=None, **kwargs):
    """
    Function to create a `sunpy.util.meta.MetaDict` from a coordinate object (`astropy.coordinates.SkyCoord` or coordinate frame)
    that is required to create a `sunpy.map.Map`

    Parameters
    ----------
    data : `~numpy.ndarray`
        Array data of Map for which a header is required.
    coordinates : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseFrame`, optional
        Coordinate object to get meta information for map header.
         If not set, default parameters are used.
    crpix : `numpy.array` or `list` of size of 2, optional
        Reference pixel along each axis.
        first index is ``crpix1``, second index is ``crpix2``.
        Defaults to the venter of data array, ``(data.shape[0] + 1)/2., (data.shape[1] + 1)/2.)``.
    cdelt : `numpy.array` or `list` of size of 2, optional
        Pixel scaling along each axis (i.e. the spatial scale of the pixels)
        first index is ``cdelt1``, second index is ``cdelt2``.
        Defaults to ``(0. 0)``.
    **kwargs:
        Additional arguments that will be put into metadict header. These must be in the
        list of acceptable meta keywords.

    Returns
    -------
    `~sunpy.util.MetaDict`
        A `sunpy.util.MetaDict` that contains the header information required for making a `sunpy.map.Map`

    Examples
    --------
    >>> from sunpy import map
    >>> from sunpy.coordinates import frames
    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> import numpy as np

    >>> data = np.random.rand(1024, 1024)
    >>> my_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime="2017-08-01", observer = 'earth', frame=frames.Helioprojective)
    >>> my_header = map.make_header_from_skycoord(data, my_coord, crval = [0,0], cdelt = [2,2])
    >>> my_map = map.Map(data, my_header)
    """

    if not isinstance(coordinate, (SkyCoord, frames.BaseCoordinateFrame)):
        raise ValueError("coordinate needs to be an Astropy coordinate frame or SkyCoord instance.")
    meta = _get_meta_from_coordinate(coordinate)

    meta['crval1'], meta['crval2'] = coordinate.Tx.value, coordinate.Ty.value
    # if no crpix options set, default to 0, 0
    if crpix is None:
        meta['crpix1'], meta['crpix2'] = (data.shape[0] + 1)/2., (data.shape[1] + 1)/2.
    else:
        meta['crpix1'], meta['crpix2'] = crpix[0], crpix[1]

    if cdelt is not None and isinstance(cdelt[0] and cdelt[1],u.Quantity):
        meta['cdelt1'], meta['cdelt2'] = cdelt[0].to(u.deg), cdelt[1].to(u.deg)

    meta_dict = MetaDict(meta)

    for key in kwargs:
        if key in acceptable_meta:
            meta_dict[key] = kwargs[key]

    return meta_dict  
    
def _get_meta_from_coordinate(coordinate):
    """
    Function to get WCS meta data from an `astropy.coordinates.SkyCoord` or frame using
    `astropy.wcs.utils.celestial_frame_to_wcs`.

    Parameters
    ----------
    coordinate : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseFrame`
        Coordinate object to pull metadata from.

    Returns
    -------
    `dict`
        Containing the meta information required for `sunpy.map.Map` header:
            * ctype1, ctype2
            * cunit1, cunit2
            * hgln_obs, hglt_obs
            * dsun_obs
            * rsun_obs
            * t_obs, date_obs
    """

    coord_meta = {}

    if not isinstance(coordinate, (SkyCoord, frames.BaseCoordinateFrame)):
        raise ValueError("Input must be an Astropy coordinate frame or SkyCoord instance.")

    if isinstance(coordinate, SkyCoord):
        if coordinate.obstime is None or coordinate.frame is None:
            raise ValueError("SkyCoord needs an observation time and a frame.")
        skycoord_wcs = astropy.wcs.utils.celestial_frame_to_wcs(coordinate.frame)
    elif isinstance(coordinate, frames.SunPyBaseCoordinateFrame):
        if coordinate.obstime is None:
            raise ValueError("Frame needs an observation time.")
        skycoord_wcs = astropy.wcs.utils.celestial_frame_to_wcs(coordinate)

    coord_meta = dict(skycoord_wcs.to_header())

    coord_meta['hgln_obs'], coord_meta['hglt_obs'] = skycoord_wcs.heliographic_observer.lon, skycoord_wcs.heliographic_observer.lat
    coord_meta['dsun_obs'] = skycoord_wcs.heliographic_observer.radius.to(u.m)
    coord_meta['rsun_ref'] = skycoord_wcs.rsun.to(u.m)
    coord_meta['rsun_obs'] = ((skycoord_wcs.rsun.to(u.AU)/skycoord_wcs.heliographic_observer.radius.to(u.AU)).decompose() * u.radian).to(u.arcsec).value

    return coord_meta


acceptable_meta = {'cunit1': 'Units of the coordinate increments along naxis1 e.g. arcsec **required',
                   'cunit2': 'Units of the coordinate increments along naxis2 e.g. arcsec **required',
                   'crval1': 'Coordinate value at reference point on naxis1 **required',
                   'crval2': 'Coordinate value at reference point on naxis2 **required',
                   'cdelt1': 'Spatial scale of pixels for naxis1, i.e. coordinate increment at reference point',
                   'cdelt2': 'Spatial scale of pixels for naxis2, i.e. coordinate increment at reference point',
                   'crpix1': 'Pixel coordinate at reference point naxis1',
                   'crpix2': 'Pixel coordinate at reference point naxis2',
                   'ctype1': 'Coordinate type projection along naxis1 of data e.g. HPLT-TAN',
                   'ctype2': 'Coordinate type projection along naxis2 of data e.g. HPLN-TAN',
                   'hgln_obs': 'Heliographic longitude of observation',
                   'hglt_obs': 'Heliographic latitude of observation',
                   'dsun_obs': 'distance to Sun from observation in metres',
                   'rsun_obs': 'radius of Sun in meters from observation',
                   'date-obs': 'date of observation e.g. 2013-10-28 00:00',
                   'date_obs': 'date of observation e.g. 2013-10-28 00:00',
                   'rsun_ref': 'reference radius of Sun in meters',
                   'solar_r': 'radius of Sun in meters from observation',
                   'radius': 'radius of Sun in meters from observation',
                   'crln_obs': 'Carrington longitude of observation',
                   'crlt_obs': 'Heliographic latitude of observation',
                   'solar_b0': 'Solar B0 angle',
                   'detector': 'name of detector e.g. AIA',
                   'exptime': 'exposure time of observation, in seconds e.g 2',
                   'instrume': 'name of instrument',
                   'wavelnth': 'wavelength of observation',
                   'waveunit': 'unit for which observation is taken e.g. angstom',
                   'obsrvtry': 'name of observatory of observation',
                   'telescop': 'name of telescope of observation',
                   'lvl_num': 'FITS processing level',
                   'crota2': 'Rotation of the horizontal and vertical axes in degrees',
                   'PC1_1': 'Matrix element PCi_j describing the rotation required to align solar North with the top of the image.',
                   'PC1_2': 'Matrix element PCi_j describing the rotation required to align solar North with the top of the image.',
                   'PC2_1': 'Matrix element PCi_j describing the rotation required to align solar North with the top of the image.',
                   'PC2_2': 'Matrix element PCi_j describing the rotation required to align solar North with the top of the image.',
                   'CD1_1': 'Matrix element CDi_j describing the rotation required to align solar North with the top of the image.',
                   'CD1_2': 'Matrix element CDi_j describing the rotation required to align solar North with the top of the image.',
                   'CD2_1': 'Matrix element CDi_j describing the rotation required to align solar North with the top of the image.',
                   'CD2_2': 'Matrix element CDi_j describing the rotation required to align solar North with the top of the image.'}
