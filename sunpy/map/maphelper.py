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

__all__ = ['meta_keywords', 'make_fitswcs_header']


def meta_keywords(*key):
    """
    Prints the accepted metadata keywords that are used when
    creating a `sunpy.map.Map`.

    Examples
    --------

    * Print all accepted meta keywords that are used in a `sunpy.map.GenericMap` header:
        >>> from sunpy import map
        >>> map.meta_keywords()
        cunit1  :  Units of the coordinate increments along naxis1 e.g. arcsec **required
        cunit2  :  Units of the coordinate increments along naxis2 e.g. arcsec **required
        crval1  :  Coordinate value at reference point on naxis1 **required
        crval2  :  Coordinate value at reference point on naxis2 **required
        ...

    * Print information about certain meta keyword:
        >>> map.meta_keywords('ctype1')
        ctype1  :  Coordinate type projection along naxis1 of data e.g. HPLT-TAN

    """
    if key:
        for k in key:
            if k in map_meta_keywords.keys():
                print(k, " : ", acceptable_meta[k])
            else:
                print(k, " is not a meta keyword required for sunpy.map.Map \n The map meta keywords are : \n", list(map_meta_keywords.keys()))
        return

    for k in map_meta_keywords:
        print(k, " : ", map_meta_keywords[k])


@u.quantity_input
def make_fitswcs_header(data, coordinate, reference_pixel: u.pix = None, scale: u.arcsec/u.pix = None, **kwargs):
    """
    Function to create a `sunpy.util.meta.MetaDict` from a coordinate object (`astropy.coordinates.SkyCoord` or coordinate frame)
    that is required to create a `sunpy.map.Map`

    Parameters
    ----------
    data : `~numpy.ndarray`
        Array data of Map for which a header is required.
    coordinates : `~astropy.coordinates.SkyCoord`
        Coordinate object to get meta information for map header. 
    reference_pixel :`~astropy.units.Quantity` of size 2, optional
        Reference pixel along each axis. These are expected to be Cartestian ordered, i.e
        the first index is the x axis, second index is the y axis.
        Defaults to the center of data array, ``(data.shape[0] + 1)/2., (data.shape[1] + 1)/2.)``.
    scale : `~astropy.units.Quantity` of size 2, optional
        Pixel scaling along x and y axis (i.e. the spatial scale of the pixels (dx, dy)). These are 
        expected to be Cartestian ordered, i.e [dx, dy].
        Defaults to ``([1., 1.] arcsec)``.
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
    >>> my_header = map.make_fits_wcs(data, my_coord)
    >>> my_map = map.Map(data, my_header)
    """

    if not isinstance(coordinate, SkyCoord):
        raise ValueError("coordinate needs to be a SkyCoord instance.")

    if coordinate.obstime is None:
      raise ValueError("The coordinate needs an observation time, `obstime`.")

    coordinate = coordinate.transform_to(frames.Helioprojective)

    meta_wcs = _get_wcs_meta(coordinate)

    meta_observer = _get_observer_meta(coordinate)

    meta_instrument = _get_instrument_meta(**kwargs)

    meta_wcs.update(meta_observer)
    meta_wcs.update(meta_instrument)

    if reference_pixel is None:
        reference_pixel = u.Quantity([(data.shape[1] + 1)/2.*u.pixel, (data.shape[0] + 1)/2.*u.pixel])
    if scale is None:
        scale = u.Quantity([1.0*u.arcsec, 1.0*u.arcsec])

    meta_wcs['crval1'], meta_wcs['crval2'] = coordinate.Tx.to_value(meta_wcs['cunit1']), coordinate.Ty.to_value(meta_wcs['cunit2'])    

    meta_wcs['crpix1'], meta_wcs['crpix2'] = reference_pixel[0].to_value(u.pixel), reference_pixel[1].to_value(u.pixel)

    meta_wcs['cdelt1'], meta_wcs['cdelt2'] = scale[0]*meta_wcs['cunit1']/u.pixel, scale[1]*meta_wcs['cunit2']/u.pixel

    meta_dict = MetaDict(meta_wcs)

    for key in kwargs:
        if key in map_meta_keywords:
            meta_dict[key] = kwargs[key]

    return meta_dict


def _get_wcs_meta(coordinate):
    """
    Function to get WCS meta from the SkyCoord using 
    `astropy.wcs.utils.celestial_frame_to_wcs`
    
    Parameters
    ----------
    coordinate : ~`astropy.coordinates.SkyCoord`

    Returns
    -------
    `dict`
        Containing the WCS meta information
            * ctype1, ctype2
            * cunit1, cunit2
            * date_obs
    """

    coord_meta = {}

    skycoord_wcs = astropy.wcs.utils.celestial_frame_to_wcs(coordinate.frame)

    cunit1, cunit2 = skycoord_wcs.wcs.cunit
    coord_meta = dict(skycoord_wcs.to_header())
    coord_meta['cunit1'], coord_meta['cunit2'] = cunit1, cunit2  

    return coord_meta

def _get_observer_meta(coordinate):
    """
    Function to get observer meta from the SkyCoord using 
    `astropy.wcs.utils.celestial_frame_to_wcs`
    
    Parameters
    ----------
    coordinate : ~`astropy.coordinates.SkyCoord`

    Returns
    -------
    `dict`
        Containing the WCS meta information
            * hgln_obs, hglt_obs
            * dsun_obs
            * rsun_obs
    """    
    coord_meta = {}

    skycoord_wcs = astropy.wcs.utils.celestial_frame_to_wcs(coordinate.frame)

    coord_meta['hgln_obs'], coord_meta['hglt_obs'] = skycoord_wcs.heliographic_observer.lon.to_value(), skycoord_wcs.heliographic_observer.lat.to_value()
    coord_meta['dsun_obs'] = skycoord_wcs.heliographic_observer.radius.to(u.m).to_value()
    coord_meta['rsun_ref'] = skycoord_wcs.rsun.to(u.m).to_value()
    coord_meta['rsun_obs'] = ((skycoord_wcs.rsun/skycoord_wcs.heliographic_observer.radius).decompose() * u.radian).to(u.arcsec).to_value()

    return coord_meta

def _get_instrument_meta(**kwargs):
    """
    Function to correctly name keywords from kwargs
    """
    coord = {}

    conversion = {'instrument':'instrume', 'telescope':'telescop', 'observatory':'obsrvtry', 'wavelength':'wavelnth', 'exposure':'exptime'}

    for key in kwargs:
        if key in conversion:
            coord[conversion[key]] = kwargs[key]

    return coord


map_meta_keywords = {'cunit1': 'Units of the coordinate increments along naxis1 e.g. arcsec **required',
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
