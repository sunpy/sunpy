import warnings

import numpy as np

import astropy.units as u
import astropy.wcs
from astropy.coordinates import SkyCoord

from sunpy.coordinates import frames
from sunpy.util import MetaDict
from sunpy.util.exceptions import SunpyUserWarning

__all__ = ['meta_keywords', 'make_fitswcs_header']


def meta_keywords():
    """
    Returns the metadata keywords that are used when creating a `sunpy.map.GenericMap`.

    Examples
    --------

    Returns a dictionary of all meta keywords that are used in a `sunpy.map.GenericMap` header:
        >>> import sunpy.map
        >>> sunpy.map.meta_keywords()
        {'cunit1': 'Units of the coordinate increments along naxis1 e.g. arcsec **required',
         'cunit2': 'Units of the coordinate increments along naxis2 e.g. arcsec **required',
         ...
    """
    return _map_meta_keywords


@u.quantity_input
def make_fitswcs_header(data, coordinate, reference_pixel: u.pix = None,
                        scale: u.arcsec/u.pix = None, **kwargs):
    """
    Function to create a FITS-WCS header from a coordinate object
    (`~astropy.coordinates.SkyCoord`) that is required to
    create a `~sunpy.map.GenericMap`.

    Parameters
    ----------
    data : `~numpy.ndarray`
        Array data of Map for which a header is required.
    coordinates : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseFrame`
        Coordinate object to get meta information for map header.
    reference_pixel :`~astropy.units.Quantity` of size 2, optional
        Reference pixel along each axis. These are expected to be Cartestian ordered, i.e
        the first index is the x axis, second index is the y axis.
        Defaults to the center of data array, ``(data.shape[1] + 1)/2., (data.shape[0] + 1)/2.)``.
    scale : `~astropy.units.Quantity` of size 2, optional
        Pixel scaling along x and y axis (i.e. the spatial scale of the pixels (dx, dy)). These are
        expected to be Cartestian ordered, i.e [dx, dy].
        Defaults to ``([1., 1.] arcsec/pixel)``.
    **kwargs:
        Additional arguments that will be put into the metadict header if they
        are in the list returned by `~sunpy.map.meta_keywords`. Additional
        keyword arguments for the instrument meta can also be given in the
        following forms which will be translated to fits standard:

        | ``instrument``
        | ``telescope``
        | ``observatory``
        | ``wavelength``
        | ``exposure``

    Returns
    -------
    `~sunpy.util.MetaDict`
        The header information required for making a `sunpy.map.GenericMap`.

    Examples
    --------
    >>> import sunpy.map
    >>> from sunpy.coordinates import frames
    >>> from astropy.coordinates import SkyCoord
    >>> from astropy import units as u
    >>> import numpy as np

    >>> data = np.random.rand(1024, 1024)
    >>> my_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime="2017-08-01",
    ...                     observer = 'earth', frame=frames.Helioprojective)
    >>> my_header = sunpy.map.make_fitswcs_header(data, my_coord)
    >>> my_map = sunpy.map.Map(data, my_header)
    """

    if not isinstance(coordinate, (SkyCoord, frames.BaseCoordinateFrame)):
        raise ValueError("coordinate needs to be a coordinate frame or an SkyCoord instance.")

    if isinstance(coordinate, SkyCoord):
        coordinate = coordinate.frame

    if coordinate.obstime is None:
        raise ValueError("The coordinate needs an observation time, `obstime`.")

    if isinstance(coordinate, frames.Heliocentric):
        raise ValueError("This function does not currently support heliocentric coordinates.")

    meta_wcs = _get_wcs_meta(coordinate)

    if hasattr(coordinate, "observer") and isinstance(coordinate.observer, frames.BaseCoordinateFrame):
        meta_observer = _get_observer_meta(coordinate)
        meta_wcs.update(meta_observer)

    meta_instrument = _get_instrument_meta(**kwargs)
    meta_wcs.update(meta_instrument)

    if reference_pixel is None:
        reference_pixel = u.Quantity([(data.shape[1] + 1)/2.*u.pixel, (data.shape[0] + 1)/2.*u.pixel])
    if scale is None:
        scale = u.Quantity([1.0*u.arcsec, 1.0*u.arcsec])

    meta_wcs['crval1'], meta_wcs['crval2'] = (coordinate.spherical.lat.to_value(meta_wcs['cunit1']),
                                              coordinate.spherical.lon.to_value(meta_wcs['cunit2']))

    meta_wcs['crpix1'], meta_wcs['crpix2'] = (reference_pixel[0].to_value(u.pixel),
                                              reference_pixel[1].to_value(u.pixel))

    meta_wcs['cdelt1'], meta_wcs['cdelt2'] = (scale[0]*meta_wcs['cunit1']/u.pixel,
                                              scale[1]*meta_wcs['cunit2']/u.pixel)

    meta_dict = MetaDict(meta_wcs)

    for key in kwargs:
        if key in _map_meta_keywords:
            meta_dict[key] = kwargs[key]

    return meta_dict


def _get_wcs_meta(coordinate):
    """
    Function to get WCS meta from the SkyCoord using
    `astropy.wcs.utils.celestial_frame_to_wcs`

    Parameters
    ----------
    coordinate : ~`astropy.coordinates.BaseFrame`

    Returns
    -------
    `dict`
        Containing the WCS meta information
            * ctype1, ctype2
            * cunit1, cunit2
            * date_obs
    """

    coord_meta = {}

    skycoord_wcs = astropy.wcs.utils.celestial_frame_to_wcs(coordinate)

    cunit1, cunit2 = skycoord_wcs.wcs.cunit
    coord_meta = dict(skycoord_wcs.to_header())
    coord_meta['cunit1'], coord_meta['cunit2'] = cunit1, cunit2

    return coord_meta


def _get_observer_meta(coordinate):
    """
    Function to get observer meta from coordinate frame.

    Parameters
    ----------
    coordinate : ~`astropy.coordinates.BaseFrame`

    Returns
    -------
    `dict`
        Containing the WCS meta information
            * hgln_obs, hglt_obs
            * dsun_obs
            * rsun_obs
    """
    coord_meta = {}

    coord_meta['hgln_obs'] = coordinate.observer.lon.to_value(u.deg)
    coord_meta['hglt_obs'] = coordinate.observer.lat.to_value(u.deg)
    coord_meta['dsun_obs'] = coordinate.observer.radius.to_value(u.m)
    coord_meta['rsun_ref'] = coordinate.rsun.to_value()
    coord_meta['rsun_obs'] = np.arctan(coordinate.rsun / coordinate.observer.radius).to_value(u.arcsec)

    return coord_meta


def _get_instrument_meta(**kwargs):
    """
    Function to correctly name keywords from kwargs
    """
    coord = {}

    conversion = {'instrument': 'instrume', 'telescope': 'telescop',
                  'observatory': 'obsrvtry', 'wavelength': 'wavelnth', 'exposure': 'exptime'}

    kwargs_temp = kwargs.copy()
    for key in kwargs_temp:
        if key in conversion:
            coord[conversion[key]] = kwargs.pop(key)

    return coord


_map_meta_keywords = {
    'cunit1':
    'Units of the coordinate increments along naxis1 e.g. arcsec **required',
    'cunit2':
    'Units of the coordinate increments along naxis2 e.g. arcsec **required',
    'crval1':
    'Coordinate value at reference point on naxis1 **required',
    'crval2':
    'Coordinate value at reference point on naxis2 **required',
    'cdelt1':
    'Spatial scale of pixels for naxis1, i.e. coordinate increment at reference point',
    'cdelt2':
    'Spatial scale of pixels for naxis2, i.e. coordinate increment at reference point',
    'crpix1':
    'Pixel coordinate at reference point naxis1',
    'crpix2':
    'Pixel coordinate at reference point naxis2',
    'ctype1':
    'Coordinate type projection along naxis1 of data e.g. HPLT-TAN',
    'ctype2':
    'Coordinate type projection along naxis2 of data e.g. HPLN-TAN',
    'hgln_obs':
    'Heliographic longitude of observation',
    'hglt_obs':
    'Heliographic latitude of observation',
    'dsun_obs':
    'distance to Sun from observation in metres',
    'rsun_obs':
    'radius of Sun in meters from observation',
    'date-obs':
    'date of observation e.g. 2013-10-28 00:00',
    'date_obs':
    'date of observation e.g. 2013-10-28 00:00',
    'rsun_ref':
    'reference radius of Sun in meters',
    'solar_r':
    'radius of Sun in meters from observation',
    'radius':
    'radius of Sun in meters from observation',
    'crln_obs':
    'Carrington longitude of observation',
    'crlt_obs':
    'Heliographic latitude of observation',
    'solar_b0':
    'Solar B0 angle',
    'detector':
    'name of detector e.g. AIA',
    'exptime':
    'exposure time of observation, in seconds e.g 2',
    'instrume':
    'name of instrument',
    'wavelnth':
    'wavelength of observation',
    'waveunit':
    'unit for which observation is taken e.g. angstom',
    'obsrvtry':
    'name of observatory of observation',
    'telescop':
    'name of telescope of observation',
    'lvl_num':
    'FITS processing level',
    'crota2':
    'Rotation of the horizontal and vertical axes in degrees',
    'PC1_1':
    'Matrix element PCi_j describing the rotation required to align solar North with the top of the image.',
    'PC1_2':
    'Matrix element PCi_j describing the rotation required to align solar North with the top of the image.',
    'PC2_1':
    'Matrix element PCi_j describing the rotation required to align solar North with the top of the image.',
    'PC2_2':
    'Matrix element PCi_j describing the rotation required to align solar North with the top of the image.',
    'CD1_1':
    'Matrix element CDi_j describing the rotation required to align solar North with the top of the image.',
    'CD1_2':
    'Matrix element CDi_j describing the rotation required to align solar North with the top of the image.',
    'CD2_1':
    'Matrix element CDi_j describing the rotation required to align solar North with the top of the image.',
    'CD2_2':
    'Matrix element CDi_j describing the rotation required to align solar North with the top of the image.'
}
