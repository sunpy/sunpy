"""
This submodule provides utility functions to act on `sunpy.map.Map` instances.
"""
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

__all__ = ['MapHelper']


class MapHelper:
    """
    A helper class for knowing what is required in a `sunpy.map.Map` header.
    A `sunpy.map.Map` header consists of a dictionary that contains meta information
    about the map data. This meta is generally found in the FITS header of an image
    file. This helper is defined to assist the user in finding out what meta data
    keywords are accepted in a map header, and can produce a header MetaDict from a 
    coordinate object, or through default meta that can be then passed into `Map` 
    to create a Map object.

    `map_helper.acceptable_meta_keywords()` prints the accepted meta that is used in 
    creating a sunpy.map.Map

    `map_helper.make_header(data)` will create a map header for the numpy.ndarray data, 
    either with default meta keywords, or a SkyCoord or frame can be passed to extract
    reference meta

    Examples
    --------
    Printing meta keywords:

        * print the accepeted meta keywords that are used in a sunpy.map.Map header
          >>> from sunpy import map
          >>> map.MapHelper.acceptable_meta_keywords() 
          cunit1  :  Units of the coordinate increments along naxis1 e.g. arcsec **required
          ...
          CD2_2  :  Matrix element CDi_j describing the rotation required to align solar North with the top of the image.

        * check is a meta keyword is acceptable and print information about it
          >>> map.MapHelper.acceptable_meta_keywords('exptime')
          exptime is an accepted keyword  :  exposure time of observation, in seconds e.g 2

          >>> map.MapHelper.acceptable_meta_keywords('keyword_test')
          keyword_test  is not in the accepted meta data required for sunpy.map.Map 
    
    Creating headers from map_helper:

        * Create a default header with meta data assumed to be Earth-based observer at current time and arcsec units
          >>> from sunpy import map
          >>> from sunpy.coordinates import frames
          >>> from astropy.coordinates import SkyCoord
          >>> from astropy import units as u
          >>> import numpy as np

          >>> data = np.random.rand(1024, 1024)
          >>> header = map.MapHelper.make_header(data)
          >>> my_map = map.Map(data, header)
        
        * Create a header from a reference SkyCoord object
          >>> my_coord = SkyCoord(100*u.arcsec, 100*u.arcsec, obstime = '2013-10-28 01:00', frame = frames.Helioprojective)
          >>> header = map.MapHelper.make_header(data, my_coord)
          >>> my_map = map.Map(data, header)



    """
    def acceptable_meta_keywords(self, *key):
        """
        prints the accepted metadata keywords that are used when
        creating a `sunpy.map.Map`.
        """
        if key:
            for k in key:
                if k in acceptable_meta.keys():
                    print(k, ' : ', acceptable_meta[k])
                else:
                    print(k, ' is not in the accepted meta data required for sunpy.map.Map \n The accepted keywords are : \n', list(acceptable_meta.keys()))
            return

        for k in acceptable_meta:
            print(k, ' : ', acceptable_meta[k])

    def make_header(self, data, *coordinate, crval=None, crpix=None, cdelt=None, **kwargs):
        """
        Function to create a MetaDict from a coordinate object (SkyCoord or coordinate frame)
        that is required to create a `sunpy.map.Map`

        Parameters
        ----------
        data : `~numpy.ndarray`
            array data of Map for which a header is requried

        coordinates : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseFrame`, optional if not set default parameters are set
            coordinate object to get meta information for map header

        crval : array or list of size of 2, default (0, 0)
            coordinate of reference pixel along each axis
            first index is crval1, second index is crval2

        crpix : array or list of size of 2, default center of data array ((data.shape[0] + 1)/2., (data.shape[1] + 1)/2.)
            reference pixel along each axis
            first index is crpix1, second index is crpix2

        cdelt : array or list of size of 2, default (0, 0)
            pixel scaling along each axis (i.e. the spatial scale of the pixels)
            first index is cdelt1, second index is cdelt2

        **kwargs
            additional arguments that will be put into metadict header. These must be in the
            list of acceptable meta keywords

        Returns
        -------
        meta_dict : `~sunpy.util.MetaDict`
            a MetaDict that contains the header information required for making a Map

        Notes
        -----
        The default parameters that are set if no coordinate is given is as follows:
            * ctype1, ctype2 - HPLN-TAN, HPLT-TAN
            * cunit1, cunit2 - arcsec, arcsec
            * hgln_obs, hglt_obs - 0, 0
            * dsun_obs - reference Sun-Earth distance
            * rsun_obs - reference radius of Sun as seen from Earth
            * t_obs, date_obs - present time

        Examples
        --------
        >>> from sunpy import map
        >>> from sunpy.coordinates import frames
        >>> from astropy.coordinates import SkyCoord
        >>> from astropy import units as u
        >>> import numpy as np

        >>> data = np.random.rand(1024, 1024)
        >>> my_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime="2017-08-01", observer = 'earth', frame=frames.Helioprojective)
        >>> my_header = map.MapHelper.make_header(data, my_coord, crval = [0,0], cdelt = [2,2])
        >>> my_map = map.Map(data, my_header)

        """

        if coordinate:
            if not isinstance(coordinate[0], (SkyCoord, frames.BaseCoordinateFrame)):
                raise ValueError('coordinate needs to be an Astropy coordinate frame or SkyCoord instance.')

            meta = _get_meta_from_coordinate(coordinate[0])

        else:
            print('Using default meta assuming Earth based observer at current time with `arcsec` cunit type')
            meta = _default_meta()

        # if no crval options set, default to 0, 0
        if crval is None:
            meta['crval1'], meta['crval2'] = 0, 0
        else:
            meta['crval1'], meta['crval2'] = crval[0], crval[1]

        # if no crpix options set, default to 0, 0
        if crpix is None:
            meta['crpix1'], meta['crpix2'] = (data.shape[0] + 1)/2., (data.shape[1] + 1)/2.
        else:
            meta['crpix1'], meta['crpix2'] = crpix[0], crpix[1]

        # if no cdelt options set, default to 1, 1
        if cdelt is None:
            meta['cdelt1'], meta['cdelt2'] = 1, 1
        else:
            meta['cdelt1'], meta['cdelt2'] = cdelt[0], cdelt[1]

        meta_dict = MetaDict({
                    'ctype1': meta.get('ctype1'),
                    'ctype2': meta.get('ctype2'),
                    'cunit1': meta.get('cunit1'),
                    'cunit2': meta.get('cunit1'),
                    'crpix1': meta.get('crpix1'),
                    'crpix2': meta.get('crpix2'),
                    'cdelt1': meta.get('cdelt1'),
                    'cdelt2': meta.get('cdelt2'),
                    'crval1': meta.get('crval1'),
                    'crval2': meta.get('crval2'),
                    'hgln_obs': meta.get('hgln_obs'),
                    'hglt_obs': meta.get('hglt_obs'),
                    'date_obs': meta.get('date_obs'),
                    't_obs': meta.get('t_obs'),
                    'dsun_obs': meta.get('dsun_obs'),
                    'rsun_obs': meta.get('rsun_obs'),
                    'rsun_ref': meta.get('rsun_ref')
        })

        for key in kwargs:
            if key in acceptable_meta:
                meta_dict[key] = kwargs[key]

        return meta_dict


def _default_meta():
    """
    Creates default meta data if none is specified
    """
    coord_meta = {}

    coord_meta['ctype1'], coord_meta['ctype2'] = 'HPLN-TAN', 'HPLT-TAN'
    coord_meta['cunit1'], coord_meta['cunit2'] = 'arcsec', 'arcsec'
    coord_meta['hgln_obs'], coord_meta['hglt_obs'] = 0.0, 0.0
    coord_meta['dsun_obs'] = con.au.to(u.m).value
    coord_meta['rsun_ref'] = con.radius.to(u.m).value
    coord_meta['rsun_obs'] = ((con.radius / con.au).decompose() * u.radian).to(u.arcsec).value
    coord_meta['t_obs'] = parse_time('now').iso
    coord_meta['date_obs'] = parse_time('now').iso

    return coord_meta


def _get_meta_from_coordinate(coordinate):
    """
    Function to get WCS meta data from a SkyCoord or frame using
    `astropy.wcs.utils.celestial_frame_to_wcs`

    Parameters
    ----------
    coordinate : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseFrame`
        coordinate object to pull meta data from

    Returns
    -------
    coord_meta : dict
        containing the meta information required for `sunpy.map.Map` header
            * ctype1, ctype2
            * cunit1, cunit2
            * hgln_obs, hglt_obs
            * dsun_obs
            * rsun_obs
            * t_obs, date_obs

    """
    coord_meta = {}

    if not isinstance(coordinate, (SkyCoord, frames.BaseCoordinateFrame)):
        raise ValueError('takes an Astropy coordinate frame or SkyCoord instance.')

    if isinstance(coordinate, SkyCoord):
        if coordinate.obstime is None or coordinate.frame is None:
            raise ValueError('SkyCoord needs an observation time and a frame')
        skycoord_wcs = astropy.wcs.utils.celestial_frame_to_wcs(coordinate.frame)
    elif isinstance(coordinate, frames.SunPyBaseCoordinateFrame):
        if coordinate.obstime is None:
            raise ValueError('frame needs an observation time')
        skycoord_wcs = astropy.wcs.utils.celestial_frame_to_wcs(coordinate)

    coord_meta['ctype1'], coord_meta['ctype2'] = skycoord_wcs.wcs.ctype
    coord_meta['cunit1'], coord_meta['cunit2'] = skycoord_wcs.wcs.cunit
    coord_meta['hgln_obs'], coord_meta['hglt_obs'] = skycoord_wcs.heliographic_observer.lon, skycoord_wcs.heliographic_observer.lat
    coord_meta['dsun_obs'] = skycoord_wcs.heliographic_observer.radius.to(u.m)
    coord_meta['rsun_ref'] = skycoord_wcs.rsun.to(u.m)
    coord_meta['rsun_obs'] = ((skycoord_wcs.rsun.to(u.AU)/skycoord_wcs.heliographic_observer.radius.to(u.AU)).decompose() * u.radian).to(u.arcsec).value
    coord_meta['t_obs'], coord_meta['date_obs'] = skycoord_wcs.heliographic_observer.obstime, skycoord_wcs.heliographic_observer.obstime

    return coord_meta

MapHelper = MapHelper()

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
