
import numpy as np

import astropy.units as u
import astropy.wcs
from astropy.coordinates import SkyCoord

from sunpy.coordinates import frames, sun
from sunpy.util import MetaDict

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


@u.quantity_input(equivalencies=u.spectral())
def make_fitswcs_header(data, coordinate,
                        reference_pixel: u.pix = None,
                        scale: u.arcsec/u.pix = None,
                        rotation_angle: u.deg = None,
                        rotation_matrix=None,
                        instrument=None,
                        telescope=None,
                        observatory=None,
                        wavelength: u.angstrom = None,
                        exposure: u.s = None,
                        projection_code="TAN"):
    """
    Function to create a FITS-WCS header from a coordinate object
    (`~astropy.coordinates.SkyCoord`) that is required to
    create a `~sunpy.map.GenericMap`.

    Parameters
    ----------
    data : `~numpy.ndarray` or `tuple`
        Array data of Map for which a header is required, or the shape of the
        data array (in numpy order, i.e. ``(y_size, x_size)``).
    coordinate : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseFrame`
        The coordinate of the reference pixel.
    reference_pixel :`~astropy.units.Quantity` of size 2, optional
        Reference pixel along each axis. These are expected to be Cartestian ordered, i.e
        the first index is the x axis, second index is the y axis. Defaults to
        the center of data array, ``(data.shape[1] - 1)/2., (data.shape[0] - 1)/2.)``,
        this argument is zero indexed (Python convention) not 1 indexed (FITS
        convention).
    scale : `~astropy.units.Quantity` of size 2, optional
        Pixel scaling along x and y axis (i.e. the spatial scale of the pixels (dx, dy)). These are
        expected to be Cartestian ordered, i.e [dx, dy].
        Defaults to ``([1., 1.] arcsec/pixel)``.
    rotation_angle : `~astropy.units.Quantity`, optional
        Coordinate system rotation angle, will be converted to a rotation
        matrix and stored in the ``PCi_j`` matrix. Can not be specified with
        ``rotation_matrix``.
    rotation_matrix : `~numpy.ndarray` of dimensions 2x2, optional
        Matrix describing the rotation required to align solar North with
        the top of the image in FITS ``PCi_j`` convention. Can not be specified
        with ``rotation_angle``.
    instrument : `~str`, optional
        Name of the instrument of the observation.
    telescope : `~str`, optional
        Name of the telescope of the observation.
    observatory : `~str`, optional
        Name of the observatory of the observation.
    wavelength : `~astropy.units.Quantity`, optional
        Wavelength of the observation as an astropy quanitity, e.g. 171*u.angstrom.
        From this keyword, the meta keywords ``wavelnth`` and ``waveunit`` will be populated.
    exposure : `~astropy.units.Quantity`, optional
        Exposure time of the observation
    projection_code : `str`, optional
        The FITS standard projection code for the new header.

    Returns
    -------
    `~sunpy.util.MetaDict`
        The header information required for making a `sunpy.map.GenericMap`.

    Notes
    -----
    The observer coordinate is taken from the observer property of the ``reference_pixel``
    argument.

    Examples
    --------
    >>> import sunpy.map
    >>> from sunpy.coordinates import frames
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
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

    if hasattr(data, "shape"):
        shape = data.shape
    else:
        shape = data

    meta_wcs = _get_wcs_meta(coordinate, projection_code)

    meta_instrument = _get_instrument_meta(instrument, telescope, observatory, wavelength, exposure)
    meta_wcs.update(meta_instrument)

    if reference_pixel is None:
        reference_pixel = u.Quantity([(shape[1] - 1)/2.*u.pixel, (shape[0] - 1)/2.*u.pixel])
    if scale is None:
        scale = [1., 1.] * (u.arcsec/u.pixel)

    meta_wcs['crval1'], meta_wcs['crval2'] = (coordinate.spherical.lon.to_value(meta_wcs['cunit1']),
                                              coordinate.spherical.lat.to_value(meta_wcs['cunit2']))

    # Add 1 to go from input 0-based indexing to FITS 1-based indexing
    meta_wcs['crpix1'], meta_wcs['crpix2'] = (reference_pixel[0].to_value(u.pixel) + 1,
                                              reference_pixel[1].to_value(u.pixel) + 1)

    meta_wcs['cdelt1'], meta_wcs['cdelt2'] = (scale[0].to_value(meta_wcs['cunit1']/u.pixel),
                                              scale[1].to_value(meta_wcs['cunit2']/u.pixel))

    if rotation_angle is not None and rotation_matrix is not None:
        raise ValueError("Can not specify both rotation angle and rotation matrix.")

    if rotation_angle is not None:
        lam = meta_wcs['cdelt1'] / meta_wcs['cdelt2']
        p = np.deg2rad(rotation_angle)

        rotation_matrix = np.array([[np.cos(p), -1 * lam * np.sin(p)],
                                    [1/lam * np.sin(p), np.cos(p)]])

    if rotation_matrix is not None:
        (meta_wcs['PC1_1'], meta_wcs['PC1_2'],
         meta_wcs['PC2_1'], meta_wcs['PC2_2']) = (rotation_matrix[0, 0], rotation_matrix[0, 1],
                                                  rotation_matrix[1, 0], rotation_matrix[1, 1])

    meta_dict = MetaDict(meta_wcs)

    return meta_dict


def _get_wcs_meta(coordinate, projection_code):
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

    skycoord_wcs = astropy.wcs.utils.celestial_frame_to_wcs(coordinate, projection_code)

    cunit1, cunit2 = skycoord_wcs.wcs.cunit
    coord_meta = dict(skycoord_wcs.to_header())
    coord_meta['cunit1'], coord_meta['cunit2'] = cunit1.to_string("fits"), cunit2.to_string("fits")

    return coord_meta


@u.quantity_input
def get_observer_meta(observer, rsun: (u.Mm, None)):
    """
    Function to get observer meta from coordinate frame.

    Parameters
    ----------
    coordinate : ~`astropy.coordinates.BaseFrame`
        The coordinate of the observer, must be transformable to Heliographic
        Stonyhurst.
    rsun : `astropy.units.Quantity`
        The radius of the Sun.

    Returns
    -------
    `dict`
        Containing the WCS meta information
            * hgln_obs, hglt_obs
            * dsun_obs
            * rsun_obs
            * rsun_ref
    """
    observer = observer.transform_to(frames.HeliographicStonyhurst(obstime=observer.obstime))
    coord_meta = {}

    coord_meta['hgln_obs'] = observer.lon.to_value(u.deg)
    coord_meta['hglt_obs'] = observer.lat.to_value(u.deg)
    coord_meta['dsun_obs'] = observer.radius.to_value(u.m)
    if rsun is not None:
        coord_meta['rsun_ref'] = rsun.to_value(u.m)
        coord_meta['rsun_obs'] = sun._angular_radius(rsun, observer.radius).to_value(u.arcsec)

    return coord_meta


def _get_instrument_meta(instrument, telescope, observatory, wavelength, exposure):
    """
    Function to correctly name keywords from keyword arguments
    """
    coord = {}

    if instrument is not None:
        coord['instrume'] = str(instrument)
    if telescope is not None:
        coord['telescop'] = str(telescope)
    if observatory is not None:
        coord['obsrvtry'] = str(observatory)
    if wavelength is not None:
        coord['wavelnth'] = wavelength.to_value()
        coord['waveunit'] = wavelength.unit.to_string("fits")
    if exposure is not None:
        coord['exptime'] = exposure.to_value(u.s)

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
