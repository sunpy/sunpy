import numpy as np

import astropy.units as u
import astropy.wcs
from astropy.coordinates import SkyCoord

from sunpy import log
from sunpy.coordinates import frames, sun
from sunpy.util import MetaDict

__all__ = ['make_fitswcs_header', 'get_observer_meta', 'make_heliographic_header', 'make_hpr_header']


@u.quantity_input(equivalencies=u.spectral())
def make_fitswcs_header(data,
                        coordinate,
                        reference_pixel: u.pix = None,
                        scale: u.arcsec/u.pix = None,
                        rotation_angle: u.deg = None,
                        rotation_matrix=None,
                        instrument=None,
                        telescope=None,
                        observatory=None,
                        detector=None,
                        wavelength: u.angstrom = None,
                        exposure: u.s = None,
                        projection_code="TAN",
                        unit=None):
    """
    Function to create a FITS-WCS header from a coordinate object
    (`~astropy.coordinates.SkyCoord`) that is required to
    create a `~sunpy.map.GenericMap`.

    Parameters
    ----------
    data : `~numpy.ndarray`, `~astropy.units.Quantity`, or `tuple`
        Array data of Map for which a header is required, or the shape of the
        data array (in numpy order, i.e. ``(y_size, x_size)``).
    coordinate : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
        The coordinate of the reference pixel.
    reference_pixel : `~astropy.units.Quantity`, optional
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
        ``rotation_matrix``. Defaults to no rotation.
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
    detector : `str`, optional
        Name of the detector of the observation.
    wavelength : `~astropy.units.Quantity`, optional
        Wavelength of the observation as an astropy quantity, e.g. 171*u.angstrom.
        From this keyword, the meta keywords ``wavelnth`` and ``waveunit`` will be populated.
    exposure : `~astropy.units.Quantity`, optional
        Exposure time of the observation
    projection_code : `str`, optional
        The FITS standard projection code for the new header.
    unit : `~astropy.units.Unit`, optional
        Units of the array data of the Map. This will populate the the ``'bunit'`` meta keyword.
        If ``data`` is a `~astropy.units.Quantity`, the unit specified here will take precedence
        over the unit information attached to ``data``.

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
    coordinate = _validate_coordinate(coordinate)

    if hasattr(data, "shape"):
        shape = data.shape
    else:
        shape = data

    if hasattr(data, "unit"):
        if unit is None:
            unit = data.unit
        else:
            log.info("Overwriting data's current unit with specified unit.")

    meta_wcs = _get_wcs_meta(coordinate, projection_code)

    meta_wcs = _set_instrument_meta(meta_wcs, instrument, telescope, observatory, detector, wavelength, exposure, unit)
    meta_wcs = _set_transform_params(meta_wcs, coordinate, reference_pixel, scale, shape)
    meta_wcs = _set_rotation_params(meta_wcs, rotation_angle, rotation_matrix)

    if getattr(coordinate, 'observer', None) is not None:
        # Have to check for str, as doing == on a SkyCoord and str raises an error
        if isinstance(coordinate.observer, str) and coordinate.observer == 'self':
            dsun_obs = coordinate.radius
        else:
            dsun_obs = coordinate.observer.radius
        meta_wcs['rsun_obs'] = sun._angular_radius(coordinate.rsun, dsun_obs).to_value(u.arcsec)

    meta_dict = MetaDict(meta_wcs)
    return meta_dict


def _validate_coordinate(coordinate):
    if not isinstance(coordinate, SkyCoord | frames.BaseCoordinateFrame):
        raise ValueError("coordinate needs to be a coordinate frame or an SkyCoord instance.")

    if isinstance(coordinate, SkyCoord):
        coordinate = coordinate.frame

    if coordinate.obstime is None:
        raise ValueError("The coordinate needs an observation time, `obstime`.")

    if isinstance(coordinate, frames.Heliocentric):
        raise ValueError("This function does not currently support heliocentric coordinates.")
    return coordinate


def _set_rotation_params(meta_wcs, rotation_angle, rotation_matrix):
    if rotation_angle is not None and rotation_matrix is not None:
        raise ValueError("Can not specify both rotation angle and rotation matrix.")
    if rotation_angle is None and rotation_matrix is None:
        rotation_angle = 0 * u.deg

    if rotation_angle is not None:
        lam = meta_wcs['cdelt2'] / meta_wcs['cdelt1']
        p = np.deg2rad(rotation_angle)

        rotation_matrix = np.array([[np.cos(p), -1 * lam * np.sin(p)],
                                    [1/lam * np.sin(p), np.cos(p)]])

    if rotation_matrix is not None:
        (meta_wcs['PC1_1'], meta_wcs['PC1_2'],
         meta_wcs['PC2_1'], meta_wcs['PC2_2']) = (rotation_matrix[0, 0], rotation_matrix[0, 1],
                                                  rotation_matrix[1, 0], rotation_matrix[1, 1])
    return meta_wcs


def _set_transform_params(meta_wcs, coordinate, reference_pixel, scale, shape):
    meta_wcs['naxis'] = 2
    meta_wcs['naxis1'] = shape[1]
    meta_wcs['naxis2'] = shape[0]

    if reference_pixel is None:
        reference_pixel = u.Quantity([(shape[1] - 1)/2.*u.pixel,
                                      (shape[0] - 1)/2.*u.pixel])
    if scale is None:
        scale = [1., 1.] * (u.arcsec/u.pixel)

    meta_wcs['crval1'], meta_wcs['crval2'] = (coordinate.spherical.lon.to_value(meta_wcs['cunit1']),
                                              coordinate.spherical.lat.to_value(meta_wcs['cunit2']))

    # When the native latitude of the fiducial point is 0 degrees, which is typical for cylindrical
    # projections, the correct value of `lonpole` depends on whether the world latitude of the
    # fiducial point is greater than or less than its native latitude.
    if coordinate.spherical.lat.to_value(u.deg) < meta_wcs['theta0']:
        meta_wcs['LONPOLE'] = 180.
    else:
        meta_wcs['LONPOLE'] = 0.
    del meta_wcs['theta0']  # remove the native latitude of the fiducial point

    # Add 1 to go from input 0-based indexing to FITS 1-based indexing
    meta_wcs['crpix1'], meta_wcs['crpix2'] = (reference_pixel[0].to_value(u.pixel) + 1,
                                              reference_pixel[1].to_value(u.pixel) + 1)

    meta_wcs['cdelt1'] = scale[0].to_value(u.Unit(meta_wcs['cunit1']) / u.pixel)
    meta_wcs['cdelt2'] = scale[1].to_value(u.Unit(meta_wcs['cunit2']) / u.pixel)
    return meta_wcs


def _get_wcs_meta(coordinate, projection_code):
    """
    Function to get WCS meta from the SkyCoord using
    `astropy.wcs.utils.celestial_frame_to_wcs`

    Parameters
    ----------
    coordinate : `~astropy.coordinates.BaseCoordinateFrame`

    Returns
    -------
    `dict`
        Containing the WCS meta information
            * ctype1, ctype2
            * cunit1, cunit2
            * date_obs
            * observer auxiliary information, if set on `coordinate`
    """

    coord_meta = {}

    skycoord_wcs = astropy.wcs.utils.celestial_frame_to_wcs(coordinate, projection_code)

    cunit1, cunit2 = skycoord_wcs.wcs.cunit
    coord_meta = dict(skycoord_wcs.to_header())
    coord_meta['cunit1'], coord_meta['cunit2'] = cunit1.to_string("fits"), cunit2.to_string("fits")
    coord_meta['theta0'] = skycoord_wcs.wcs.theta0  # add the native latitude of the fiducial point

    return coord_meta


@u.quantity_input
def get_observer_meta(observer, rsun: (u.Mm, None) = None):
    """
    Function to get observer meta from coordinate frame.

    Parameters
    ----------
    observer : `~astropy.coordinates.BaseCoordinateFrame`
        The coordinate of the observer, must be transformable to Heliographic
        Stonyhurst.
    rsun : `astropy.units.Quantity`, optional
        The radius of the Sun. If ``None``, the RSUN_OBS and RSUN_REF keys are
        not set.

    Returns
    -------
    coord_meta : `dict`
        WCS metadata, with the keys ``['hgln_obs', 'hglt_obs', 'dsun_obs']``,
        and additionally if ``rsun`` is given ``['rsun_obs', 'rsun_ref']``.
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


def _set_instrument_meta(meta_wcs, instrument, telescope, observatory, detector, wavelength, exposure, unit):
    """
    Function to correctly name keywords from keyword arguments
    """
    if instrument is not None:
        meta_wcs['instrume'] = str(instrument)
    if telescope is not None:
        meta_wcs['telescop'] = str(telescope)
    if observatory is not None:
        meta_wcs['obsrvtry'] = str(observatory)
    if detector is not None:
        meta_wcs['detector'] = str(detector)
    if wavelength is not None:
        meta_wcs['wavelnth'] = wavelength.to_value()
        meta_wcs['waveunit'] = wavelength.unit.to_string("fits")
    if exposure is not None:
        meta_wcs['exptime'] = exposure.to_value(u.s)
    if unit is not None:
        _unit = u.Unit(unit)
        # NOTE: This conditional is because 'DN' is not part of the
        # FITS standard and thus cannot be serialized to a string
        # using the FITS format
        unit_format = 'generic' if u.DN in _unit.bases else 'fits'
        meta_wcs['bunit'] = _unit.to_string(format=unit_format)
    return meta_wcs


@u.quantity_input
def make_heliographic_header(date, observer_coordinate, shape, *, frame, projection_code="CAR",
                             map_center_longitude: u.Quantity[u.deg] = 0.0*u.deg):
    """
    Construct a FITS-WCS header for a full-Sun heliographic (Carrington or Stonyhurst) coordinate frame.

    The date-time and observer coordinate of the new coordinate frame
    are taken from the input map. The resulting WCS covers the full surface
    of the Sun, and has a reference coordinate at (0, 0) degrees Longitude/Latitude.

    Parameters
    ----------
    date :
        Date for the output header.
    observer_coordinate :
        Observer coordinate for the output header.
    shape : [int, int]
        Output map shape, number of pixels in (latitude, longitude).
    frame : {'carrington', 'stonyhurst'}
        Coordinate frame.
    projection_code : {'CAR', 'CEA'}
        Projection to use for the latitude coordinate.
    map_center_longitude : `~astropy.units.Quantity`
        Heliographic longitude of the map center

    Returns
    -------
    `~sunpy.util.MetaDict`

    See Also
    --------
    sunpy.map.header_helper.make_fitswcs_header : A more generic header helper that can be used if more customisation is required.

    Examples
    --------
    >>> from sunpy.map.header_helper import make_heliographic_header
    >>> from sunpy.coordinates import get_earth
    >>>
    >>> date = '2020-01-01 12:00:00'
    >>> observer = get_earth(date)
    >>> header = make_heliographic_header(date, observer, [90, 180], frame='carrington')
    >>> header
    MetaDict([('wcsaxes': '2')
    ('crpix1': '90.5')
    ('crpix2': '45.5')
    ('cdelt1': '2.0')
    ('cdelt2': '2.0')
    ('cunit1': 'deg')
    ('cunit2': 'deg')
    ('ctype1': 'CRLN-CAR')
    ('ctype2': 'CRLT-CAR')
    ('crval1': '0.0')
    ('crval2': '0.0')
    ('lonpole': '0.0')
    ('latpole': '90.0')
    ('mjdref': '0.0')
    ('date-obs': '2020-01-01T12:00:00.000')
    ('rsun_ref': '695700000.0')
    ('dsun_obs': '147096975776.97')
    ('hgln_obs': '0.0')
    ('hglt_obs': '-3.0011725838606')
    ('naxis': '2')
    ('naxis1': '180')
    ('naxis2': '90')
    ('pc1_1': '1.0')
    ('pc1_2': '-0.0')
    ('pc2_1': '0.0')
    ('pc2_2': '1.0')
    ('rsun_obs': '975.53984320334...

    .. minigallery:: sunpy.map.make_heliographic_header
    """
    valid_codes = {"CAR", "CEA"}
    if projection_code not in valid_codes:
        raise ValueError(f"projection_code must be one of {valid_codes}")

    valid_frames = {'carrington', 'stonyhurst'}
    if frame not in valid_frames:
        raise ValueError(f"frame must be one of {valid_frames}")

    frame_out = SkyCoord(
        map_center_longitude,
        0 * u.deg,
        frame=f"heliographic_{frame}",
        obstime=date,
        observer=observer_coordinate,
        rsun=getattr(observer_coordinate, "rsun", None),
    )

    if projection_code == "CAR":
        scale = [360 / int(shape[1]), 180 / int(shape[0])] * u.deg / u.pix
    elif projection_code == "CEA":
        # Using the cylindrical equal-area (CEA) projection,
        # scale needs to be to 180/pi times the sin(latitude) spacing
        # See Section 5.5, Thompson 2006
        scale = [
            360 / int(shape[1]),
            (180 / np.pi) / (int(shape[0]) / 2)
        ] * u.deg / u.pix

    header = make_fitswcs_header(shape, frame_out, scale=scale, projection_code=projection_code)
    return header


@u.quantity_input
def make_hpr_header(observer_coordinate, shape, theta_binsize: u.Quantity[u.arcsec], *,
                    theta_min: u.Quantity[u.arcsec] = 0.0*u.arcsec,
                    psi_center: u.Quantity[u.deg] = 180.0*u.deg):
    """
    Construct a FITS-WCS header for a helioprojective radial coordinate frame.

    The header uses the plate carrée projection ("CAR" in FITS-WCS) with the
    position angle (``psi``) on the horizontal axis and impact angle (``theta``) on
    the vertical axis. To accommodate the FITS-WCS machinery, the vertical axis is
    actually the declination (``delta``), which is the impact angle minus 90 degrees.

    Parameters
    ----------
    observer_coordinate :
        Observer coordinate, with corresponding observation time.
    shape : [int, int]
        Output map shape, number of pixels in (theta, psi).
    theta_binsize : `~astropy.units.Quantity`
        The size of each pixel in the impact-angle direction.
    theta_min : `~astropy.units.Quantity`
        The minimum impact angle. Defaults to 0 arcsec.
    psi_center : `~astropy.units.Quantity`
        The center of the map in position angle. Defaults to 180 degrees.

    Returns
    -------
    `~sunpy.util.MetaDict`

    See Also
    --------
    sunpy.map.header_helper.make_fitswcs_header : A more generic header helper that can be used if more customisation is required.
    sunpy.coordinates.HelioprojectiveRadial : The coordinate-frame class

    Examples
    --------
    >>> import astropy.units as u
    >>> from sunpy.coordinates import get_earth

    >>> observer = get_earth('2024-09-16 01:02:03')
    >>> make_hpr_header(observer, (180, 360), theta_binsize=10*u.arcsec)
    MetaDict([('wcsaxes': '2')
    ('crpix1': '180.5')
    ('crpix2': '32400.5')
    ('cdelt1': '1.0')
    ('cdelt2': '10.0')
    ('cunit1': 'deg')
    ('cunit2': 'arcsec')
    ('ctype1': 'HRLN-CAR')
    ('ctype2': 'HRLT-CAR')
    ('crval1': '180.0')
    ('crval2': '0.0')
    ('lonpole': '0.0')
    ('latpole': '90.0')
    ('mjdref': '0.0')
    ('date-obs': '2024-09-16T01:02:03.000')
    ('rsun_ref': '695700000.0')
    ('dsun_obs': '150399816551.15')
    ('hgln_obs': '0.0')
    ('hglt_obs': '7.1884873184853')
    ('naxis': '2')
    ('naxis1': '360')
    ('naxis2': '180')
    ('pc1_1': '1.0')
    ('pc1_2': '-0.0')
    ('pc2_1': '0.0')
    ('pc2_2': '1.0')
    ('rsun_obs': '954.1164393469891')])

    .. minigallery:: sunpy.map.make_hpr_header
    """
    # The reference coord must be on the equator (delta = 0 deg)
    reference_coord = SkyCoord(
        psi_center,
        0*u.deg,
        frame="helioprojectiveradial",
        obstime=observer_coordinate.obstime,
        observer=observer_coordinate,
        rsun=getattr(observer_coordinate, "rsun", None),
    )

    # Calculate the correct location for this reference coord, which will likely be far off the data array
    reference_pixel = [
        shape[1] / 2. - 0.5,
        ((90*u.deg - theta_min) / theta_binsize).to_value(u.one) - 0.5
    ] * u.pix

    scale = u.Quantity([360*u.deg / int(shape[1]), theta_binsize]) / u.pix

    header = make_fitswcs_header(shape, reference_coord, reference_pixel=reference_pixel,
                                 scale=scale, projection_code="CAR")
    return header
