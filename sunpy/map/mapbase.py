"""
Map is a generic Map class from which all other Map classes inherit from.
"""
import copy
import warnings
from collections import namedtuple
import textwrap

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import astropy.wcs
import astropy.units as u
from astropy.visualization import AsymmetricPercentileInterval
from astropy.visualization.wcsaxes import WCSAxes
from astropy.coordinates import SkyCoord, UnitSphericalRepresentation

import sunpy.io as io
# The next two are not used but are called to register functions with external modules
import sunpy.coordinates
import sunpy.visualization.colormaps
from sunpy import config
from sunpy.visualization import wcsaxes_compat, axis_labels_from_ctype, peek_show
from sunpy.sun import constants
from sunpy.coordinates import sun
from sunpy.time import parse_time, is_time
from sunpy.image.resample import reshape_image_to_4d_superpixel
from sunpy.image.resample import resample as sunpy_image_resample
from sunpy.coordinates import get_earth
from sunpy.util import expand_list
from sunpy.util.exceptions import SunpyUserWarning

from astropy.nddata import NDData

TIME_FORMAT = config.get("general", "time_format")
PixelPair = namedtuple('PixelPair', 'x y')
SpatialPair = namedtuple('SpatialPair', 'axis1 axis2')

__all__ = ['GenericMap']


class MapMetaValidationError(AttributeError):
    pass


class GenericMap(NDData):
    """
    A Generic spatially-aware 2D data array

    Parameters
    ----------
    data : `numpy.ndarray`, list
        A 2d list or ndarray containing the map data.
    header : dict
        A dictionary of the original image header tags.
    plot_settings : dict, optional
        Plot settings.

    Other Parameters
    ----------------
    **kwargs :
        Additional keyword arguments are passed to `~astropy.nddata.NDData`
        init.

    Examples
    --------
    >>> import sunpy.map
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> aia   # doctest: +REMOTE_DATA
    SunPy Map
    ---------
    Observatory:		 SDO
    Instrument:		 AIA 3
    Detector:		 AIA
    Measurement:		 171.0 Angstrom
    Wavelength:		 171.0 Angstrom
    Observation Date:	 2011-06-07 06:33:02
    Exposure Time:		 0.234256 s
    Dimension:		 [1024. 1024.] pix
    Coordinate System:	 helioprojective
    Scale:			 [2.402792 2.402792] arcsec / pix
    Reference Pixel:	 [512.5 512.5] pix
    Reference Coord:	 [3.22309951 1.38578135] arcsec
    array([[ -95.92475  ,    7.076416 ,   -1.9656711, ..., -127.96519  ,
            -127.96519  , -127.96519  ],
           [ -96.97533  ,   -5.1167884,    0.       , ...,  -98.924576 ,
            -104.04137  , -127.919716 ],
           [ -93.99607  ,    1.0189276,   -4.0757103, ...,   -5.094638 ,
             -37.95505  , -127.87541  ],
           ...,
           [-128.01454  , -128.01454  , -128.01454  , ..., -128.01454  ,
            -128.01454  , -128.01454  ],
           [-127.899666 , -127.899666 , -127.899666 , ..., -127.899666 ,
            -127.899666 , -127.899666 ],
           [-128.03072  , -128.03072  , -128.03072  , ..., -128.03072  ,
            -128.03072  , -128.03072  ]], dtype=float32)

    >>> aia.spatial_units   # doctest: +REMOTE_DATA
    SpatialPair(axis1=Unit("arcsec"), axis2=Unit("arcsec"))
    >>> aia.peek()   # doctest: +SKIP

    Notes
    -----

    A number of the properties of this class are returned as two-value named
    tuples that can either be indexed by position ([0] or [1]) or be accessed
    by the names (.x and .y) or (.axis1 and .axis2). Things that refer to pixel
    axes use the ``.x``, ``.y`` convention, where x and y refer to the FITS
    axes (x for columns y for rows). Spatial axes use ``.axis1`` and ``.axis2``
    which correspond to the first and second axes in the header. ``axis1``
    corresponds to the coordinate axis for ``x`` and ``axis2`` corresponds to
    ``y``.

    This class makes some assumptions about the WCS information contained in
    the meta data. The first and most extensive assumption is that it is
    FITS-like WCS information as defined in the FITS WCS papers.

    Within this scope it also makes some other assumptions.

    * In the case of APIS convention headers where the CROTAi/j arguments are
      provided it assumes that these can be converted to the standard PCi_j
      notation using equations 32 in Thompson (2006).

    * If a CDi_j matrix is provided it is assumed that it can be converted to a
      PCi_j matrix and CDELT keywords as described in
      `Greisen & Calabretta (2002) <https://doi.org/10.1051/0004-6361:20021327>`_

    * The 'standard' FITS keywords that are used by this class are the PCi_j
      matrix and CDELT, along with the other keywords specified in the WCS
      papers. All subclasses of this class must convert their header
      information to this formalism. The CROTA to PCi_j conversion is done in
      this class.

    .. warning::
        This class currently assumes that a header with the CDi_j matrix
        information also includes the CDELT keywords, without these keywords
        this class will not process the WCS.
        Also the rotation_matrix does not work if the CDELT1 and CDELT2
        keywords are exactly equal.
        Also, if a file with more than two dimensions is feed into the class,
        only the first two dimensions (NAXIS1, NAXIS2) will be loaded and the
        rest will be discarded.
    """

    _registry = dict()

    def __init_subclass__(cls, **kwargs):
        """
        An __init_subclass__ hook initializes all of the subclasses of a given class.
        So for each subclass, it will call this block of code on import.
        This replicates some metaclass magic without the need to be aware of metaclasses.
        Here we use this to register each subclass in a dict that has the
        `is_datasource_for` attribute.
        This is then passed into the Map Factory so we can register them.
        """
        super().__init_subclass__(**kwargs)
        if hasattr(cls, 'is_datasource_for'):
            cls._registry[cls] = cls.is_datasource_for

    def __init__(self, data, header, plot_settings=None, **kwargs):
        # If the data has more than two dimensions, the first dimensions
        # (NAXIS1, NAXIS2) are used and the rest are discarded.
        ndim = data.ndim
        if ndim > 2:
            # We create a slice that removes all but the 'last' two
            # dimensions. (Note dimensions in ndarray are in reverse order)

            new_2d_slice = [0]*(ndim-2)
            new_2d_slice.extend([slice(None), slice(None)])
            data = data[tuple(new_2d_slice)]
            # Warn the user that the data has been truncated
            warnings.warn("This file contains more than 2 dimensions. "
                          "Data will be truncated to the first two dimensions.", SunpyUserWarning)

        super().__init__(data, meta=header, **kwargs)

        # Correct possibly missing meta keywords
        self._fix_date()
        self._fix_naxis()

        # Setup some attributes
        self._nickname = None
        # These are palceholders for default attributes, which are only set
        # once if their data isn't present in the map metadata.
        self._default_time = None
        self._default_dsun = None
        self._default_carrington_longitude = None
        self._default_heliographic_latitude = None
        self._default_heliographic_longitude = None

        # Validate header
        # TODO: This should be a function of the header, not of the map
        self._validate_meta()
        self._shift = SpatialPair(0 * u.arcsec, 0 * u.arcsec)

        if self.dtype == np.uint8:
            norm = None
        else:
            # Put import here to reduce sunpy.map import time
            from matplotlib import colors
            norm = colors.Normalize()

        # Visualization attributes
        self.plot_settings = {'cmap': 'gray',
                              'norm': norm,
                              'interpolation': 'nearest',
                              'origin': 'lower'
                              }
        if plot_settings:
            self.plot_settings.update(plot_settings)

    def __getitem__(self, key):
        """ This should allow indexing by physical coordinate """
        raise NotImplementedError(
            "The ability to index Map by physical"
            " coordinate is not yet implemented.")

    def __repr__(self):
        return textwrap.dedent("""\
                   SunPy Map
                   ---------
                   Observatory:\t\t {obs}
                   Instrument:\t\t {inst}
                   Detector:\t\t {det}
                   Measurement:\t\t {meas}
                   Wavelength:\t\t {wave}
                   Observation Date:\t {date}
                   Exposure Time:\t\t {dt:f}
                   Dimension:\t\t {dim}
                   Coordinate System:\t {coord}
                   Scale:\t\t\t {scale}
                   Reference Pixel:\t {refpix}
                   Reference Coord:\t {refcoord}
                   """).format(obs=self.observatory, inst=self.instrument, det=self.detector,
                               meas=self.measurement, wave=self.wavelength,
                               date=self.date.strftime(TIME_FORMAT),
                               dt=self.exposure_time,
                               dim=u.Quantity(self.dimensions),
                               scale=u.Quantity(self.scale),
                               coord=self._coordinate_frame_name,
                               refpix=u.Quantity(self.reference_pixel),
                               refcoord=u.Quantity((self.reference_coordinate.data.lon,
                                                    self.reference_coordinate.data.lat)),
                               tmf=TIME_FORMAT) + self.data.__repr__()

    @classmethod
    def _new_instance(cls, data, meta, plot_settings=None, **kwargs):
        """
        Instantiate a new instance of this class using given data.
        This is a shortcut for ``type(self)(data, meta, plot_settings)``.
        """
        return cls(data, meta, plot_settings=plot_settings, **kwargs)

    def _get_lon_lat(self, frame):
        """
        Given a coordinate frame, extract the lon and lat by casting to
        SphericalRepresentation first.
        """
        r = frame.represent_as(UnitSphericalRepresentation)
        return r.lon.to(self.spatial_units[0]), r.lat.to(self.spatial_units[1])

    @property
    def wcs(self):
        """
        The `~astropy.wcs.WCS` property of the map.
        """
        # Construct the WCS based on the FITS header, but don't "do_set" which
        # analyses the FITS header for correctness.
        with warnings.catch_warnings():
            # Ignore warnings we may raise when constructing the fits header about dropped keys.
            warnings.simplefilter("ignore", SunpyUserWarning)
            try:
                w2 = astropy.wcs.WCS(header=self.fits_header, _do_set=False)
            except Exception:
                warnings.warn("Unable to treat `.meta` as a FITS header, assuming a simple WCS.")
                w2 = astropy.wcs.WCS(naxis=2)

        # If the FITS header is > 2D pick the first 2 and move on.
        # This will require the FITS header to be valid.
        if w2.naxis > 2:
            # We have to change this or else the WCS doesn't parse properly, even
            # though we don't care about the third dimension. This applies to both
            # EIT and IRIS data, it is here to reduce the chances of silly errors.
            if self.meta.get('cdelt3', None) == 0:
                w2.wcs.cdelt[2] = 1e-10

            w2 = w2.sub([1, 2])

        w2.wcs.crpix = u.Quantity(self.reference_pixel)
        # Make these a quantity array to prevent the numpy setting element of
        # array with sequence error.
        w2.wcs.cdelt = u.Quantity(self.scale)
        w2.wcs.crval = u.Quantity([self._reference_longitude, self._reference_latitude])
        w2.wcs.ctype = self.coordinate_system
        w2.wcs.pc = self.rotation_matrix
        w2.wcs.cunit = self.spatial_units
        w2.wcs.dateobs = self.date.iso
        w2.rsun = self.rsun_meters

        # Astropy WCS does not understand the SOHO default of "solar-x" and
        # "solar-y" ctypes.  This overrides the default assignment and
        # changes it to a ctype that is understood.  See Thompson, 2006, A.&A.,
        # 449, 791.
        if w2.wcs.ctype[0].lower() in ("solar-x", "solar_x"):
            w2.wcs.ctype[0] = 'HPLN-TAN'

        if w2.wcs.ctype[1].lower() in ("solar-y", "solar_y"):
            w2.wcs.ctype[1] = 'HPLT-TAN'

        # GenericMap.coordinate_frame is implemented using this method, so we
        # need to do this only based on .meta.
        ctypes = {c[:4] for c in w2.wcs.ctype}
        # Check that the ctypes contains one of these two pairs of axes.
        if {'HPLN', 'HPLT'} <= ctypes or {'SOLX', 'SOLY'} <= ctypes:
            w2.heliographic_observer = self.observer_coordinate

        # Validate the WCS here.
        w2.wcs.set()
        return w2

    @property
    def coordinate_frame(self):
        """
        An `astropy.coordinates.BaseFrame` instance created from the coordinate
        information for this Map, or None if the frame cannot be determined.
        """
        try:
            return astropy.wcs.utils.wcs_to_celestial_frame(self.wcs)
        except ValueError as e:
            warnings.warn(f'Could not determine coordinate frame from map metadata',
                          SunpyUserWarning)
            return None

    @property
    def _coordinate_frame_name(self):
        if self.coordinate_frame is None:
            return 'Unknown'
        return self.coordinate_frame.name

    def _as_mpl_axes(self):
        """
        Compatibility hook for Matplotlib and WCSAxes.
        This functionality requires the WCSAxes package to work. The reason
        we include this here is that it allows users to use WCSAxes without
        having to explicitly import WCSAxes
        With this method, one can do::

            import matplotlib.pyplot as plt
            import sunpy.map
            amap = sunpy.map.Map('filename.fits')
            fig = plt.figure()
            ax = plt.subplot(projection=amap)
            ...

        and this will generate a plot with the correct WCS coordinates on the
        axes. See https://wcsaxes.readthedocs.io for more information.
        """
        # This code is reused from Astropy

        return WCSAxes, {'wcs': self.wcs}

    # Some numpy extraction
    @property
    def dimensions(self):
        """
        The dimensions of the array (x axis first, y axis second).
        """
        return PixelPair(*u.Quantity(np.flipud(self.data.shape), 'pixel'))

    @property
    def dtype(self):
        """
        The `numpy.dtype` of the array of the map.
        """
        return self.data.dtype

    @property
    def size(self):
        """
        The number of pixels in the array of the map.
        """
        return u.Quantity(self.data.size, 'pixel')

    @property
    def ndim(self):
        """
        The value of `numpy.ndarray.ndim` of the data array of the map.
        """
        return self.data.ndim

    def std(self, *args, **kwargs):
        """
        Calculate the standard deviation of the data array.
        """
        return self.data.std(*args, **kwargs)

    def mean(self, *args, **kwargs):
        """
        Calculate the mean of the data array.
        """
        return self.data.mean(*args, **kwargs)

    def min(self, *args, **kwargs):
        """
        Calculate the minimum value of the data array.
        """
        return self.data.min(*args, **kwargs)

    def max(self, *args, **kwargs):
        """
        Calculate the maximum value of the data array.
        """
        return self.data.max(*args, **kwargs)

# #### Keyword attribute and other attribute definitions #### #

    def _base_name(self):
        """Abstract the shared bit between name and latex_name"""
        return "{nickname} {{measurement}} {date}".format(
            nickname=self.nickname,
            date=parse_time(self.date).strftime(TIME_FORMAT)
        )

    @property
    def name(self):
        """Human-readable description of the Map."""
        return self._base_name().format(measurement=self.measurement)

    @property
    def latex_name(self):
        """LaTeX formatted description of the Map."""
        if isinstance(self.measurement, u.Quantity):
            return self._base_name().format(measurement=self.measurement._repr_latex_())
        else:
            return self.name

    @property
    def nickname(self):
        """An abbreviated human-readable description of the map-type; part of
        the Helioviewer data model."""
        return self._nickname if self._nickname else self.detector

    @nickname.setter
    def nickname(self, n):
        self._nickname = n

    @property
    def date(self):
        """Image observation time."""
        time = self.meta.get('date-obs', None)
        if time is None:
            if self._default_time is None:
                warnings.warn("Missing metadata for observation time: setting observation time to current time.",
                              SunpyUserWarning)
                self._default_time = parse_time('now')
            time = self._default_time
        return parse_time(time)

    @property
    def detector(self):
        """Detector name."""
        return self.meta.get('detector', "")

    @property
    def exposure_time(self):
        """Exposure time of the image in seconds."""
        return self.meta.get('exptime', 0.0) * u.s

    @property
    def instrument(self):
        """Instrument name."""
        return self.meta.get('instrume', "").replace("_", " ")

    @property
    def measurement(self):
        """Measurement name, defaults to the wavelength of image."""
        return u.Quantity(self.meta.get('wavelnth', 0),
                          self.waveunit)

    @property
    def waveunit(self):
        """The `~astropy.units.Unit` of the wavelength of this observation."""
        unit = self.meta.get("waveunit")
        if unit is None:
            return u.one
        return u.Unit(unit)

    @property
    def wavelength(self):
        """Wavelength of the observation."""
        return u.Quantity(self.meta.get('wavelnth', 0),
                          self.waveunit)

    @property
    def observatory(self):
        """Observatory or Telescope name."""
        return self.meta.get('obsrvtry',
                             self.meta.get('telescop', "")).replace("_", " ")

    @property
    def processing_level(self):
        """
        Returns the FITS processing level if present.
        """
        return self.meta.get('lvl_num', None)

    @property
    def bottom_left_coord(self):
        """
        The physical coordinate at the center of the bottom left ([0, 0]) pixel.
        """
        return self.pixel_to_world(0*u.pix, 0*u.pix)

    @property
    def top_right_coord(self):
        """
        The physical coordinate at the center of the the top left ([-1, -1]) pixel.
        """
        return self.pixel_to_world(*self.dimensions)

    @property
    def center(self):
        """
        Return a coordinate object for the center pixel of the array.
        """
        center = u.Quantity(self.dimensions) / 2.
        return self.pixel_to_world(*center)

    @property
    def shifted_value(self):
        """The total shift applied to the reference coordinate by past applications of
        `~sunpy.map.GenericMap.shift`."""
        return self._shift

    @u.quantity_input
    def shift(self, axis1: u.deg, axis2: u.deg):
        """
        Returns a map shifted by a specified amount to, for example, correct
        for a bad map location. These values are applied directly to the
        `~sunpy.map.GenericMap.reference_coordinate`. To check how much shift
        has already been applied see `~sunpy.map.GenericMap.shifted_value`

        Parameters
        ----------
        axis1 : `~astropy.units.Quantity`
            The shift to apply to the Longitude (solar-x) coordinate.

        axis2 : `~astropy.units.Quantity`
            The shift to apply to the Latitude (solar-y) coordinate

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new shifted Map.
        """
        new_meta = self.meta.copy()

        # Update crvals
        new_meta['crval1'] = ((self.meta['crval1'] *
                               self.spatial_units[0] + axis1).to(self.spatial_units[0])).value
        new_meta['crval2'] = ((self.meta['crval2'] *
                               self.spatial_units[1] + axis2).to(self.spatial_units[1])).value

        # Create new map with the modification
        new_map = self._new_instance(self.data, new_meta, self.plot_settings)

        new_map._shift = SpatialPair(self.shifted_value[0] + axis1,
                                     self.shifted_value[1] + axis2)

        return new_map

    @property
    def rsun_meters(self):
        """Radius of the sun in meters."""
        return u.Quantity(self.meta.get('rsun_ref', constants.radius), 'meter')

    @property
    def rsun_obs(self):
        """Radius of the Sun."""
        rsun_arcseconds = self.meta.get('rsun_obs',
                                        self.meta.get('solar_r',
                                                      self.meta.get('radius',
                                                                    None)))

        if rsun_arcseconds is None:
            warnings.warn("Missing metadata for solar radius: assuming photospheric limb as seen from Earth.",
                          SunpyUserWarning)
            rsun_arcseconds = sun.angular_radius(self.date).to('arcsec').value

        return u.Quantity(rsun_arcseconds, 'arcsec')

    @property
    def coordinate_system(self):
        """Coordinate system used for x and y axes (ctype1/2)."""
        return SpatialPair(self.meta.get('ctype1', 'HPLN-   '),
                           self.meta.get('ctype2', 'HPLT-   '))

    @property
    def _supported_observer_coordinates(self):
        """
        A list of supported coordinate systems.

        This is a list so it can easily maintain a strict order. The list of
        two element tuples, the first item in the tuple is the keys that need
        to be in the header to use this coordinate system and the second is the
        kwargs to SkyCoord.
        """
        return [(('hgln_obs', 'hglt_obs', 'dsun_obs'), {'lon': self.meta.get('hgln_obs'),
                                                        'lat': self.meta.get('hglt_obs'),
                                                        'radius': self.meta.get('dsun_obs'),
                                                        'unit': (u.deg, u.deg, u.m),
                                                        'frame': "heliographic_stonyhurst"}),
                (('crln_obs', 'crlt_obs', 'dsun_obs'), {'lon': self.meta.get('crln_obs'),
                                                        'lat': self.meta.get('crlt_obs'),
                                                        'radius': self.meta.get('dsun_obs'),
                                                        'unit': (u.deg, u.deg, u.m),
                                                        'frame': "heliographic_carrington"}),]

    def _remove_existing_observer_location(self):
        """
        Remove all keys that this map might use for observer location.
        """
        all_keys = expand_list([e[0] for e in self._supported_observer_coordinates])
        for key in all_keys:
            self.meta.pop(key)

    @property
    def observer_coordinate(self):
        """
        The Heliographic Stonyhurst Coordinate of the observer.
        """
        missing_meta = {}
        for keys, kwargs in self._supported_observer_coordinates:
            meta_list = [k in self.meta for k in keys]
            if all(meta_list):
                return SkyCoord(obstime=self.date, **kwargs).heliographic_stonyhurst
            elif any(meta_list) and not set(keys).isdisjoint(self.meta.keys()):
                if not isinstance(kwargs['frame'], str):
                    kwargs['frame'] = kwargs['frame'].name
                missing_meta[kwargs['frame']] = set(keys).difference(self.meta.keys())

        warning_message = "".join([f"For frame '{frame}' the following metadata is missing: {','.join(keys)}\n" for frame, keys in missing_meta.items()])
        warning_message = "Missing metadata for observer: assuming Earth-based observer.\n" + warning_message
        warnings.warn(warning_message, SunpyUserWarning)

        return get_earth(self.date)

    @property
    def heliographic_latitude(self):
        """Observer heliographic latitude."""
        return self.observer_coordinate.lat

    @property
    def heliographic_longitude(self):
        """Observer heliographic longitude."""
        return self.observer_coordinate.lon

    @property
    def carrington_latitude(self):
        """Observer Carrington latitude."""
        return self.observer_coordinate.heliographic_carrington.lat

    @property
    def carrington_longitude(self):
        """Observer Carrington longitude."""
        return self.observer_coordinate.heliographic_carrington.lon

    @property
    def dsun(self):
        """Observer distance from the center of the Sun."""
        return self.observer_coordinate.radius.to('m')

    @property
    def _reference_longitude(self):
        """
        FITS-WCS compatible longitude. Used in self.wcs and
        self.reference_coordinate.
        """
        return self.meta.get('crval1', 0.) * self.spatial_units[0]

    @property
    def _reference_latitude(self):
        return self.meta.get('crval2', 0.) * self.spatial_units[1]

    @property
    def reference_coordinate(self):
        """Reference point WCS axes in data units (i.e. crval1, crval2). This value
        includes a shift if one is set."""
        return SkyCoord(self._reference_longitude,
                        self._reference_latitude,
                        frame=self.coordinate_frame)

    @property
    def reference_pixel(self):
        """Reference point axes in pixels (i.e. crpix1, crpix2)."""
        return PixelPair(self.meta.get('crpix1',
                                       (self.meta.get('naxis1') + 1) / 2.) * u.pixel,
                         self.meta.get('crpix2',
                                       (self.meta.get('naxis2') + 1) / 2.) * u.pixel)

    @property
    def scale(self):
        """
        Image scale along the x and y axes in units/pixel
        (i.e. cdelt1, cdelt2).
        """
        # TODO: Fix this if only CDi_j matrix is provided
        return SpatialPair(self.meta.get('cdelt1', 1.) * self.spatial_units[0] / u.pixel,
                           self.meta.get('cdelt2', 1.) * self.spatial_units[1] / u.pixel)

    @property
    def spatial_units(self):
        """
        Image coordinate units along the x and y axes (i.e. cunit1, cunit2).
        """
        return SpatialPair(u.Unit(self.meta.get('cunit1')),
                           u.Unit(self.meta.get('cunit2')))

    @property
    def rotation_matrix(self):
        """
        Matrix describing the rotation required to align solar North with
        the top of the image.
        """
        if 'PC1_1' in self.meta:
            return np.array([[self.meta['PC1_1'], self.meta['PC1_2']],
                             [self.meta['PC2_1'], self.meta['PC2_2']]])

        elif 'CD1_1' in self.meta:
            cd = np.array([[self.meta['CD1_1'], self.meta['CD1_2']],
                           [self.meta['CD2_1'], self.meta['CD2_2']]])

            cdelt = u.Quantity(self.scale).value

            return cd / cdelt
        else:
            return self._rotation_matrix_from_crota()

    def _rotation_matrix_from_crota(self):
        """
        This method converts the deprecated CROTA FITS kwargs to the new
        PC rotation matrix.

        This method can be overriden if an instruments header does not use this
        conversion.
        """
        lam = self.scale[0] / self.scale[1]
        p = np.deg2rad(self.meta.get('CROTA2', 0))

        return np.array([[np.cos(p), -1 * lam * np.sin(p)],
                         [1/lam * np.sin(p), np.cos(p)]])

    @property
    def fits_header(self):
        """
        A `~astropy.io.fits.Header` representation of the ``meta`` attribute.
        """
        return sunpy.io.fits.header_to_fits(self.meta)

# #### Miscellaneous #### #

    def _fix_date(self):
        # Check commonly used but non-standard FITS keyword for observation
        # time and correct the keyword if we can. Keep updating old one for
        # backwards compatibility.
        if is_time(self.meta.get('date_obs', None)):
            self.meta['date-obs'] = self.meta['date_obs']

    def _fix_naxis(self):
        # If naxis is not specified, get it from the array shape
        if 'naxis1' not in self.meta:
            self.meta['naxis1'] = self.data.shape[1]
        if 'naxis2' not in self.meta:
            self.meta['naxis2'] = self.data.shape[0]
        if 'naxis' not in self.meta:
            self.meta['naxis'] = self.ndim

    def _fix_bitpix(self):
        # Bit-depth
        #
        #   8    Character or unsigned binary integer
        #  16    16-bit twos-complement binary integer
        #  32    32-bit twos-complement binary integer
        # -32    IEEE single precision floating point
        # -64    IEEE double precision floating point
        #
        if 'bitpix' not in self.meta:
            float_fac = -1 if self.dtype.kind == "f" else 1
            self.meta['bitpix'] = float_fac * 8 * self.dtype.itemsize

    def _get_cmap_name(self):
        """Build the default color map name."""
        cmap_string = (self.observatory + self.detector +
                       str(int(self.wavelength.to('angstrom').value)))
        return cmap_string.lower()

    def _validate_meta(self):
        """
        Validates the meta-information associated with a Map.

        This method includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.map.sources package.

        Allows for default unit assignment for:
            CUNIT1, CUNIT2, WAVEUNIT

        """
        msg = ('Image coordinate units for axis {} not present in metadata.')
        err_message = []
        for i in [1, 2]:
            if self.meta.get(f'cunit{i}') is None:
                err_message.append(msg.format(i, i))

        if err_message:
            err_message.append(
                'See https://docs.sunpy.org/en/stable/code_ref/map.html#fixing-map-metadata` for '
                'instructions on how to add missing metadata.')
            raise MapMetaValidationError('\n'.join(err_message))

        for meta_property in ('waveunit', ):
            if (self.meta.get(meta_property) and
                u.Unit(self.meta.get(meta_property),
                       parse_strict='silent').physical_type == 'unknown'):
                warnings.warn(f"Unknown value for {meta_property.upper()}.", SunpyUserWarning)

        if (self.coordinate_system[0].startswith(('SOLX', 'SOLY')) or
            self.coordinate_system[1].startswith(('SOLX', 'SOLY'))):
            warnings.warn("SunPy Map does not support three dimensional data "
                          "and therefore cannot represent heliocentric coordinates. Proceed at your own risk.",
                          SunpyUserWarning)

# #### Data conversion routines #### #
    def world_to_pixel(self, coordinate, origin=0):
        """
        Convert a world (data) coordinate to a pixel coordinate by using
        `~astropy.wcs.WCS.wcs_world2pix`.

        Parameters
        ----------
        coordinate : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseFrame`
            The coordinate object to convert to pixel coordinates.

        origin : int
            Origin of the top-left corner. i.e. count from 0 or 1.
            Normally, origin should be 0 when passing numpy indices, or 1 if
            passing values from FITS header or map attributes.
            See `~astropy.wcs.WCS.wcs_world2pix` for more information.

        Returns
        -------
        x : `~astropy.units.Quantity`
            Pixel coordinate on the CTYPE1 axis.

        y : `~astropy.units.Quantity`
            Pixel coordinate on the CTYPE2 axis.
        """
        if not isinstance(coordinate, (SkyCoord,
                                       astropy.coordinates.BaseCoordinateFrame)):
            raise ValueError("world_to_pixel takes a Astropy coordinate frame or SkyCoord instance.")

        native_frame = coordinate.transform_to(self.coordinate_frame)
        lon, lat = u.Quantity(self._get_lon_lat(native_frame)).to(u.deg)
        x, y = self.wcs.wcs_world2pix(lon, lat, origin)

        return PixelPair(x * u.pixel, y * u.pixel)

    @u.quantity_input
    def pixel_to_world(self, x: u.pixel, y: u.pixel, origin=0):
        """
        Convert a pixel coordinate to a data (world) coordinate by using
        `~astropy.wcs.WCS.wcs_pix2world`.

        Parameters
        ----------

        x : `~astropy.units.Quantity`
            Pixel coordinate of the CTYPE1 axis. (Normally solar-x).

        y : `~astropy.units.Quantity`
            Pixel coordinate of the CTYPE2 axis. (Normally solar-y).

        origin : int
            Origin of the top-left corner. i.e. count from 0 or 1.
            Normally, origin should be 0 when passing numpy indices, or 1 if
            passing values from FITS header or map attributes.
            See `~astropy.wcs.WCS.wcs_pix2world` for more information.

        Returns
        -------

        coord : `astropy.coordinates.SkyCoord`
            A coordinate object representing the output coordinate.

        """

        # Hold the WCS instance here so we can inspect the output units after
        # the pix2world call
        temp_wcs = self.wcs

        x, y = temp_wcs.wcs_pix2world(x, y, origin)

        out_units = list(map(u.Unit, temp_wcs.wcs.cunit))

        x = u.Quantity(x, out_units[0])
        y = u.Quantity(y, out_units[1])

        return SkyCoord(x, y, frame=self.coordinate_frame)

# #### I/O routines #### #

    def save(self, filepath, filetype='auto', **kwargs):
        """Saves the SunPy Map object to a file.

        Currently SunPy can only save files in the FITS format. In the future
        support will be added for saving to other formats.

        Parameters
        ----------
        filepath : str
            Location to save file to.
        filetype : str
            'auto' or any supported file extension.
        hdu_type: None, `~fits.CompImageHDU`
            `None` will return a normal FITS file.
            `~fits.CompImageHDU` will rice compress the FITS file.
        kwargs :
            Any additional keyword arguments are passed to
            `~sunpy.io.write_file`.
        """
        io.write_file(filepath, self.data, self.meta, filetype=filetype,
                      **kwargs)

# #### Image processing routines #### #

    @u.quantity_input
    def resample(self, dimensions: u.pixel, method='linear'):
        """Returns a new Map that has been resampled up or down

        Arbitrary resampling of the Map to new dimension sizes.

        Uses the same parameters and creates the same co-ordinate lookup points
        as IDL''s congrid routine, which apparently originally came from a
        VAX/VMS routine of the same name.

        Parameters
        ----------
        dimensions : `~astropy.units.Quantity`
            Pixel dimensions that new Map should have.
            Note: the first argument corresponds to the 'x' axis and the second
            argument corresponds to the 'y' axis.
        method : {'neighbor' | 'nearest' | 'linear' | 'spline'}
            Method to use for resampling interpolation.
                * neighbor - Closest value from original data
                * nearest and linear - Uses n x 1-D interpolations using
                  scipy.interpolate.interp1d
                * spline - Uses ndimage.map_coordinates

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new Map which has been resampled to the desired dimensions.

        References
        ----------
        * `Rebinning <https://scipy-cookbook.readthedocs.io/items/Rebinning.html>`_
        """

        # Note: because the underlying ndarray is transposed in sense when
        #   compared to the Map, the ndarray is transposed, resampled, then
        #   transposed back
        # Note: "center" defaults to True in this function because data
        #   coordinates in a Map are at pixel centers

        # Make a copy of the original data and perform resample
        new_data = sunpy_image_resample(self.data.copy().T, dimensions,
                                        method, center=True)
        new_data = new_data.T

        scale_factor_x = float(self.dimensions[0] / dimensions[0])
        scale_factor_y = float(self.dimensions[1] / dimensions[1])

        # Update image scale and number of pixels
        new_meta = self.meta.copy()

        # Update metadata
        new_meta['cdelt1'] *= scale_factor_x
        new_meta['cdelt2'] *= scale_factor_y
        if 'CD1_1' in new_meta:
            new_meta['CD1_1'] *= scale_factor_x
            new_meta['CD2_1'] *= scale_factor_x
            new_meta['CD1_2'] *= scale_factor_y
            new_meta['CD2_2'] *= scale_factor_y
        new_meta['crpix1'] = (dimensions[0].value + 1) / 2.
        new_meta['crpix2'] = (dimensions[1].value + 1) / 2.
        lon, lat = self._get_lon_lat(self.center.frame)
        new_meta['crval1'] = lon.value
        new_meta['crval2'] = lat.value

        # Create new map instance
        new_map = self._new_instance(new_data, new_meta, self.plot_settings)
        return new_map

    @u.quantity_input
    def rotate(self, angle: u.deg = None, rmatrix=None, order=4, scale=1.0,
               recenter=False, missing=0.0, use_scipy=False):
        """
        Returns a new rotated and rescaled map.

        Specify either a rotation angle or a rotation matrix, but not both. If
        neither an angle or a rotation matrix are specified, the map will be
        rotated by the rotation angle in the metadata.

        The map will be rotated around the reference coordinate defined in the
        meta data.

        This method also updates the ``rotation_matrix`` attribute and any
        appropriate header data so that they correctly describe the new map.

        Parameters
        ----------
        angle : `~astropy.units.Quantity`
            The angle (degrees) to rotate counterclockwise.
        rmatrix : 2x2
            Linear transformation rotation matrix.
        order : int 0-5
            Interpolation order to be used. When using scikit-image this
            parameter is passed into :func:`skimage.transform.warp` (e.g., 4
            corresponds to bi-quartic interpolation).
            When using scipy it is passed into
            :func:`scipy.ndimage.interpolation.affine_transform` where it
            controls the order of the spline. Faster performance may be
            obtained at the cost of accuracy by using lower values.
            Default: 4
        scale : float
            A scale factor for the image, default is no scaling
        recenter : bool
            If True, position the axis of rotation at the center of the new map
            Default: False
        missing : float
            The numerical value to fill any missing points after rotation.
            Default: 0.0
        use_scipy : bool
            If True, forces the rotation to use
            :func:`scipy.ndimage.interpolation.affine_transform`, otherwise it
            uses the :func:`skimage.transform.warp`.
            Default: False, unless scikit-image can't be imported

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new Map instance containing the rotated and rescaled data of the
            original map.

        See Also
        --------
        sunpy.image.transform.affine_transform : The routine this method calls
        for the rotation.

        Notes
        -----
        This function will remove old CROTA keywords from the header.
        This function will also convert a CDi_j matrix to a PCi_j matrix.

        See :func:`sunpy.image.transform.affine_transform` for details on the
        transformations, situations when the underlying data is modified prior
        to rotation, and differences from IDL's rot().
        """
        # Put the import here to reduce sunpy.map import time
        from sunpy.image.transform import affine_transform

        if angle is not None and rmatrix is not None:
            raise ValueError("You cannot specify both an angle and a rotation matrix.")
        elif angle is None and rmatrix is None:
            rmatrix = self.rotation_matrix

        if order not in range(6):
            raise ValueError("Order must be between 0 and 5.")

        # The FITS-WCS transform is by definition defined around the
        # reference coordinate in the header.
        lon, lat = self._get_lon_lat(self.reference_coordinate.frame)
        rotation_center = u.Quantity([lon, lat])

        # Copy meta data
        new_meta = self.meta.copy()
        if angle is not None:
            # Calculate the parameters for the affine_transform
            c = np.cos(np.deg2rad(angle))
            s = np.sin(np.deg2rad(angle))
            rmatrix = np.array([[c, -s],
                                [s, c]])

        # Calculate the shape in pixels to contain all of the image data
        extent = np.max(np.abs(np.vstack((self.data.shape @ rmatrix,
                                          self.data.shape @ rmatrix.T))), axis=0)

        # Calculate the needed padding or unpadding
        diff = np.asarray(np.ceil((extent - self.data.shape) / 2), dtype=int).ravel()
        # Pad the image array
        pad_x = int(np.max((diff[1], 0)))
        pad_y = int(np.max((diff[0], 0)))

        new_data = np.pad(self.data,
                          ((pad_y, pad_y), (pad_x, pad_x)),
                          mode='constant',
                          constant_values=(missing, missing))
        new_meta['crpix1'] += pad_x
        new_meta['crpix2'] += pad_y

        # All of the following pixel calculations use a pixel origin of 0

        pixel_array_center = (np.flipud(new_data.shape) - 1) / 2.0

        # Create a temporary map so we can use it for the data to pixel calculation.
        temp_map = self._new_instance(new_data, new_meta, self.plot_settings)

        # Convert the axis of rotation from data coordinates to pixel coordinates
        pixel_rotation_center = u.Quantity(temp_map.world_to_pixel(self.reference_coordinate,
                                                                   origin=0)).value
        del temp_map

        if recenter:
            pixel_center = pixel_rotation_center
        else:
            pixel_center = pixel_array_center

        # Apply the rotation to the image data
        new_data = affine_transform(new_data.T,
                                    np.asarray(rmatrix),
                                    order=order, scale=scale,
                                    image_center=np.flipud(pixel_center),
                                    recenter=recenter, missing=missing,
                                    use_scipy=use_scipy).T

        if recenter:
            new_reference_pixel = pixel_array_center
        else:
            # Calculate new pixel coordinates for the rotation center
            new_reference_pixel = pixel_center + np.dot(rmatrix,
                                                        pixel_rotation_center - pixel_center)
            new_reference_pixel = np.array(new_reference_pixel).ravel()

        # Define the new reference_pixel
        new_meta['crval1'] = rotation_center[0].value
        new_meta['crval2'] = rotation_center[1].value
        new_meta['crpix1'] = new_reference_pixel[0] + 1  # FITS pixel origin is 1
        new_meta['crpix2'] = new_reference_pixel[1] + 1  # FITS pixel origin is 1

        # Unpad the array if necessary
        unpad_x = -np.min((diff[1], 0))
        if unpad_x > 0:
            new_data = new_data[:, unpad_x:-unpad_x]
            new_meta['crpix1'] -= unpad_x
        unpad_y = -np.min((diff[0], 0))
        if unpad_y > 0:
            new_data = new_data[unpad_y:-unpad_y, :]
            new_meta['crpix2'] -= unpad_y

        # Calculate the new rotation matrix to store in the header by
        # "subtracting" the rotation matrix used in the rotate from the old one
        # That being calculate the dot product of the old header data with the
        # inverse of the rotation matrix.
        pc_C = np.dot(self.rotation_matrix, np.linalg.inv(rmatrix))
        new_meta['PC1_1'] = pc_C[0, 0]
        new_meta['PC1_2'] = pc_C[0, 1]
        new_meta['PC2_1'] = pc_C[1, 0]
        new_meta['PC2_2'] = pc_C[1, 1]

        # Update pixel size if image has been scaled.
        if scale != 1.0:
            new_meta['cdelt1'] = (self.scale[0] / scale).value
            new_meta['cdelt2'] = (self.scale[1] / scale).value

        # Remove old CROTA kwargs because we have saved a new PCi_j matrix.
        new_meta.pop('CROTA1', None)
        new_meta.pop('CROTA2', None)
        # Remove CDi_j header
        new_meta.pop('CD1_1', None)
        new_meta.pop('CD1_2', None)
        new_meta.pop('CD2_1', None)
        new_meta.pop('CD2_2', None)

        # Create new map with the modification
        new_map = self._new_instance(new_data, new_meta, self.plot_settings)

        return new_map

    def submap(self, bottom_left, top_right=None):
        """
        Returns a submap of the map defined by the rectangle given by the
        ``[bottom_left, top_right]`` coordinates.

        Parameters
        ----------
        bottom_left : `astropy.units.Quantity` or `~astropy.coordinates.SkyCoord`
            The bottom_left coordinate of the rectangle. If a `SkyCoord` it can
            have shape ``(2,)`` and also define ``top_right``. If specifying
            pixel coordinates it must be given as an `~astropy.units.Quantity`
            object with units of `~astropy.units.pixel`.
        top_right : `astropy.units.Quantity` or `~astropy.coordinates.SkyCoord`
            The top_right coordinate of the rectangle. Can only be omitted if
            ``bottom_left`` has shape ``(2,)``.

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new map instance is returned representing to specified
            sub-region.

        Examples
        --------
        >>> import astropy.units as u
        >>> from astropy.coordinates import SkyCoord
        >>> import sunpy.map
        >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
        >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
        >>> bl = SkyCoord(-300*u.arcsec, -300*u.arcsec, frame=aia.coordinate_frame)  # doctest: +REMOTE_DATA
        >>> tr = SkyCoord(500*u.arcsec, 500*u.arcsec, frame=aia.coordinate_frame)  # doctest: +REMOTE_DATA
        >>> aia.submap(bl, tr)   # doctest: +REMOTE_DATA
        SunPy Map
        ---------
        Observatory:		 SDO
        Instrument:		 AIA 3
        Detector:		 AIA
        Measurement:		 171.0 Angstrom
        Wavelength:		 171.0 Angstrom
        Observation Date:	 2011-06-07 06:33:02
        Exposure Time:		 0.234256 s
        Dimension:		 [334. 334.] pix
        Coordinate System:	 helioprojective
        Scale:			 [2.402792 2.402792] arcsec / pix
        Reference Pixel:	 [127.5 126.5] pix
        Reference Coord:	 [3.22309951 1.38578135] arcsec
        array([[ 450.4546 ,  565.81494,  585.0416 , ..., 1178.3234 , 1005.28284,
                977.8161 ],
            [ 474.20004,  516.1865 ,  555.7032 , ..., 1024.9636 , 1010.1449 ,
                1010.1449 ],
            [ 548.1609 ,  620.9256 ,  620.9256 , ...,  933.8139 , 1074.4924 ,
                1108.4492 ],
                ...,
            [ 203.58617,  195.52335,  225.75891, ...,  612.7742 ,  580.52295,
                560.3659 ],
            [ 206.00058,  212.1806 ,  232.78065, ...,  650.96185,  622.12177,
                537.6615 ],
            [ 229.32516,  236.07002,  222.5803 , ...,  517.1058 ,  586.8026 ,
                591.2992 ]], dtype=float32)

        >>> aia.submap([0,0]*u.pixel, [5,5]*u.pixel)   # doctest: +REMOTE_DATA
        SunPy Map
        ---------
        Observatory:		 SDO
        Instrument:		 AIA 3
        Detector:		 AIA
        Measurement:		 171.0 Angstrom
        Wavelength:		 171.0 Angstrom
        Observation Date:	 2011-06-07 06:33:02
        Exposure Time:		 0.234256 s
        Dimension:		 [5. 5.] pix
        Coordinate System:	 helioprojective
        Scale:			 [2.402792 2.402792] arcsec / pix
        Reference Pixel:	 [512.5 512.5] pix
        Reference Coord:	 [3.22309951 1.38578135] arcsec
        array([[-95.92475   ,   7.076416  ,  -1.9656711 ,  -2.9485066 ,
                -0.98283553],
            [-96.97533   ,  -5.1167884 ,   0.        ,   0.        ,
                0.9746264 ],
            [-93.99607   ,   1.0189276 ,  -4.0757103 ,   2.0378551 ,
                -2.0378551 ],
            [-96.97533   ,  -8.040668  ,  -2.9238791 ,  -5.1167884 ,
                -0.9746264 ],
            [-95.92475   ,   6.028058  ,  -4.9797    ,  -1.0483578 ,
                -3.9313421 ]], dtype=float32)
        """

        if isinstance(bottom_left, (astropy.coordinates.SkyCoord,
                                    astropy.coordinates.BaseCoordinateFrame)):
            if not top_right:
                if bottom_left.shape[0] != 2:
                    raise ValueError("If top_right is not specified bottom_left must have length two.")
                else:
                    lon, lat = self._get_lon_lat(bottom_left)
                    top_right = u.Quantity([lon[1], lat[1]])
                    bottom_left = u.Quantity([lon[0], lat[0]])
            else:
                bottom_left = u.Quantity(self._get_lon_lat(bottom_left))
                top_right = u.Quantity(self._get_lon_lat(top_right))

            top_left = u.Quantity([bottom_left[0], top_right[1]])
            bottom_right = u.Quantity([top_right[0], bottom_left[1]])

            corners = u.Quantity([bottom_left, bottom_right, top_left, top_right])
            coord = SkyCoord(corners, frame=self.coordinate_frame)
            pixel_corners = self.world_to_pixel(coord)

            # Round the pixel values, we use floor+1 so that we always have at
            # least one pixel width of data.
            x_pixels = u.Quantity([np.min(pixel_corners.x), np.max(pixel_corners.x)]).value
            x_pixels[0] = np.ceil(x_pixels[0])
            x_pixels[1] = np.floor(x_pixels[1] + 1)
            y_pixels = u.Quantity([np.min(pixel_corners.y), np.max(pixel_corners.y)]).value
            y_pixels[0] = np.ceil(y_pixels[0])
            y_pixels[1] = np.floor(y_pixels[1] + 1)

        elif (isinstance(bottom_left, u.Quantity) and bottom_left.unit.is_equivalent(u.pix) and
              isinstance(top_right, u.Quantity) and top_right.unit.is_equivalent(u.pix)):
            x_pixels = u.Quantity([bottom_left[0], top_right[0]]).value
            y_pixels = u.Quantity([top_right[1], bottom_left[1]]).value

        else:
            raise ValueError("Invalid input, bottom_left and top_right must either be SkyCoord or Quantity in pixels.")

        # Sort the pixel values so we always slice in the correct direction
        x_pixels.sort()
        y_pixels.sort()

        x_pixels = np.array(x_pixels)
        y_pixels = np.array(y_pixels)

        # Clip pixel values to max of array, prevents negative
        # indexing
        x_pixels[np.less(x_pixels, 0)] = 0
        x_pixels[np.greater(x_pixels, self.data.shape[1])] = self.data.shape[1]

        y_pixels[np.less(y_pixels, 0)] = 0
        y_pixels[np.greater(y_pixels, self.data.shape[0])] = self.data.shape[0]

        # Get ndarray representation of submap
        xslice = slice(int(x_pixels[0]), int(x_pixels[1]))
        yslice = slice(int(y_pixels[0]), int(y_pixels[1]))
        new_data = self.data[yslice, xslice].copy()

        # Make a copy of the header with updated centering information
        new_meta = self.meta.copy()
        new_meta['crpix1'] = self.reference_pixel.x.value - x_pixels[0]
        new_meta['crpix2'] = self.reference_pixel.y.value - y_pixels[0]
        new_meta['naxis1'] = new_data.shape[1]
        new_meta['naxis2'] = new_data.shape[0]

        # Create new map instance
        if self.mask is not None:
            new_mask = self.mask[yslice, xslice].copy()
            # Create new map with the modification
            new_map = self._new_instance(new_data, new_meta, self.plot_settings, mask=new_mask)
            return new_map
        # Create new map with the modification
        new_map = self._new_instance(new_data, new_meta, self.plot_settings)
        return new_map

    @u.quantity_input
    def superpixel(self, dimensions: u.pixel, offset: u.pixel=(0, 0)*u.pixel, func=np.sum):
        """Returns a new map consisting of superpixels formed by applying
        'func' to the original map data.

        Parameters
        ----------
        dimensions : tuple
            One superpixel in the new map is equal to (dimension[0],
            dimension[1]) pixels of the original map.
            Note: the first argument corresponds to the 'x' axis and the second
            argument corresponds to the 'y' axis.
        offset : tuple
            Offset from (0,0) in original map pixels used to calculate where
            the data used to make the resulting superpixel map starts.
        func : function applied to the original data
            The function 'func' must take a numpy array as its first argument,
            and support the axis keyword with the meaning of a numpy axis
            keyword (see the description of `~numpy.sum` for an example.)
            The default value of 'func' is `~numpy.sum`; using this causes
            superpixel to sum over (dimension[0], dimension[1]) pixels of the
            original map.

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new Map which has superpixels of the required size.

        References
        ----------
        | `Summarizing blocks of an array using a moving window <https://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html>`_
        """

        # Note: because the underlying ndarray is transposed in sense when
        #   compared to the Map, the ndarray is transposed, resampled, then
        #   transposed back.
        # Note: "center" defaults to True in this function because data
        #   coordinates in a Map are at pixel centers.

        if (offset.value[0] < 0) or (offset.value[1] < 0):
            raise ValueError("Offset is strictly non-negative.")

        # Make a copy of the original data, perform reshaping, and apply the
        # function.
        if self.mask is not None:
            reshaped = reshape_image_to_4d_superpixel(np.ma.array(self.data.copy(), mask=self.mask),
                                                      [dimensions.value[1], dimensions.value[0]],
                                                      [offset.value[1], offset.value[0]])
        else:
            reshaped = reshape_image_to_4d_superpixel(self.data.copy(),
                                                      [dimensions.value[1], dimensions.value[0]],
                                                      [offset.value[1], offset.value[0]])
        new_array = func(func(reshaped, axis=3), axis=1)

        # Update image scale and number of pixels

        # create copy of new meta data
        new_meta = self.meta.copy()

        new_nx = new_array.shape[1]
        new_ny = new_array.shape[0]

        # Update metadata
        new_meta['cdelt1'] = (dimensions[0] * self.scale[0]).value
        new_meta['cdelt2'] = (dimensions[1] * self.scale[1]).value
        if 'CD1_1' in new_meta:
            new_meta['CD1_1'] *= dimensions[0].value
            new_meta['CD2_1'] *= dimensions[0].value
            new_meta['CD1_2'] *= dimensions[1].value
            new_meta['CD2_2'] *= dimensions[1].value
        new_meta['crpix1'] = (new_nx + 1) / 2.
        new_meta['crpix2'] = (new_ny + 1) / 2.
        lon, lat = self._get_lon_lat(self.center.frame)
        new_meta['crval1'] = lon.to(self.spatial_units[0]).value + 0.5*(offset[0]*self.scale[0]).to(self.spatial_units[0]).value
        new_meta['crval2'] = lat.to(self.spatial_units[1]).value + 0.5*(offset[1]*self.scale[1]).to(self.spatial_units[1]).value

        # Create new map instance
        if self.mask is not None:
            new_data = np.ma.getdata(new_array)
            new_mask = np.ma.getmask(new_array)
        else:
            new_data = new_array
            new_mask = None

        # Create new map with the modified data
        new_map = self._new_instance(new_data, new_meta, self.plot_settings, mask=new_mask)
        return new_map

# #### Visualization #### #

    @property
    def cmap(self):
        """
        Return the `matplotlib.colors.Colormap` instance this map uses.
        """
        cmap = self.plot_settings['cmap']
        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)
            # Set the colormap to be this specific instance so we are not
            # returning a copy
            self.plot_settings['cmap'] = cmap
        return cmap

    @u.quantity_input
    def draw_grid(self, axes=None, grid_spacing: u.deg = 15*u.deg, annotate=True, **kwargs):
        """
        Draws a coordinate overlay on the plot in the Heliographic Stonyhurst
        coordinate system.

        To overlay other coordinate systems see the `WCSAxes Documentation
        <https://docs.astropy.org/en/stable/visualization/wcsaxes/overlaying_coordinate_systems.html>`_

        Parameters
        ----------
        axes: `~matplotlib.axes` or None
        Axes to plot limb on or None to use current axes.

        grid_spacing: `~astropy.units.Quantity`
            Spacing for longitude and latitude grid, if length two it specifies
            (lon, lat) spacing.

        annotate : `bool`
            Passing `False` disables the axes labels and the ticks on the top and right axes.

        Returns
        -------
        overlay: `~astropy.visualization.wcsaxes.coordinates_map.CoordinatesMap`
            The wcsaxes coordinate overlay instance.

        Notes
        -----
        Keyword arguments are passed onto the `sunpy.visualization.wcsaxes_compat.wcsaxes_heliographic_overlay` function.
        """

        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)
        if not wcsaxes_compat.is_wcsaxes(axes):
            raise TypeError("Overlay grids can only be plotted on WCSAxes plots.")
        return wcsaxes_compat.wcsaxes_heliographic_overlay(axes,
                                                           grid_spacing=grid_spacing,
                                                           annotate=annotate,
                                                           **kwargs)

    def draw_limb(self, axes=None, **kwargs):
        """
        Draws a circle representing the solar limb

        Parameters
        ----------
        axes: `~matplotlib.axes` or None
            Axes to plot limb on or None to use current axes.

        Returns
        -------
        circ: list
            A list containing the `~matplotlib.patches.Circle` object that
            has been added to the axes.

        Notes
        -----
        Keyword arguments are passed onto `matplotlib.patches.Circle`.
        """
        # Put import here to reduce sunpy.map import time
        from matplotlib import patches

        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        transform = wcsaxes_compat.get_world_transform(axes)
        if wcsaxes_compat.is_wcsaxes(axes):
            radius = self.rsun_obs.to(u.deg).value
        else:
            radius = self.rsun_obs.value
        c_kw = {'radius': radius,
                'fill': False,
                'color': 'white',
                'zorder': 100,
                'transform': transform
                }
        c_kw.update(kwargs)

        circ = patches.Circle([0, 0], **c_kw)
        axes.add_artist(circ)

        return [circ]

    @u.quantity_input
    def draw_rectangle(self, bottom_left, width: u.deg, height: u.deg, axes=None, **kwargs):
        """
        Draw a rectangle defined in world coordinates on the plot.

        Parameters
        ----------

        bottom_left : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The bottom left corner of the rectangle.

        width : `astropy.units.Quantity`
            The width of the rectangle.

        height : `astropy.units.Quantity`
            The height of the rectangle.

        axes : `matplotlib.axes.Axes`
            The axes on which to plot the rectangle, defaults to the current
            axes.

        Returns
        -------

        rect : `list`
            A list containing the `~matplotlib.patches.Rectangle` object, after
            it has been added to ``axes``.

        Notes
        -----

        Extra keyword arguments to this function are passed through to the
        `~matplotlib.patches.Rectangle` instance.

        """

        if not axes:
            axes = plt.gca()

        if wcsaxes_compat.is_wcsaxes(axes):
            axes_unit = u.deg
        else:
            axes_unit = self.spatial_units[0]

        coord = bottom_left.transform_to(self.coordinate_frame)
        bottom_left = u.Quantity((coord.data.lon, coord.data.lat),
                                 unit=axes_unit).value

        width = width.to(axes_unit).value
        height = height.to(axes_unit).value

        kwergs = {'transform': wcsaxes_compat.get_world_transform(axes),
                  'color': 'white',
                  'fill': False}
        kwergs.update(kwargs)

        rect = plt.Rectangle(bottom_left, width, height, **kwergs)

        axes.add_artist(rect)

        return [rect]

    @u.quantity_input
    def draw_contours(self, levels: u.percent, axes=None, **contour_args):
        """
        Draw contours of the data.

        Parameters
        ----------
        levels : `~astropy.units.Quantity`
            A list of numbers indicating the level curves to draw given in
            percent.

        axes : `matplotlib.axes.Axes`
            The axes on which to plot the rectangle, defaults to the current
            axes.

        Returns
        -------
        cs : `list`
            The `~matplotlib.QuadContourSet` object, after it has been added to
            ``axes``.

        Notes
        -----
        Extra keyword arguments to this function are passed through to the
        `~matplotlib.pyplot.contour` function.

        """
        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        # TODO: allow for use of direct input of contours but requires units of
        # map flux which is not yet implemented

        cs = axes.contour(self.data, 0.01 * levels.to('percent').value * self.data.max(),
                          **contour_args)
        return cs

    @peek_show
    def peek(self, draw_limb=False, draw_grid=False,
             colorbar=True, **matplot_args):
        """
        Displays a graphical overview of the data in this object for user evaluation.
        For the creation of plots, users should instead use the `~sunpy.map.GenericMap.plot`
        method and Matplotlib's pyplot framework.

        Parameters
        ----------
        draw_limb : bool
            Whether the solar limb should be plotted.

        draw_grid : bool or `~astropy.units.Quantity`
            Whether solar meridians and parallels are plotted.
            If `~astropy.units.Quantity` then sets degree difference between
            parallels and meridians.
        colorbar : bool
            Whether to display a colorbar next to the plot.
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting.
        """
        figure = plt.figure()
        axes = wcsaxes_compat.gca_wcs(self.wcs)

        im = self.plot(axes=axes, **matplot_args)

        if colorbar:
            if draw_grid:
                pad = 0.12  # Pad to compensate for ticks and axes labels
            else:
                pad = 0.05  # Default value for vertical colorbar
            figure.colorbar(im, pad=pad)

        if draw_limb:
            self.draw_limb(axes=axes)

        if isinstance(draw_grid, bool):
            if draw_grid:
                self.draw_grid(axes=axes)
        elif isinstance(draw_grid, u.Quantity):
            self.draw_grid(axes=axes, grid_spacing=draw_grid)
        else:
            raise TypeError("draw_grid should be a bool or an astropy Quantity.")

        return figure

    @u.quantity_input
    def plot(self, annotate=True, axes=None, title=True,
             clip_interval: u.percent = None, **imshow_kwargs):
        """
        Plots the map object using matplotlib, in a method equivalent
        to plt.imshow() using nearest neighbour interpolation.

        Parameters
        ----------
        annotate : `bool`, optional
            If `True`, the data is plotted at its natural scale; with
            title and axis labels.

        axes: `~matplotlib.axes` or None
            If provided the image will be plotted on the given axes. Else the
            current matplotlib axes will be used.

        title : `bool`, optional
            If `True`, include the title.

        clip_interval : two-element `~astropy.units.Quantity`, optional
            If provided, the data will be clipped to the percentile interval bounded by the two
            numbers.

        **imshow_kwargs  : `dict`
            Any additional imshow arguments that should be used
            when plotting.

        Examples
        --------
        >>> # Simple Plot with color bar
        >>> aia.plot()   # doctest: +SKIP
        >>> plt.colorbar()   # doctest: +SKIP

        >>> # Add a limb line and grid
        >>> aia.plot()   # doctest: +SKIP
        >>> aia.draw_limb()   # doctest: +SKIP
        >>> aia.draw_grid()   # doctest: +SKIP

        """
        # Get current axes
        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        if not wcsaxes_compat.is_wcsaxes(axes):
            warnings.warn("WCSAxes not being used as the axes object for this plot."
                          " Plots may have unexpected behaviour. To fix this pass "
                          "'projection=map' when creating the axes",
                          SunpyUserWarning)
            # Check if the image is properly oriented
            if not np.array_equal(self.rotation_matrix, np.identity(2)):
                warnings.warn("The axes of this map are not aligned to the pixel grid. Plot axes may be incorrect.",
                              SunpyUserWarning)

        # Normal plot
        imshow_args = copy.deepcopy(self.plot_settings)
        if 'title' in imshow_args:
            plot_settings_title = imshow_args.pop('title')
        else:
            plot_settings_title = self.latex_name

        if annotate:
            if title is True:
                title = plot_settings_title

            if title:
                axes.set_title(title)

            axes.set_xlabel(axis_labels_from_ctype(self.coordinate_system[0],
                                                   self.spatial_units[0]))
            axes.set_ylabel(axis_labels_from_ctype(self.coordinate_system[1],
                                                   self.spatial_units[1]))

        if not wcsaxes_compat.is_wcsaxes(axes):
            bl = self._get_lon_lat(self.bottom_left_coord)
            tr = self._get_lon_lat(self.top_right_coord)
            x_range = list(u.Quantity([bl[0], tr[0]]).to(self.spatial_units[0]).value)
            y_range = list(u.Quantity([bl[1], tr[1]]).to(self.spatial_units[1]).value)
            imshow_args.update({'extent': x_range + y_range})
        imshow_args.update(imshow_kwargs)

        if clip_interval is not None:
            if len(clip_interval) == 2:
                clip_percentages = clip_interval.to('%').value
                vmin, vmax = AsymmetricPercentileInterval(*clip_percentages).get_limits(self.data)
            else:
                raise ValueError("Clip percentile interval must be specified as two numbers.")

            imshow_args['vmin'] = vmin
            imshow_args['vmax'] = vmax

        if self.mask is None:
            ret = axes.imshow(self.data, **imshow_args)
        else:
            ret = axes.imshow(np.ma.array(np.asarray(self.data), mask=self.mask), **imshow_args)

        if wcsaxes_compat.is_wcsaxes(axes):
            wcsaxes_compat.default_wcs_grid(axes)

        # Set current image (makes colorbar work)
        plt.sca(axes)
        plt.sci(ret)

        return ret


class InvalidHeaderInformation(ValueError):
    """Exception to raise when an invalid header tag value is encountered for a
    FITS/JPEG 2000 file."""
    pass
