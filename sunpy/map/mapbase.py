"""
Map is a generic Map class from which all other Map classes inherit from.
"""
from __future__ import absolute_import, division, print_function
from sunpy.extern.six.moves import range

#pylint: disable=E1101,E1121,W0404,W0613
__authors__ = ["Russell Hewett, Stuart Mumford, Keith Hughitt, Steven Christe", "Jack Ireland"]
__email__ = "stuart@mumford.me.uk"

import warnings
import inspect
from abc import ABCMeta
from copy import deepcopy
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches, cm, colors

import astropy.wcs
from astropy.coordinates import Longitude, Latitude

from sunpy.image.transform import affine_transform
from .nddata_compat import NDDataCompat as NDData

import sunpy.io as io
import sunpy.wcs as wcs
import sunpy.coordinates
from sunpy.visualization import toggle_pylab, wcsaxes_compat
from sunpy.sun import constants
from sunpy.sun import sun
from sunpy.time import parse_time, is_time
from sunpy.image.rescale import reshape_image_to_4d_superpixel
from sunpy.image.rescale import resample as sunpy_image_resample

from sunpy.extern import six

import astropy.units as u

from collections import namedtuple
Pair = namedtuple('Pair', 'x y')


from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

__all__ = ['GenericMap']


"""
Questions
---------
* Should we use Helioviewer or VSO's data model? (e.g. map.meas, map.wavelength
or something else?)
* Should 'center' be renamed to 'offset' and crpix1 & 2 be used for 'center'?
"""

# GenericMap subclass registry.
MAP_CLASSES = OrderedDict()


class GenericMapMeta(ABCMeta):
    """
    Registration metaclass for `~sunpy.map.GenericMap`.

    This class checks for the existance of a method named ``is_datasource_for``
    when a subclass of `GenericMap` is defined. If it exists it will add that
    class to the registry.
    """

    _registry = MAP_CLASSES

    def __new__(mcls, name, bases, members):
        cls = super(GenericMapMeta, mcls).__new__(mcls, name, bases, members)

        # The registry contains the class as the key and the validation method
        # as the item.
        if 'is_datasource_for' in members:
            mcls._registry[cls] = cls.is_datasource_for

        return cls


@six.add_metaclass(GenericMapMeta)
class GenericMap(NDData):
    """
    A Generic spatially-aware 2D data array

    Parameters
    ----------
    data : `~numpy.ndarray`, list
        A 2d list or ndarray containing the map data
    meta : dict
        A dictionary of the original image header tags

    Examples
    --------
    >>> import sunpy.map
    >>> import sunpy.data
    >>> sunpy.data.download_sample_data(overwrite=False)   # doctest: +SKIP
    >>> import sunpy.data.sample
    >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> aia   # doctest: +NORMALIZE_WHITESPACE
    SunPy AIAMap
    ---------
    Observatory:         SDO
    Instrument:  AIA 3
    Detector:    AIA
    Measurement:         171.0 Angstrom
    Wavelength:  171.0 Angstrom
    Obs Date:    2011-03-19 10:54:00
    dt:          1.999601 s
    Dimension:   [ 1024.  1024.] pix
    scale:               [ 2.4  2.4] arcsec / pix
    <BLANKLINE>
    array([[ 0.3125, -0.0625, -0.125 , ...,  0.625 , -0.625 ,  0.    ],
           [ 1.    ,  0.1875, -0.8125, ...,  0.625 , -0.625 ,  0.    ],
           [-1.1875,  0.375 , -0.5   , ..., -0.125 , -0.625 , -1.1875],
           ...,
           [-0.625 ,  0.0625, -0.3125, ...,  0.125 ,  0.125 ,  0.125 ],
           [ 0.5625,  0.0625,  0.5625, ..., -0.0625, -0.0625,  0.    ],
           [ 0.5   , -0.125 ,  0.4375, ...,  0.6875,  0.6875,  0.6875]])


    >>> aia.units
    Pair(x=Unit("arcsec"), y=Unit("arcsec"))
    >>> aia.peek()   # doctest: +SKIP

    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    | http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    | http://docs.scipy.org/doc/numpy/reference/ufuncs.html
    | http://www.scipy.org/Subclasses

    Notes
    -----

    A number of the properties of this class are returned as two-value named
    tuples that can either be indexed by position ([0] or [1]) or be accessed by
    name (.x or .y).  The names "x" and "y" here refer to the first and second
    axes of the map, and may not necessarily correspond to any similarly named
    axes in the coordinate system.

    This class makes some assumptions about the WCS information contained in
    the meta data. The first and most extensive assumption is that it is
    FITS-like WCS information as defined in the FITS WCS papers.

    Within this scope it also makes some other assumptions.

    * In the case of APIS convention headers where the CROTAi/j arguments are
      provided it assumes that these can be converted to the standard PCi_j
      notation using equations 32 in Thompson (2006).

    * If a CDi_j matrix is provided it is assumed that it can be converted to a
      PCi_j matrix and CDELT keywords as described in Greisen & Calabretta (2002).

    * The 'standard' FITS keywords that are used by this class are the PCi_j
      matrix and CDELT, along with the other keywords specified in the WCS papers.
      All subclasses of this class must convert their header information to
      this formalism. The CROTA to PCi_j conversion is done in this class.

    .. warning::
        This class currently assumes that a header with the CDi_j matrix
        information also includes the CDELT keywords, without these keywords
        this class will not process the WCS information. This will be fixed.
        Also the rotation_matrix does not work if the CDELT1 and CDELT2
        keywords are exactly equal.
    """

    def __init__(self, data, header, **kwargs):

        super(GenericMap, self).__init__(data, meta=header, **kwargs)

        # Correct possibly missing meta keywords
        self._fix_date()
        self._fix_naxis()

        # Setup some attributes
        self._nickname = self.detector

        # Validate header
        # TODO: This should be a function of the header, not of the map
        self._validate_meta()
        self._shift = Pair(0 * u.arcsec, 0 * u.arcsec)

        if self.dtype == np.uint8:
            norm = None
        else:
            norm = colors.Normalize()
        # Visualization attributes
        self.plot_settings = {'cmap': cm.gray,
                              'norm': norm,
                              'interpolation': 'nearest',
                              'origin': 'lower'
                              }

    def __getitem__(self, key):
        """ This should allow indexing by physical coordinate """
        raise NotImplementedError(
            "The ability to index Map by physical coordinate is not yet implemented.")

    def __repr__(self):
        if not self.observatory:
            return self.data.__repr__()
        return (
"""SunPy {dtype!s}
---------
Observatory:\t {obs}
Instrument:\t {inst}
Detector:\t {det}
Measurement:\t {meas}
Wavelength:\t {wave}
Obs Date:\t {date:{tmf}}
dt:\t\t {dt:f}
Dimension:\t {dim}
scale:\t\t {scale}

""".format(dtype=self.__class__.__name__,
           obs=self.observatory, inst=self.instrument, det=self.detector,
           meas=self.measurement, wave=self.wavelength, date=self.date, dt=self.exposure_time,
           dim=u.Quantity(self.dimensions),
           scale=u.Quantity(self.scale),
           tmf=TIME_FORMAT)
+ self.data.__repr__())

    @property
    def wcs(self):
        """
        The `~astropy.wcs.WCS` property of the map.
        """
        w2 = astropy.wcs.WCS(naxis=2)
        w2.wcs.crpix = u.Quantity(self.reference_pixel)
        # Make these a quantity array to prevent the numpy setting element of
        # array with sequence error.
        w2.wcs.cdelt = u.Quantity(self.scale)
        w2.wcs.crval = u.Quantity(self.reference_coordinate)
        w2.wcs.ctype = self.coordinate_system
        w2.wcs.pc = self.rotation_matrix
        w2.wcs.cunit = self.spatial_units
        w2.wcs.dateobs = self.date.isoformat()
        w2.heliographic_latitude = self.heliographic_latitude
        w2.heliographic_longitude = self.heliographic_longitude
        w2.dsun = self.dsun

        # Astropy WCS does not understand the SOHO default of "solar-x" and
        # "solar-y" ctypes.  This overrides the default assignment and
        # changes it to a ctype that is understood.  See Thompson, 2006, A.&A.,
        # 449, 791.
        if w2.wcs.ctype[0].lower() in ("solar-x", "solar_x"):
            w2.wcs.ctype[0] = 'HPLN-TAN'

        if w2.wcs.ctype[1].lower() in ("solar-y", "solar_y"):
            w2.wcs.ctype[1] = 'HPLT-TAN'

        return w2

    # Some numpy extraction
    @property
    def dimensions(self):
        """
        The dimensions of the array (x axis first, y axis second).
        """
        return Pair(*u.Quantity(np.flipud(self.data.shape), 'pixel'))

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

    @property
    def name(self):
        """Human-readable description of map-type"""
        return "{obs} {detector} {measurement} {date:{tmf}}".format(obs=self.observatory,
                                                                detector=self.detector,
                                                                measurement=self.measurement,
                                                                date=parse_time(self.date),
                                                                tmf=TIME_FORMAT)

    @property
    def nickname(self):
        """An abbreviated human-readable description of the map-type; part of
        the Helioviewer data model"""
        return self._nickname

    @nickname.setter
    def nickname(self, n):
        self._nickname = n

    @property
    def date(self):
        """Image observation time"""
        time = parse_time(self.meta.get('date-obs', 'now'))
        if time is None:
            warnings.warn_explicit("Missing metadata for observation time. Using current time.",
                                       Warning, __file__, inspect.currentframe().f_back.f_lineno)
        return parse_time(time)

#    @date.setter
#    def date(self, new_date):
#        self.meta['date-obs'] = new_date
#        #propagate change to malformed FITS keywords
#        if is_time(self.meta.get('date_obs', None)):
#            self.meta['date_obs'] = new_date

    @property
    def detector(self):
        """Detector name"""
        return self.meta.get('detector', "")

    @property
    def dsun(self):
        """The observer distance from the Sun."""
        dsun = self.meta.get('dsun_obs', None)

        if dsun is None:
            warnings.warn_explicit("Missing metadata for Sun-spacecraft separation: assuming Sun-Earth distance",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
            dsun = sun.sunearth_distance(self.date).to(u.m)

        return u.Quantity(dsun, 'm')

    @property
    def exposure_time(self):
        """Exposure time of the image in seconds."""
        return self.meta.get('exptime', 0.0) * u.s

    @property
    def instrument(self):
        """Instrument name"""
        return self.meta.get('instrume', "").replace("_", " ")

    @property
    def measurement(self):
        """Measurement name, defaults to the wavelength of image"""
        return u.Quantity(self.meta.get('wavelnth', 0), self.meta.get('waveunit', ""))

    @property
    def wavelength(self):
        """wavelength of the observation"""
        return u.Quantity(self.meta.get('wavelnth', 0), self.meta.get('waveunit', ""))

    @property
    def observatory(self):
        """Observatory or Telescope name"""
        return self.meta.get('obsrvtry', self.meta.get('telescop', "")).replace("_", " ")

    @property
    def xrange(self):
        """Return the X range of the image from edge to edge."""
        #TODO: This should be reading from the WCS object
        xmin = self.center.x - self.dimensions[0] / 2. * self.scale.x
        xmax = self.center.x + self.dimensions[0] / 2. * self.scale.x
        return u.Quantity([xmin, xmax])

    @property
    def yrange(self):
        """Return the Y range of the image from edge to edge."""
        #TODO: This should be reading from the WCS object
        ymin = self.center.y - self.dimensions[1] / 2. * self.scale.y
        ymax = self.center.y + self.dimensions[1] / 2. * self.scale.y
        return u.Quantity([ymin, ymax])

    @property
    def center(self):
        """The offset between the center of the Sun and the center of the map."""
        return Pair(wcs.get_center(self.dimensions[0], self.scale.x,
                                   self.reference_pixel.x,
                                   self.reference_coordinate.x),
                    wcs.get_center(self.dimensions[1], self.scale.y,
                                   self.reference_pixel.y,
                                   self.reference_coordinate.y))

    @property
    def shifted_value(self):
        """The total shift applied to the reference coordinate by past applications of
        `~sunpy.map.GenericMap.shift`."""
        return self._shift

    @u.quantity_input(x=u.deg, y=u.deg)
    def shift(self, x, y):
        """Returns a map shifted by a specified amount to, for example, correct for a bad
        map location. These values are applied directly to the `~sunpy.map.GenericMap.reference_coordinate`.
        To check how much shift has already been applied see `~sunpy.map.GenericMap.shifted_value`

        Parameters
        ----------

        x : `~astropy.units.Quantity`
            The shift to apply to the X coordinate.

        y : `~astropy.units.Quantity`
            The shift to apply to the Y coordinate

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new shifted Map.
        """
        new_map = deepcopy(self)
        new_map._shift = Pair(self.shifted_value.x + x, self.shifted_value.y + y)

        new_meta = self.meta.copy()

        # Update crvals
        new_meta['crval1'] = ((self.meta['crval1'] * self.spatial_units.x + x).to(self.spatial_units.x)).value
        new_meta['crval2'] = ((self.meta['crval2'] * self.spatial_units.y + y).to(self.spatial_units.y)).value

        new_map.meta = new_meta

        return new_map

    @property
    def rsun_meters(self):
        """Radius of the sun in meters"""
        return u.Quantity(self.meta.get('rsun_ref', constants.radius), 'meter')

    @property
    def rsun_obs(self):
        """Radius of the Sun."""
        rsun_arcseconds = self.meta.get('rsun_obs',
                                        self.meta.get('solar_r',
                                                      self.meta.get('radius', None)))

        if rsun_arcseconds is None:
            warnings.warn_explicit("Missing metadata for solar radius: assuming photospheric limb as seen from Earth",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
            rsun_arcseconds = sun.solar_semidiameter_angular_size(self.date).to('arcsec').value

        return u.Quantity(rsun_arcseconds, 'arcsec')

    @property
    def coordinate_system(self):
        """Coordinate system used for x and y axes (ctype1/2)"""
        return Pair(self.meta.get('ctype1', 'HPLN-TAN'),
                    self.meta.get('ctype2', 'HPLT-TAN'))

    @property
    def carrington_longitude(self):
        """Carrington longitude (crln_obs)"""
        carrington_longitude = self.meta.get('crln_obs', None)

        if carrington_longitude is None:
            warnings.warn_explicit("Missing metadata for Carrington longitude: assuming Earth-based observer",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
            carrington_longitude = (sun.heliographic_solar_center(self.date))[0]

        return u.Quantity(carrington_longitude, 'deg')

    @property
    def heliographic_latitude(self):
        """Heliographic latitude"""
        heliographic_latitude = self.meta.get('hglt_obs',
                                              self.meta.get('crlt_obs',
                                                            self.meta.get('solar_b0', None)))

        if heliographic_latitude is None:
            warnings.warn_explicit("Missing metadata for heliographic latitude: assuming Earth-based observer",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
            heliographic_latitude = (sun.heliographic_solar_center(self.date))[1]

        return u.Quantity(heliographic_latitude, 'deg')

    @property
    def heliographic_longitude(self):
        """Heliographic longitude"""
        return u.Quantity(self.meta.get('hgln_obs', 0.), 'deg')

    @property
    def reference_coordinate(self):
        """Reference point WCS axes in data units (i.e. crval1, crval2). This value
        includes a shift if one is set."""
        return Pair(self.meta.get('crval1', 0.) * self.spatial_units.x,
                    self.meta.get('crval2', 0.) * self.spatial_units.y)

    @property
    def reference_pixel(self):
        """Reference point axes in pixels (i.e. crpix1, crpix2)"""
        return Pair(self.meta.get('crpix1', (self.meta.get('naxis1') + 1) / 2.) * u.pixel,
                    self.meta.get('crpix2', (self.meta.get('naxis2') + 1) / 2.) * u.pixel)

    @property
    def scale(self):
        """Image scale along the x and y axes in units/pixel (i.e. cdelt1, cdelt2)"""
        #TODO: Fix this if only CDi_j matrix is provided
        return Pair(self.meta.get('cdelt1', 1.) * self.spatial_units.x / u.pixel,
                    self.meta.get('cdelt2', 1.) * self.spatial_units.y / u.pixel)

    @property
    def spatial_units(self):
        """Image coordinate units along the x and y axes (i.e. cunit1, cunit2)."""
        return Pair(u.Unit(self.meta.get('cunit1', 'arcsec')),
                    u.Unit(self.meta.get('cunit2', 'arcsec')))

    @property
    def rotation_matrix(self):
        """Matrix describing the rotation required to align solar North with
        the top of the image."""
        if 'PC1_1' in self.meta:
            return np.matrix([[self.meta['PC1_1'], self.meta['PC1_2']],
                              [self.meta['PC2_1'], self.meta['PC2_2']]])

        elif 'CD1_1' in self.meta:
            cd = np.matrix([[self.meta['CD1_1'], self.meta['CD1_2']],
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
        lam = self.scale.y / self.scale.x
        p = np.deg2rad(self.meta.get('CROTA2', 0))

        return np.matrix([[np.cos(p), -1 * lam * np.sin(p)],
                          [1/lam * np.sin(p), np.cos(p)]])

# #### Miscellaneous #### #

    def _fix_date(self):
        # Check commonly used but non-standard FITS keyword for observation time
        # and correct the keyword if we can.  Keep updating old one for
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
        cmap_string = self.observatory + self.meta['detector'] + str(int(self.wavelength.to('angstrom').value))
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

        warnings.simplefilter('always', Warning)

        for meta_property in ('cunit1', 'cunit2', 'waveunit'):
            if (self.meta.get(meta_property) and
                u.Unit(self.meta.get(meta_property),
                       parse_strict='silent').physical_type == 'unknown'):

                warnings.warn("Unknown value for "+meta_property.upper(), Warning)


# #### Data conversion routines #### #

    @u.quantity_input(x=u.deg, y=u.deg)
    def data_to_pixel(self, x, y, origin=0):
        """
        Convert a data (world) coordinate to a pixel coordinate by using
        `~astropy.wcs.WCS.wcs_world2pix`.

        Parameters
        ----------

        x : `~astropy.units.Quantity`
            Data coordinate of the CTYPE1 axis. (Normally solar-x).

        y : `~astropy.units.Quantity`
            Data coordinate of the CTYPE2 axis. (Normally solar-y).

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
        x, y = self.wcs.wcs_world2pix(x.to(u.deg).value, y.to(u.deg).value, origin)

        return x * u.pixel, y * u.pixel

    @u.quantity_input(x=u.pixel, y=u.pixel)
    def pixel_to_data(self, x, y, origin=0):
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

        x : `~astropy.units.Quantity`
            Coordinate of the CTYPE1 axis. (Normally solar-x).

        y : `~astropy.units.Quantity`
            Coordinate of the CTYPE2 axis. (Normally solar-y).
        """
        x, y = self.wcs.wcs_pix2world(x, y, origin)

        # If the wcs is celestial it is output in degress
        if self.wcs.is_celestial:
            x = u.Quantity(x, u.deg)
            y = u.Quantity(y, u.deg)
        else:
            x = u.Quantity(x, self.spatial_units.x)
            y = u.Quantity(y, self.spatial_units.y)

        x = Longitude(x, wrap_angle=180*u.deg)
        y = Latitude(y)

        return x.to(self.spatial_units.x), y.to(self.spatial_units.y)


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
            'auto' or any supported file extension
        """
        io.write_file(filepath, self.data, self.meta, filetype=filetype,
                      **kwargs)

# #### Image processing routines #### #

    @u.quantity_input(dimensions=u.pixel)
    def resample(self, dimensions, method='linear'):
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
        * `Rebinning <http://www.scipy.org/Cookbook/Rebinning>`_ (Original source, 2011/11/19)
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

        new_map = deepcopy(self)
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
        new_meta['crval1'] = self.center.x.value
        new_meta['crval2'] = self.center.y.value

        # Create new map instance
        new_map.data = new_data
        new_map.meta = new_meta
        return new_map

    def rotate(self, angle=None, rmatrix=None, order=4, scale=1.0,
               recenter=False, missing=0.0, use_scipy=False):
        """
        Returns a new rotated and rescaled map.  Specify either a rotation
        angle or a rotation matrix, but not both.  If neither an angle or a
        rotation matrix are specified, the map will be rotated by the rotation
        angle in the metadata.

        The map will be rotated around the reference coordinate defined in the
        meta data.

        Also updates the rotation_matrix attribute and any appropriate header
        data so that they correctly describe the new map.

        Parameters
        ----------
        angle : `~astropy.units.Quantity`
            The angle (degrees) to rotate counterclockwise.
        rmatrix : 2x2
            Linear transformation rotation matrix.
        order : int 0-5
            Interpolation order to be used. When using scikit-image this parameter
            is passed into :func:`skimage.transform.warp` (e.g., 4 corresponds to
            bi-quartic interpolation).
            When using scipy it is passed into
            :func:`scipy.ndimage.interpolation.affine_transform` where it controls
            the order of the spline.
            Faster performance may be obtained at the cost of accuracy by using lower values.
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
        sunpy.image.transform.affine_transform : The routine this method calls for the rotation.

        Notes
        -----
        This function will remove old CROTA keywords from the header.
        This function will also convert a CDi_j matrix to a PCi_j matrix.

        See :func:`sunpy.image.transform.affine_transform` for details on the
        transformations, situations when the underlying data is modified prior to rotation,
        and differences from IDL's rot().
        """
        if angle is not None and rmatrix is not None:
            raise ValueError("You cannot specify both an angle and a matrix")
        elif angle is None and rmatrix is None:
            rmatrix = self.rotation_matrix

        # This is out of the quantity_input decorator. To allow the angle=None
        # case. See https://github.com/astropy/astropy/issues/3734
        if angle:
            try:
                equivalent = angle.unit.is_equivalent(u.deg)

                if not equivalent:
                    raise u.UnitsError("Argument '{0}' to function '{1}'"
                                       " must be in units convertable to"
                                       " '{2}'.".format('angle', 'rotate',
                                                      u.deg.to_string()))

            # Either there is no .unit or no .is_equivalent
            except AttributeError:
                if hasattr(angle, "unit"):
                    error_msg = "a 'unit' attribute without an 'is_equivalent' method"
                else:
                    error_msg = "no 'unit' attribute"
                raise TypeError("Argument '{0}' to function '{1}' has {2}. "
                      "You may want to pass in an astropy Quantity instead."
                         .format('angle', 'rotate', error_msg))

        # Interpolation parameter sanity
        if order not in range(6):
            raise ValueError("Order must be between 0 and 5")

        # The FITS-WCS transform is by definition defined around the
        # reference coordinate in the header.
        rotation_center = u.Quantity([self.reference_coordinate.x,
                                      self.reference_coordinate.y])

        # Copy Map
        new_map = deepcopy(self)

        if angle is not None:
            # Calculate the parameters for the affine_transform
            c = np.cos(np.deg2rad(angle))
            s = np.sin(np.deg2rad(angle))
            rmatrix = np.matrix([[c, -s], [s, c]])

        # Calculate the shape in pixels to contain all of the image data
        extent = np.max(np.abs(np.vstack((new_map.data.shape * rmatrix,
                                          new_map.data.shape * rmatrix.T))), axis=0)
        # Calculate the needed padding or unpadding
        diff = np.asarray(np.ceil((extent - new_map.data.shape) / 2)).ravel()
        # Pad the image array
        pad_x = int(np.max((diff[1], 0)))
        pad_y = int(np.max((diff[0], 0)))
        new_map.data = np.pad(new_map.data,
                              ((pad_y, pad_y), (pad_x, pad_x)),
                              mode='constant',
                              constant_values=(missing, missing))
        new_map.meta['crpix1'] += pad_x
        new_map.meta['crpix2'] += pad_y

        # All of the following pixel calculations use a pixel origin of 0

        pixel_array_center = (np.flipud(new_map.data.shape) - 1) / 2.0

        # Convert the axis of rotation from data coordinates to pixel coordinates
        pixel_rotation_center = u.Quantity(new_map.data_to_pixel(*rotation_center,
                                                                 origin=0)).value
        if recenter:
            pixel_center = pixel_rotation_center
        else:
            pixel_center = pixel_array_center

        # Apply the rotation to the image data
        new_map.data = affine_transform(new_map.data.T,
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
        new_map.meta['crval1'] = rotation_center[0].value
        new_map.meta['crval2'] = rotation_center[1].value
        new_map.meta['crpix1'] = new_reference_pixel[0] + 1 # FITS pixel origin is 1
        new_map.meta['crpix2'] = new_reference_pixel[1] + 1 # FITS pixel origin is 1

        # Unpad the array if necessary
        unpad_x = -np.min((diff[1], 0))
        if unpad_x > 0:
            new_map.data = new_map.data[:, unpad_x:-unpad_x]
            new_map.meta['crpix1'] -= unpad_x
        unpad_y = -np.min((diff[0], 0))
        if unpad_y > 0:
            new_map.data = new_map.data[unpad_y:-unpad_y, :]
            new_map.meta['crpix2'] -= unpad_y

        # Calculate the new rotation matrix to store in the header by
        # "subtracting" the rotation matrix used in the rotate from the old one
        # That being calculate the dot product of the old header data with the
        # inverse of the rotation matrix.
        pc_C = np.dot(new_map.rotation_matrix, rmatrix.I)
        new_map.meta['PC1_1'] = pc_C[0,0]
        new_map.meta['PC1_2'] = pc_C[0,1]
        new_map.meta['PC2_1'] = pc_C[1,0]
        new_map.meta['PC2_2'] = pc_C[1,1]

        # Update pixel size if image has been scaled.
        if scale != 1.0:
            new_map.meta['cdelt1'] = (new_map.scale.x / scale).value
            new_map.meta['cdelt2'] = (new_map.scale.y / scale).value

        # Remove old CROTA kwargs because we have saved a new PCi_j matrix.
        new_map.meta.pop('CROTA1', None)
        new_map.meta.pop('CROTA2', None)
        # Remove CDi_j header
        new_map.meta.pop('CD1_1', None)
        new_map.meta.pop('CD1_2', None)
        new_map.meta.pop('CD2_1', None)
        new_map.meta.pop('CD2_2', None)

        return new_map

    def submap(self, range_a, range_b):
        """
        Returns a submap of the map with the specified range.

        Parameters
        ----------
        range_a : `astropy.units.Quantity`
            The range of the Map to select across either the x axis.
            Can be either in data units (normally arcseconds) or pixel units.
        range_b : `astropy.units.Quantity`
            The range of the Map to select across either the y axis.
            Can be either in data units (normally arcseconds) or pixel units.

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new map instance is returned representing to specified sub-region

        Examples
        --------
        >>> import astropy.units as u
        >>> import sunpy.map
        >>> import sunpy.data
        >>> sunpy.data.download_sample_data(overwrite=False)   # doctest: +SKIP
        >>> import sunpy.data.sample
        >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
        >>> aia.submap([-5,5]*u.arcsec, [-5,5]*u.arcsec)   # doctest: +NORMALIZE_WHITESPACE
        SunPy AIAMap
        ---------
        Observatory:         SDO
        Instrument:  AIA 3
        Detector:    AIA
        Measurement:         171.0 Angstrom
        Wavelength:  171.0 Angstrom
        Obs Date:    2011-03-19 10:54:00
        dt:          1.999601 s
        Dimension:   [ 4.  4.] pix
        scale:               [ 2.4  2.4] arcsec / pix
        <BLANKLINE>
        array([[ 273.4375,  247.4375,  303.5   ,  305.3125],
               [ 302.3125,  298.125 ,  299.    ,  261.5   ],
               [ 289.75  ,  269.25  ,  256.375 ,  242.3125],
               [ 241.75  ,  248.8125,  263.0625,  249.0625]])

        >>> aia.submap([0,5]*u.pixel, [0,5]*u.pixel)   # doctest: +NORMALIZE_WHITESPACE
        SunPy AIAMap
        ---------
        Observatory:         SDO
        Instrument:  AIA 3
        Detector:    AIA
        Measurement:         171.0 Angstrom
        Wavelength:  171.0 Angstrom
        Obs Date:    2011-03-19 10:54:00
        dt:          1.999601 s
        Dimension:   [ 5.  5.] pix
        scale:               [ 2.4  2.4] arcsec / pix
        <BLANKLINE>
        array([[ 0.3125, -0.0625, -0.125 ,  0.    , -0.375 ],
               [ 1.    ,  0.1875, -0.8125,  0.125 ,  0.3125],
               [-1.1875,  0.375 , -0.5   ,  0.25  , -0.4375],
               [-0.6875, -0.3125,  0.8125,  0.0625,  0.1875],
               [-0.875 ,  0.25  ,  0.1875,  0.    , -0.6875]])
        """

        # Do manual Quantity input validation to allow for two unit options
        if ((isinstance(range_a, u.Quantity) and isinstance(range_b, u.Quantity)) or
            (hasattr(range_a, 'unit') and hasattr(range_b, 'unit'))):

            if (range_a.unit.is_equivalent(self.spatial_units.x) and
                range_b.unit.is_equivalent(self.spatial_units.y)):
                units = 'data'
            elif range_a.unit.is_equivalent(u.pixel) and range_b.unit.is_equivalent(u.pixel):
                units = 'pixels'
            else:
                raise u.UnitsError("range_a and range_b but be "
                                   "in units convertable to {} or {}".format(self.spatial_units['x'],
                                                                             u.pixel))
        else:
            raise TypeError("Arguments range_a and range_b to function submap "
                            "have an invalid unit attribute "
                            "You may want to pass in an astropy Quantity instead.")

        if units is "data":
            # Check edges (e.g. [:512,..] or [:,...])
            if range_a[0] is None:
                range_a[0] = self.xrange[0]
            if range_a[1] is None:
                range_a[1] = self.xrange[1]
            if range_b[0] is None:
                range_b[0] = self.yrange[0]
            if range_b[1] is None:
                range_b[1] = self.yrange[1]

            x1, y1 = np.ceil(u.Quantity(self.data_to_pixel(range_a[0], range_b[0]))).value
            x2, y2 = np.floor(u.Quantity(self.data_to_pixel(range_a[1], range_b[1])) + 1*u.pix).value

            x_pixels = [x1, x2]
            y_pixels = [y1, y2]

        elif units is "pixels":
            # Check edges
            if range_a[0] is None:
                range_a[0] = 0
            if range_a[1] is None:
                range_a[1] = self.data.shape[1]
            if range_b[0] is None:
                range_b[0] = 0
            if range_b[1] is None:
                range_b[1] = self.data.shape[0]

            x_pixels = range_a.value
            y_pixels = range_b.value
        else:
            raise ValueError("Invalid unit. Must be one of 'data' or 'pixels'")

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
        new_map = deepcopy(self)
        new_map.meta['crpix1'] = self.reference_pixel.x.value - x_pixels[0]
        new_map.meta['crpix2'] = self.reference_pixel.y.value - y_pixels[0]
        new_map.meta['naxis1'] = new_data.shape[1]
        new_map.meta['naxis2'] = new_data.shape[0]

        # Create new map instance
        new_map.data = new_data
        if self.mask is not None:
            new_map.mask = self.mask[yslice, xslice].copy()

        return new_map

    @u.quantity_input(dimensions=u.pixel, offset=u.pixel)
    def superpixel(self, dimensions, offset=(0, 0)*u.pixel, func=np.sum):
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
        | `Summarizing blocks of an array using a moving window <http://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html>`_
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
        new_map = deepcopy(self)
        new_meta = new_map.meta

        new_nx = new_array.shape[1]
        new_ny = new_array.shape[0]

        # Update metadata
        new_meta['cdelt1'] = (dimensions[0] * self.scale.x).value
        new_meta['cdelt2'] = (dimensions[1] * self.scale.y).value
        if 'CD1_1' in new_meta:
            new_meta['CD1_1'] *= dimensions[0].value
            new_meta['CD2_1'] *= dimensions[0].value
            new_meta['CD1_2'] *= dimensions[1].value
            new_meta['CD2_2'] *= dimensions[1].value
        new_meta['crpix1'] = (new_nx + 1) / 2.
        new_meta['crpix2'] = (new_ny + 1) / 2.
        new_meta['crval1'] = self.center.x.to(self.spatial_units.x).value + 0.5*(offset[0]*self.scale.x).to(self.spatial_units.x).value
        new_meta['crval2'] = self.center.y.to(self.spatial_units.y).value + 0.5*(offset[1]*self.scale.y).to(self.spatial_units.y).value

        # Create new map instance
        if self.mask is not None:
            new_map.data = np.ma.getdata(new_array)
            new_map.mask = np.ma.getmask(new_array)
        else:
            new_map.data = new_array
            new_map.mask = None
        return new_map

# #### Visualization #### #

    @u.quantity_input(grid_spacing=u.deg)
    def draw_grid(self, axes=None, grid_spacing=15*u.deg, **kwargs):
        """Draws a grid over the surface of the Sun

        Parameters
        ----------
        axes: `~matplotlib.axes` or None
        Axes to plot limb on or None to use current axes.

        grid_spacing: float
            Spacing (in degrees) for longitude and latitude grid.

        Returns
        -------
        lines: list
            A list of `matplotlib.lines.Line2D` objects that have been plotted.

        Notes
        -----
        keyword arguments are passed onto matplotlib.pyplot.plot
        """

        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        lines = []

        # Do not automatically rescale axes when plotting the overlay
        axes.set_autoscale_on(False)

        transform = wcsaxes_compat.get_world_transform(axes)

        XX, YY = np.meshgrid(np.arange(self.data.shape[0]),
                             np.arange(self.data.shape[1]))
        x, y = self.pixel_to_data(XX*u.pix, YY*u.pix)
        dsun = self.dsun

        b0 = self.heliographic_latitude.to(u.deg).value
        l0 = self.heliographic_longitude.to(u.deg).value
        units = self.spatial_units

        # Prep the plot kwargs
        plot_kw = {'color': 'white',
                   'linestyle': 'dotted',
                   'zorder': 100,
                   'transform': transform}
        plot_kw.update(kwargs)

        hg_longitude_deg = np.linspace(-180, 180, num=361) + l0
        hg_latitude_deg = np.arange(-90, 90, grid_spacing.to(u.deg).value)

        # draw the latitude lines
        for lat in hg_latitude_deg:
            x, y = wcs.convert_hg_hpc(hg_longitude_deg, lat * np.ones(361),
                                      b0_deg=b0, l0_deg=l0, dsun_meters=dsun,
                                      angle_units=units.x, occultation=True)
            valid = np.logical_and(np.isfinite(x), np.isfinite(y))
            x = x[valid]
            y = y[valid]
            if wcsaxes_compat.is_wcsaxes(axes):
                x = (x*u.arcsec).to(u.deg).value
                y = (y*u.arcsec).to(u.deg).value
            lines += axes.plot(x, y, **plot_kw)

        hg_longitude_deg = np.arange(-180, 180, grid_spacing.to(u.deg).value) + l0
        hg_latitude_deg = np.linspace(-90, 90, num=181)

        # draw the longitude lines
        for lon in hg_longitude_deg:
            x, y = wcs.convert_hg_hpc(lon * np.ones(181), hg_latitude_deg,
                                      b0_deg=b0, l0_deg=l0, dsun_meters=dsun,
                                      angle_units=units[0], occultation=True)
            valid = np.logical_and(np.isfinite(x), np.isfinite(y))
            x = x[valid]
            y = y[valid]
            if wcsaxes_compat.is_wcsaxes(axes):
                x = (x*u.arcsec).to(u.deg).value
                y = (y*u.arcsec).to(u.deg).value
            lines += axes.plot(x, y, **plot_kw)

        # Turn autoscaling back on.
        axes.set_autoscale_on(True)
        return lines

    def draw_limb(self, axes=None, **kwargs):
        """Draws a circle representing the solar limb

            Parameters
            ----------
            axes: `~matplotlib.axes` or None
                Axes to plot limb on or None to use current axes.

            Returns
            -------
            circ: list
                A list containing the `matplotlib.patches.Circle` object that
                has been added to the axes.

            Notes
            -----
            keyword arguments are passed onto the Circle Patch, see:
            http://matplotlib.org/api/artist_api.html#matplotlib.patches.Patch
            http://matplotlib.org/api/artist_api.html#matplotlib.patches.Circle
        """

        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        transform = wcsaxes_compat.get_world_transform(axes)
        if wcsaxes_compat.is_wcsaxes(axes):
            radius = self.rsun_obs.to(u.deg).value
        else:
            radius = self.rsun_obs.value
        c_kw = {'radius':radius,
                'fill':False,
                'color':'white',
                'zorder':100,
                'transform': transform
                }
        c_kw.update(kwargs)

        circ = patches.Circle([0, 0], **c_kw)
        axes.add_artist(circ)

        return [circ]

    @u.quantity_input(bottom_left=u.deg, width=u.deg, height=u.deg)
    def draw_rectangle(self, bottom_left, width, height, axes=None, **kwargs):
        """
        Draw a rectangle defined in world coordinates on the plot.

        Parameters
        ----------

        bottom_left : `astropy.units.Quantity`
            The bottom left corner of the rectangle.

        width : `astropy.units.Quantity`
            The width of the rectangle.

        height : `astropy.units.Quantity`
            The height of the rectangle.

        axes : `matplotlib.axes.Axes`
            The axes on which to plot the rectangle, defaults to the current axes.

        Returns
        -------

        rect : `list`
            A list containing the `~matplotlib.patches.Rectangle` object, after it has been added to ``axes``.

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

        bottom_left = bottom_left.to(axes_unit).value
        width = width.to(axes_unit).value
        height = height.to(axes_unit).value

        kwergs = {'transform':wcsaxes_compat.get_world_transform(axes),
                  'color':'white',
                  'fill':False}
        kwergs.update(kwargs)
        rect = plt.Rectangle(bottom_left, width, height, **kwergs)

        axes.add_artist(rect)

        return [rect]

    @u.quantity_input(levels=u.percent)
    def draw_contours(self, levels, axes=None, **contour_args):
        """
        Draw contours of the data

        Parameters
        ----------

        levels : `~astropy.units.Quantity`
            A list of numbers indicating the level curves to draw given in percent.

        axes : `matplotlib.axes.Axes`
            The axes on which to plot the rectangle, defaults to the current axes.

        Returns
        -------

        cs : `list`
            The `~matplotlib.QuadContourSet` object, after it has been added to ``axes``.

        Notes
        -----

        Extra keyword arguments to this function are passed through to the
        `~matplotlib.pyplot.contour` function.

        """
        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        #TODO: allow for use of direct input of contours but requires units of map flux which is not yet implemented

        cs = axes.contour(self.data, 0.01 * levels.to('percent').value * self.data.max(), **contour_args)
        return cs

    @toggle_pylab
    def peek(self, draw_limb=False, draw_grid=False,
                   colorbar=True, basic_plot=False, **matplot_args):
        """Displays the map in a new figure

        Parameters
        ----------
        draw_limb : bool
            Whether the solar limb should be plotted.

        draw_grid : bool or `~astropy.units.Quantity`
            Whether solar meridians and parallels are plotted.
            If `~astropy.units.Quantity` then sets degree difference between
            parallels and meridians.
        gamma : float
            Gamma value to use for the color map
        colorbar : bool
            Whether to display a colorbar next to the plot
        basic_plot : bool
            If true, the data is plotted by itself at it's natural scale; no
            title, labels, or axes are shown.
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting.
        """

        # Create a figure and add title and axes
        figure = plt.figure(frameon=not basic_plot)

        # Basic plot
        if basic_plot:
            axes = plt.Axes(figure, [0., 0., 1., 1.])
            axes.set_axis_off()
            figure.add_axes(axes)
            matplot_args.update({'annotate':False})

        # Normal plot
        else:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        im = self.plot(axes=axes, **matplot_args)

        if colorbar and not basic_plot:
            figure.colorbar(im)

        if draw_limb:
            self.draw_limb(axes=axes)

        if isinstance(draw_grid, bool):
            if draw_grid:
                self.draw_grid(axes=axes)
        elif isinstance(draw_grid, u.Quantity):
            self.draw_grid(axes=axes, grid_spacing=draw_grid)
        else:
            raise TypeError("draw_grid should be a bool or an astropy Quantity.")

        figure.show()

    @toggle_pylab
    def plot(self, annotate=True, axes=None, title=True, **imshow_kwargs):
        """ Plots the map object using matplotlib, in a method equivalent
        to plt.imshow() using nearest neighbour interpolation.

        Parameters
        ----------
        annotate : bool
            If True, the data is plotted at it's natural scale; with
            title and axis labels.

        axes: `~matplotlib.axes` or None
            If provided the image will be plotted on the given axes. Else the
            current matplotlib axes will be used.

        **imshow_kwargs  : dict
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

        # Check that the image is properly oriented
        if (not wcsaxes_compat.is_wcsaxes(axes) and
            not np.array_equal(self.rotation_matrix, np.matrix(np.identity(2)))):
            warnings.warn("This map is not properly oriented. Plot axes may be incorrect",
                          Warning)

        # Normal plot
        imshow_args = deepcopy(self.plot_settings)
        if 'title' in imshow_args:
            plot_settings_title = imshow_args.pop('title')
        else:
            plot_settings_title = self.name

        if annotate:
            if title is True:
                title = plot_settings_title

            if title:
                axes.set_title(title)

            # x-axis label
            if self.coordinate_system.x == 'HG':
                xlabel = 'Longitude [{lon}]'.format(lon=self.spatial_units.x)
            else:
                xlabel = 'X-position [{xpos}]'.format(xpos=self.spatial_units.x)

            # y-axis label
            if self.coordinate_system.y == 'HG':
                ylabel = 'Latitude [{lat}]'.format(lat=self.spatial_units.y)
            else:
                ylabel = 'Y-position [{ypos}]'.format(ypos=self.spatial_units.y)

            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)

        if not wcsaxes_compat.is_wcsaxes(axes):
            imshow_args.update({'extent': list(self.xrange.value) + list(self.yrange.value)})
        imshow_args.update(imshow_kwargs)

        if self.mask is None:
            ret = axes.imshow(self.data, **imshow_args)
        else:
            ret = axes.imshow(np.ma.array(np.asarray(self.data), mask=self.mask), **imshow_args)

        if wcsaxes_compat.is_wcsaxes(axes):
            wcsaxes_compat.default_wcs_grid(axes)

        #Set current image (makes colorbar work)
        plt.sca(axes)
        plt.sci(ret)
        return ret


class InvalidHeaderInformation(ValueError):
    """Exception to raise when an invalid header tag value is encountered for a
    FITS/JPEG 2000 file."""
    pass
