"""
Map is a generic Map class from which all other Map classes inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404,W0613
__authors__ = ["Russell Hewett, Stuart Mumford, Keith Hughitt, Steven Christe"]
__email__ = "stuart@mumford.me.uk"

import warnings
import inspect
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import cm

import astropy.nddata
from sunpy.image.transform import affine_transform

import sunpy.io as io
import sunpy.wcs as wcs
from sunpy.visualization import toggle_pylab
from sunpy.sun import constants
from sunpy.sun import sun
from sunpy.time import parse_time, is_time
from sunpy.image.rescale import reshape_image_to_4d_superpixel
from sunpy.image.rescale import resample as sunpy_image_resample

__all__ = ['GenericMap']

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

"""
Questions
---------
* Should we use Helioviewer or VSO's data model? (e.g. map.meas, map.wavelength
or something else?)
* Should 'center' be renamed to 'offset' and crpix1 & 2 be used for 'center'?
"""

class GenericMap(astropy.nddata.NDData):
    """
    A Generic spatially-aware 2D data array

    Parameters
    ----------
    data : numpy.ndarray, list
        A 2d list or ndarray containing the map data
    header : dict
        A dictionary of the original image header tags

    Attributes
    ----------
    cmap : matplotlib.colors.Colormap
        A color map used for plotting with matplotlib.
    mpl_color_normalizer : matplotlib.colors.Normalize
        A matplotlib normalizer used to scale the image plot.

    Examples
    --------
    >>> import sunpy.map
    >>> aia = sunpy.map.Map(sunpy.AIA_171_IMAGE)
    >>> aia.T
    AIAMap([[ 0.3125,  1.    , -1.1875, ..., -0.625 ,  0.5625,  0.5   ],
    [-0.0625,  0.1875,  0.375 , ...,  0.0625,  0.0625, -0.125 ],
    [-0.125 , -0.8125, -0.5   , ..., -0.3125,  0.5625,  0.4375],
    ...,
    [ 0.625 ,  0.625 , -0.125 , ...,  0.125 , -0.0625,  0.6875],
    [-0.625 , -0.625 , -0.625 , ...,  0.125 , -0.0625,  0.6875],
    [ 0.    ,  0.    , -1.1875, ...,  0.125 ,  0.    ,  0.6875]])
    >>> aia.units['x']
    'arcsec'
    >>> aia.peek()

    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    | http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    | http://docs.scipy.org/doc/numpy/reference/ufuncs.html
    | http://www.scipy.org/Subclasses

    Notes
    -----

    This class makes some assumptions about the WCS information contained in
    the meta data. The first and most extensive assumption is that it is
    FITS-like WCS information as defined in the FITS WCS papers.

    Within this scope it also makes some other assumptions.

    * In the case of APIS convention headers where the CROTAi/j arguments are
      provided it assumes that these can be converted to the standard PCi_j
      notation using equations 32 in Thompson (2006).

    * If a CDi_j matrix is provided it is assumed that it can be converted to a
      PCi_j matrix and CDELT keywords as descirbed in Greisen & Calabretta (2002).

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

        astropy.nddata.NDData.__init__(self, data, meta=header, **kwargs)

        # Correct possibly missing meta keywords
        self._fix_date()
        self._fix_naxis()

        # Setup some attributes
        self._name = self.observatory + " " + str(self.measurement)
        self._nickname = self.detector

        # Visualization attributes
        self.cmap = cm.gray

        # Validate header
        # TODO: This should be a function of the header, not of the map
        self._validate()

        # Set mpl.colors.Normalize instance for plot scaling
        self.mpl_color_normalizer = self._get_mpl_normalizer()

    def __getitem__(self, key):
        """ This should allow indexing by physical coordinate """
        raise NotImplementedError(
    "The ability to index Map by physical coordinate is not yet implemented.")

    def __repr__(self):
        if not self.observatory:
            return self.data.__repr__()
        return (
"""SunPy %s
---------
Observatory:\t %s
Instrument:\t %s
Detector:\t %s
Measurement:\t %s
Obs Date:\t %s
dt:\t\t %f
Dimension:\t [%d, %d]
[dx, dy] =\t [%f, %f]

""" % (self.__class__.__name__,
       self.observatory, self.instrument, self.detector, self.measurement,
       self.date, self.exposure_time,
       self.data.shape[1], self.data.shape[0], self.scale['x'], self.scale['y'])
     + self.data.__repr__())


    #Some numpy extraction
    @property
    def shape(self):
        return self.data.shape

    @property
    def dtype(self):
        return self.data.dtype

    @property
    def size(self):
        return self.data.size

    @property
    def ndim(self):
        return self.data.ndim

    def std(self, *args, **kwargs):
        return self.data.std(*args, **kwargs)

    def mean(self, *args, **kwargs):
        return self.data.mean(*args, **kwargs)

    def min(self, *args, **kwargs):
        return self.data.min(*args, **kwargs)

    def max(self, *args, **kwargs):
        return self.data.max(*args, **kwargs)

# #### Keyword attribute and other attribute definitions #### #

    @property
    def name(self):
        """Human-readable description of map-type"""
        return self._name
    @name.setter
    def name(self, n):
        self._name = n

    @property
    def nickname(self):
        """An abbreviated human-readable description of the map-type; part of the Helioviewer data model"""
        return self._nickname
    @nickname.setter
    def nickname(self, n):
        self._nickname = n

    @property
    def date(self):
        """Image observation time"""
        return self.meta.get('date-obs', 'now')
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
            dsun = sun.sunearth_distance(self.date) * constants.au.si.value

        return dsun

    @property
    def exposure_time(self):
        """Exposure time of the image in seconds."""
        return self.meta.get('exptime', 0.0)

    @property
    def instrument(self):
        """Instrument name"""
        return self.meta.get('instrume', "")

    @property
    def measurement(self):
        """Measurement name, defaults to the wavelength of image"""
        return self.meta.get('wavelnth', "")

    @property
    def wavelength(self):
        """wavelength of the observation"""
        return self.meta.get('wavelnth', "")

    @property
    def observatory(self):
        """Observatory or Telescope name"""
        return self.meta.get('obsrvtry', self.meta.get('telescop', ""))

    @property
    def xrange(self):
        """Return the X range of the image in arcsec from edge to edge."""
        xmin = self.center['x'] - self.shape[1] / 2. * self.scale['x']
        xmax = self.center['x'] + self.shape[1] / 2. * self.scale['x']
        return [xmin, xmax]

    @property
    def yrange(self):
        """Return the Y range of the image in arcsec from edge to edge."""
        ymin = self.center['y'] - self.shape[0] / 2. * self.scale['y']
        ymax = self.center['y'] + self.shape[0] / 2. * self.scale['y']
        return [ymin, ymax]

    @property
    def center(self):
        """Returns the offset between the center of the Sun and the center of
        the map."""
        return {'x': wcs.get_center(self.shape[1], self.scale['x'],
                                    self.reference_pixel['x'],
                                    self.reference_coordinate['x']),
                'y': wcs.get_center(self.shape[0], self.scale['y'],
                                    self.reference_pixel['y'],
                                    self.reference_coordinate['y']),}

    @property
    def rsun_meters(self):
        """Radius of the sun in meters"""
        return self.meta.get('rsun_ref', constants.radius)

    @property
    def rsun_arcseconds(self):
        """Radius of the sun in arcseconds"""
        rsun_arcseconds = self.meta.get('rsun_obs',
                                        self.meta.get('solar_r',
                                                      self.meta.get('radius', None)))

        if rsun_arcseconds is None:
            warnings.warn_explicit("Missing metadata for solar radius: assuming photospheric limb as seen from Earth",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
            rsun_arcseconds = sun.solar_semidiameter_angular_size(self.date).value

        return rsun_arcseconds

    @property
    def coordinate_system(self):
        """Coordinate system used for x and y axes (ctype1/2)"""
        return {'x': self.meta.get('ctype1', 'HPLN-TAN'),
                'y': self.meta.get('ctype2', 'HPLT-TAN'),}

    @property
    def carrington_longitude(self):
        """Carrington longitude (crln_obs)"""
        carrington_longitude = self.meta.get('crln_obs', None)

        if carrington_longitude is None:
            warnings.warn_explicit("Missing metadata for Carrington longitude: assuming Earth-based observer",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
            carrington_longitude = (sun.heliographic_solar_center(self.date))[0]

        return carrington_longitude

    @property
    def heliographic_latitude(self):
        """Heliographic latitude in degrees"""
        heliographic_latitude = self.meta.get('hglt_obs',
                                              self.meta.get('crlt_obs',
                                                            self.meta.get('solar_b0', None)))

        if heliographic_latitude is None:
            warnings.warn_explicit("Missing metadata for heliographic latitude: assuming Earth-based observer",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
            heliographic_latitude = (sun.heliographic_solar_center(self.date))[1]

        return heliographic_latitude

    @property
    def heliographic_longitude(self):
        """Heliographic longitude in degrees"""
        return self.meta.get('hgln_obs', 0.)

    @property
    def reference_coordinate(self):
        """Reference point WCS axes in data units (crval1/2)"""
        return {'x': self.meta.get('crval1', 0.),
                'y': self.meta.get('crval2', 0.),}

    @property
    def reference_pixel(self):
        """Reference point axes in pixels (crpix1/2)"""
        return {'x': self.meta.get('crpix1', (self.meta.get('naxis1') + 1) / 2.),
                'y': self.meta.get('crpix2', (self.meta.get('naxis2') + 1) / 2.),}

    @property
    def scale(self):
        """Image scale along the x and y axes in units/pixel (cdelt1/2)"""
        #TODO: Fix this if only CDi_j matrix is provided
        return {'x': self.meta.get('cdelt1', 1.),
                'y': self.meta.get('cdelt2', 1.),}

    @property
    def units(self):
        """Image coordinate units along the x and y axes (cunit1/2)."""
        return {'x': self.meta.get('cunit1', 'arcsec'),
                'y': self.meta.get('cunit2', 'arcsec'),}

    @property
    def rotation_matrix(self):
        """Matrix describing the rotation required to align solar North with
        the top of the image."""
        if self.meta.get('PC1_1', None) is not None:
            return np.matrix([[self.meta['PC1_1'], self.meta['PC1_2']],
                              [self.meta['PC2_1'], self.meta['PC2_2']]])

        elif self.meta.get('CD1_1', None) is not None:
            div = 1. / (self.scale['x'] - self.scale['y'])

            deltm = np.matrix([[self.scale['y']/div, 0],
                               [0, self.scale['x']/ div]])

            cd = np.matrix([[self.meta['CD1_1'], self.meta['CD1_2']],
                            [self.meta['CD2_1'], self.meta['CD2_2']]])

            return deltm * cd
        else:
            return self._rotation_matrix_from_crota()

    def _rotation_matrix_from_crota(self):
        """
        This method converts the deprecated CROTA FITS kwargs to the new
        PC rotation matrix.

        This method can be overriden if an instruments header does not use this
        conversion.
        """
        lam = self.scale['y'] / self.scale['x']
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
            self.meta['naxis1'] = self.shape[1]
        if 'naxis2' not in self.meta:
            self.meta['naxis2'] = self.shape[0]
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

    def _validate(self):
        """Validates the meta-information associated with a Map.

        This function includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.map.sources package."""
#        if (self.dsun <= 0 or self.dsun >= 40 * constants.au):
#            raise InvalidHeaderInformation("Invalid value for DSUN")
        pass

# #### Data conversion routines #### #

    def data_to_pixel(self, value, dim):
        """Convert pixel-center data coordinates to pixel values"""
        #TODO: This function should be renamed. It is confusing as data
        # coordinates are in something like arcsec but this function just changes how you
        # count pixels
        if dim not in ['x', 'y']:
            raise ValueError("Invalid dimension. Must be one of 'x' or 'y'.")

        size = self.shape[dim == 'x']  # 1 if dim == 'x', 0 if dim == 'y'.

        return (value - self.center[dim]) / self.scale[dim] + ((size - 1) / 2.)

    def pixel_to_data(self, x=None, y=None):
        """Convert from pixel coordinates to data coordinates (e.g. arcsec)"""
        width = self.shape[1]
        height = self.shape[0]

        if (x is not None) & (x > width-1):
            raise ValueError("X pixel value larger than image width (%s)." % width)
        if (x is not None) & (y > height-1):
            raise ValueError("Y pixel value larger than image height (%s)." % height)
        if (x is not None) & (x < 0):
            raise ValueError("X pixel value cannot be less than 0.")
        if (x is not None) & (y < 0):
            raise ValueError("Y pixel value cannot be less than 0.")

        scale = np.array([self.scale['x'], self.scale['y']])
        crpix = np.array([self.reference_pixel['x'], self.reference_pixel['y']])
        crval = np.array([self.reference_coordinate['x'], self.reference_coordinate['y']])
        coordinate_system = [self.coordinate_system['x'], self.coordinate_system['y']]
        x,y = wcs.convert_pixel_to_data(self.shape, scale, crpix, crval, x = x, y = y)

        return x, y

# #### I/O routines #### #

    def save(self, filepath, filetype='auto', **kwargs):
        """Saves the SunPy Map object to a file.

        Currently SunPy can only save files in the FITS format. In the future
        support will be added for saving to other formats.

        Parameters
        ----------
        filepath : string
            Location to save file to.

        filetype : string
            'auto' or any supported file extension
        """
        io.write_file(filepath, self.data, self.meta, filetype=filetype,
                      **kwargs)

# #### Image processing routines #### #

    def resample(self, dimensions, method='linear'):
        """Returns a new Map that has been resampled up or down

        Arbitrary resampling of the Map to new dimension sizes.

        Uses the same parameters and creates the same co-ordinate lookup points
        as IDL''s congrid routine, which apparently originally came from a
        VAX/VMS routine of the same name.

        Parameters
        ----------
        dimensions : tuple
            Dimensions that new Map should have.
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
        out : Map
            A new Map which has been resampled to the desired dimensions.

        References
        ----------
        | http://www.scipy.org/Cookbook/Rebinning (Original source, 2011/11/19)
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

        # Note that 'x' and 'y' correspond to 1 and 0 in self.shape,
        # respectively
        scale_factor_x = (float(self.shape[1]) / dimensions[0])
        scale_factor_y = (float(self.shape[0]) / dimensions[1])

        new_map = deepcopy(self)
        # Update image scale and number of pixels
        new_meta = self.meta.copy()

        # Update metadata
        new_meta['cdelt1'] *= scale_factor_x
        new_meta['cdelt2'] *= scale_factor_y
        new_meta['crpix1'] = (dimensions[0] + 1) / 2.
        new_meta['crpix2'] = (dimensions[1] + 1) / 2.
        new_meta['crval1'] = self.center['x']
        new_meta['crval2'] = self.center['y']

        # Create new map instance
        new_map.data = new_data
        new_map.meta = new_meta
        return new_map

    def rotate(self, angle=None, rmatrix=None, order=3, scale=1.0,
               image_center=(0,0), recenter=False, missing=0.0, use_scipy=False):
        """
        Returns a new rotated and rescaled map.  Specify either a rotation
        angle or a rotation matrix, but not both.  If neither an angle or a
        rotation matrix are specified, the map will be rotated by the rotation
        angle in the metadata.

        Also updates the rotation_matrix attribute and any appropriate header
        data so that they correctly describe the new map.

        Parameters
        ----------
        angle : float
            The angle (degrees) to rotate counterclockwise.
        rmatrix : 2x2
            Linear transformation rotation matrix.
        order : int 0-5
            Interpolation order to be used. When using scikit-image this parameter
            is passed into :func:`skimage.transform.warp`.
            When using scipy it is passed into
            :func:`scipy.ndimage.interpolation.affine_transform` where it controls
            the order of the spline.
            Higher accuracy may be obtained at the cost of performance by using
            higher values.
        scale : float
            A scale factor for the image, default is no scaling
        image_center : tuple
            The axis of rotation in data coordinates
            Default: the origin in the data coordinate system
        recenter : bool
            If True, position the axis of rotation at the center of the new map
            Default: False
        missing : float
            The numerical value to fill any missing points after rotation.
            Default: 0.0
        use_scipy : bool
            If True, forces the rotation to use
            :func:`scipy.ndimage.interpolation.affine_transform`, otherwise it
            uses the :class:`skimage.transform.AffineTransform` class and
            :func:`skimage.transform.warp`.
            The function will also automatically fall back to
            :func:`scipy.ndimage.interpolation.affine_transform` if scikit-image
            can't be imported.
            Default: False

        Returns
        -------
        out : Map
            A new Map instance containing the rotated and rescaled data of the
            original map.

        See Also
        --------
        sunpy.image.transform.affine_transform : The routine this method calls for the rotation.

        Notes
        -----
        This function will remove old CROTA keywords from the header.
        This function will also convert a CDi_j matrix to a PCi_j matrix.

        The scikit-image and scipy affine_transform routines do not use the same algorithm,
        see :func:`sunpy.image.transform.affine_transform` for details.

        This function is not numerically equalivalent to IDL's rot() see the
        :func:`sunpy.image.transform.affine_transform` documentation for a
        detailed description of the differences.
        """
        if angle is not None and rmatrix is not None:
            raise ValueError("You cannot specify both an angle and a matrix")
        elif angle is None and rmatrix is None:
            rmatrix = self.rotation_matrix

        # Interpolation parameter sanity
        if order not in range(6):
            raise ValueError("Order must be between 0 and 5")

        # Copy Map
        new_map = deepcopy(self)

        if angle is not None:
            #Calulate the parameters for the affine_transform
            c = np.cos(np.deg2rad(angle))
            s = np.sin(np.deg2rad(angle))
            rmatrix = np.matrix([[c, -s], [s, c]])

        # map_center is swapped compared to the x-y convention
        array_center = (np.array(self.data.shape)-1)/2.0

        # rotation_center is swapped compared to the x-y convention
        if recenter:
            # Convert the axis of rotation from data coordinates to pixel coordinates
            x = self.data_to_pixel(image_center[0], 'x')
            y = self.data_to_pixel(image_center[1], 'y')
            rotation_center = (y, x)
        else:
            rotation_center = array_center

        #Return a new map
        #Copy Header
        new_map = deepcopy(self)

        new_map.data = affine_transform(new_map.data.T,
                                        np.asarray(rmatrix),
                                        order=order, scale=scale,
                                        image_center=rotation_center,
                                        recenter=recenter, missing=missing,
                                        use_scipy=use_scipy).T


        # Calculate new reference pixel and coordinate at the center of the
        # image.
        if recenter:
            new_center = image_center
        else:
            # Retrieve old coordinates for the center of the array
            old_center = np.asarray(self.pixel_to_data(array_center[1], array_center[0]))

            # Calculate new coordinates for the center of the array
            new_center = image_center - np.dot(rmatrix, image_center - old_center)
            new_center = np.asarray(new_center)[0]

        # Define a new reference pixel in the rotated space
        new_map.meta['crval1'] = new_center[0]
        new_map.meta['crval2'] = new_center[1]
        new_map.meta['crpix1'] = array_center[1] + 1 # FITS counts pixels from 1
        new_map.meta['crpix2'] = array_center[0] + 1 # FITS counts pixels from 1

        # Calculate the new rotation matrix to store in the header by
        # "subtracting" the rotation matrix used in the rotate from the old one
        # That being calculate the dot product of the old header data with the
        # inverse of the rotation matrix.
        pc_C = np.dot(self.rotation_matrix, rmatrix.I)
        new_map.meta['PC1_1'] = pc_C[0,0]
        new_map.meta['PC1_2'] = pc_C[0,1]
        new_map.meta['PC2_1'] = pc_C[1,0]
        new_map.meta['PC2_2'] = pc_C[1,1]

        # Update pixel size if image has been scaled.
        if scale != 1.0:
            new_map.meta['cdelt1'] = self.scale['x'] / scale
            new_map.meta['cdelt2'] = self.scale['y'] / scale

        # Remove old CROTA kwargs because we have saved a new PCi_j matrix.
        new_map.meta.pop('CROTA1', None)
        new_map.meta.pop('CROTA2', None)
        # Remove CDi_j header
        new_map.meta.pop('CD1_1', None)
        new_map.meta.pop('CD1_2', None)
        new_map.meta.pop('CD2_1', None)
        new_map.meta.pop('CD2_2', None)

        return new_map

    def submap(self, range_a, range_b, units="data"):
        """Returns a submap of the map with the specified range

        Parameters
        ----------
        range_a : list
            The range of the Map to select across either the x axis.
        range_b : list
            The range of the Map to select across either the y axis.
        units : {'data' | 'pixels'}, optional
            The units for the supplied ranges.

        Returns
        -------
        out : Map
            A new map instance is returned representing to specified sub-region

        Examples
        --------
        >>> aia.submap([-5,5],[-5,5])
        AIAMap([[ 341.3125,  266.5   ,  329.375 ,  330.5625,  298.875 ],
        [ 347.1875,  273.4375,  247.4375,  303.5   ,  305.3125],
        [ 322.8125,  302.3125,  298.125 ,  299.    ,  261.5   ],
        [ 334.875 ,  289.75  ,  269.25  ,  256.375 ,  242.3125],
        [ 273.125 ,  241.75  ,  248.8125,  263.0625,  249.0625]])

        >>> aia.submap([0,5],[0,5], units='pixels')
        AIAMap([[ 0.3125, -0.0625, -0.125 ,  0.    , -0.375 ],
        [ 1.    ,  0.1875, -0.8125,  0.125 ,  0.3125],
        [-1.1875,  0.375 , -0.5   ,  0.25  , -0.4375],
        [-0.6875, -0.3125,  0.8125,  0.0625,  0.1875],
        [-0.875 ,  0.25  ,  0.1875,  0.    , -0.6875]])
        """
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

            #x_pixels = [self.data_to_pixel(elem, 'x') for elem in range_a]
            x_pixels = [np.ceil(self.data_to_pixel(range_a[0], 'x')),
                        np.floor(self.data_to_pixel(range_a[1], 'x')) + 1]
            #y_pixels = [self.data_to_pixel(elem, 'y') for elem in range_b]
            y_pixels = [np.ceil(self.data_to_pixel(range_b[0], 'y')),
                        np.floor(self.data_to_pixel(range_b[1], 'y')) + 1]
        elif units is "pixels":
            # Check edges
            if range_a[0] is None:
                range_a[0] = 0
            if range_a[1] is None:
                range_a[1] = self.shape[1]
            if range_b[0] is None:
                range_b[0] = 0
            if range_b[1] is None:
                range_b[1] = self.shape[0]

            x_pixels = range_a
            y_pixels = range_b
        else:
            raise ValueError(
                "Invalid unit. Must be one of 'data' or 'pixels'")


        # Get ndarray representation of submap
        xslice = slice(x_pixels[0], x_pixels[1])
        yslice = slice(y_pixels[0], y_pixels[1])
        new_data = self.data[yslice, xslice].copy()

        # Make a copy of the header with updated centering information
        new_map = deepcopy(self)
        new_map.meta['crpix1'] = self.reference_pixel['x'] - x_pixels[0]
        new_map.meta['crpix2'] = self.reference_pixel['y'] - y_pixels[0]
        new_map.meta['naxis1'] = new_data.shape[1]
        new_map.meta['naxis2'] = new_data.shape[0]

        # Create new map instance
        new_map.data = new_data
        return new_map

    def superpixel(self, dimensions, method='sum'):
        """Returns a new map consisting of superpixels formed from the
        original data.  Useful for increasing signal to noise ratio in images.

        Parameters
        ----------
        dimensions : tuple
            One superpixel in the new map is equal to (dimension[0],
            dimension[1]) pixels of the original map
            Note: the first argument corresponds to the 'x' axis and the second
            argument corresponds to the 'y' axis.
        method : {'sum' | 'average'}
            What each superpixel represents compared to the original data
                * sum - add up the original data
                * average - average the sum over the number of original pixels

        Returns
        -------
        out : Map
            A new Map which has superpixels of the required size.

        References
        ----------
        | http://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html
        """

        # Note: because the underlying ndarray is transposed in sense when
        #   compared to the Map, the ndarray is transposed, resampled, then
        #   transposed back
        # Note: "center" defaults to True in this function because data
        #   coordinates in a Map are at pixel centers

        # Make a copy of the original data and perform reshaping
        reshaped = reshape_image_to_4d_superpixel(self.data.copy().T,
                                                  dimensions)
        if method == 'sum':
            new_data = reshaped.sum(axis=3).sum(axis=1)
        elif method == 'average':
            new_data = ((reshaped.sum(axis=3).sum(axis=1)) /
                    np.float32(dimensions[0] * dimensions[1]))
        new_data = new_data.T


        # Update image scale and number of pixels
        new_map = deepcopy(self)
        new_meta = new_map.meta

        # Note that 'x' and 'y' correspond to 1 and 0 in self.shape,
        # respectively
        new_nx = self.shape[1] / dimensions[0]
        new_ny = self.shape[0] / dimensions[1]

        # Update metadata
        new_meta['cdelt1'] = dimensions[0] * self.scale['x']
        new_meta['cdelt2'] = dimensions[1] * self.scale['y']
        new_meta['crpix1'] = (new_nx + 1) / 2.
        new_meta['crpix2'] = (new_ny + 1) / 2.
        new_meta['crval1'] = self.center['x']
        new_meta['crval2'] = self.center['y']

        # Create new map instance
        new_map.data = new_data
        return new_map

# #### Visualization #### #

    def draw_grid(self, axes=None, grid_spacing=15, **kwargs):
        """Draws a grid over the surface of the Sun

        Parameters
        ----------
        axes: matplotlib.axes object or None
        Axes to plot limb on or None to use current axes.

        grid_spacing: float
            Spacing (in degrees) for longitude and latitude grid.

        Returns
        -------
        matplotlib.axes object

        Notes
        -----
        keyword arguments are passed onto matplotlib.pyplot.plot
        """

        if not axes:
            axes = plt.gca()

        x, y = self.pixel_to_data()
        dsun = self.dsun

        b0 = self.heliographic_latitude
        l0 = self.heliographic_longitude
        units = [self.units['x'], self.units['y']]

        #Prep the plot kwargs
        plot_kw = {'color':'white',
                   'linestyle':'dotted',
                   'zorder':100}
        plot_kw.update(kwargs)

        hg_longitude_deg = np.linspace(-180, 180, num=361) + self.heliographic_longitude
        hg_latitude_deg = np.arange(-90, 90, grid_spacing)

        # draw the latitude lines
        for lat in hg_latitude_deg:
            x, y = wcs.convert_hg_hpc(hg_longitude_deg, lat * np.ones(361),
                                      b0_deg=b0, l0_deg=l0, dsun_meters=dsun,
                                      angle_units=units[0], occultation=True)
            valid = np.logical_and(np.isfinite(x), np.isfinite(y))
            x = x[valid]
            y = y[valid]
            axes.plot(x, y, **plot_kw)

        hg_longitude_deg = np.arange(-180, 180, grid_spacing) + self.heliographic_longitude
        hg_latitude_deg = np.linspace(-90, 90, num=181)

        # draw the longitude lines
        for lon in hg_longitude_deg:
            x, y = wcs.convert_hg_hpc(lon * np.ones(181), hg_latitude_deg,
                                      b0_deg=b0, l0_deg=l0, dsun_meters=dsun,
                                      angle_units=units[0], occultation=True)
            valid = np.logical_and(np.isfinite(x), np.isfinite(y))
            x = x[valid]
            y = y[valid]
            axes.plot(x, y, **plot_kw)

        axes.set_ylim(self.yrange)
        axes.set_xlim(self.xrange)

        return axes

    def draw_limb(self, axes=None, **kwargs):
        """Draws a circle representing the solar limb

            Parameters
            ----------
            axes: matplotlib.axes object or None
                Axes to plot limb on or None to use current axes.

            Returns
            -------
            matplotlib.axes object

            Notes
            -----
            keyword arguments are passed onto the Circle Patch, see:
            http://matplotlib.org/api/artist_api.html#matplotlib.patches.Patch
            http://matplotlib.org/api/artist_api.html#matplotlib.patches.Circle
        """

        if not axes:
            axes = plt.gca()

        c_kw = {'radius':self.rsun_arcseconds,
                'fill':False,
                'color':'white',
                'zorder':100
                }
        c_kw.update(kwargs)

        circ = patches.Circle([0, 0], **c_kw)
        axes.add_artist(circ)

        return axes

    @toggle_pylab
    def peek(self, draw_limb=False, draw_grid=False, gamma=None,
                   colorbar=True, basic_plot=False, **matplot_args):
        """Displays the map in a new figure

        Parameters
        ----------
        draw_limb : bool
            Whether the solar limb should be plotted.
        draw_grid : bool or number
            Whether solar meridians and parallels are plotted. If float then sets
            degree difference between parallels and meridians.
        gamma : float
            Gamma value to use for the color map
        colorbar : bool
            Whether to display a colorbar next to the plot
        basic_plot : bool
            If true, the data is plotted by itself at it's natural scale; no
            title, labels, or axes are shown.
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting the image.
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
            axes = figure.gca()

        im = self.plot(axes=axes,**matplot_args)

        if colorbar and not basic_plot:
            figure.colorbar(im)

        if draw_limb:
            self.draw_limb(axes=axes)

        if isinstance(draw_grid, bool):
            if draw_grid:
                self.draw_grid(axes=axes)
        elif isinstance(draw_grid, (int, long, float)):
            self.draw_grid(axes=axes, grid_spacing=draw_grid)
        else:
            raise TypeError("draw_grid should be bool, int, long or float")

        figure.show()

        return figure

    @toggle_pylab
    def plot(self, gamma=None, annotate=True, axes=None, **imshow_args):
        """ Plots the map object using matplotlib, in a method equivalent
        to plt.imshow() using nearest neighbour interpolation.

        Parameters
        ----------
        gamma : float
            Gamma value to use for the color map

        annotate : bool
            If true, the data is plotted at it's natural scale; with
            title and axis labels.

        axes: matplotlib.axes object or None
            If provided the image will be plotted on the given axes. Else the
            current matplotlib axes will be used.

        **imshow_args : dict
            Any additional imshow arguments that should be used
            when plotting the image.

        Examples
        --------
        #Simple Plot with color bar
        plt.figure()
        aiamap.plot()
        plt.colorbar()

        #Add a limb line and grid
        aia.plot()
        aia.draw_limb()
        aia.draw_grid()
        """
        # Check that the image is properly aligned
        if not np.array_equal(self.rotation_matrix, np.matrix(np.identity(2))):
            warnings.warn("This map is not aligned. Plot axes may be incorrect",
                          Warning)
        #Get current axes
        if not axes:
            axes = plt.gca()

        # Normal plot
        if annotate:
            axes.set_title("%s %s" % (self.name, parse_time(self.date).strftime(TIME_FORMAT)))

            # x-axis label
            if self.coordinate_system['x'] == 'HG':
                xlabel = 'Longitude [%s]' % self.units['x']
            else:
                xlabel = 'X-position [%s]' % self.units['x']

            # y-axis label
            if self.coordinate_system['y'] == 'HG':
                ylabel = 'Latitude [%s]' % self.units['y']
            else:
                ylabel = 'Y-position [%s]' % self.units['y']

            axes.set_xlabel(xlabel)
            axes.set_ylabel(ylabel)

        # Determine extent
        extent = self.xrange + self.yrange

        cmap = deepcopy(self.cmap)
        if gamma is not None:
            cmap.set_gamma(gamma)

        # make imshow kwargs a dict
        kwargs = {'origin':'lower',
                  'cmap':cmap,
                  'norm':self.mpl_color_normalizer,
                  'extent':extent,
                  'interpolation':'nearest'}
        kwargs.update(imshow_args)

        ret = axes.imshow(self.data, **kwargs)

        #Set current image (makes colorbar work)
        plt.sci(ret)
        return ret

    def _get_mpl_normalizer(self):
        """
        Returns a default mpl.colors.Normalize instance for plot scaling.

        Not yet implemented.
        """
        return None


class InvalidHeaderInformation(ValueError):
    """Exception to raise when an invalid header tag value is encountered for a
    FITS/JPEG 2000 file."""
    pass
