"""
Map is a generic Map class from which all other Map classes inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404,W0613
__authors__ = ["Russell Hewett, Stuart Mumford, Keith Hughitt, Steven Christe"]
__email__ = "stuart@mumford.me.uk"

import os
from copy import deepcopy, copy
import warnings

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.interpolation
from matplotlib import patches
from matplotlib import colors
from matplotlib import cm

import astropy.nddata

try:
    import sunpy.image.Crotate as Crotate
except ImportError:
    pass

import sunpy.io as io
import sunpy.wcs as wcs
from sunpy.util import to_signed, Deprecated
from sunpy.visualization import toggle_pylab
# from sunpy.io import read_file, read_file_header
from sunpy.sun import constants
from sunpy.time import parse_time, is_time
from sunpy.image.rescale import reshape_image_to_4d_superpixel
from sunpy.image.rescale import resample as sunpy_image_resample

#from sunpy.util.cond_dispatch import ConditionalDispatch
#from sunpy.util.create import Parent

__all__ = ['GenericMap']

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

    Examples
    --------
    >>> aia = sunpy.make_map(sunpy.AIA_171_IMAGE)
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

        self.norm = self._get_norm()

    def __getitem__(self, key):
        """ This should allow indexing by physical coordinate """
        raise NotImplementedError(
    "The ability to index Map by physical coordinate is not yet implemented.")

    def __repr__(self):
        if not hasattr(self, 'observatory'):
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
        return self.meta.get('date-obs', None)
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
        return self.meta.get('dsun_obs', constants.au)

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
        return self.meta.get('rsun_obs', self.meta.get('solar_r',
                                         self.meta.get('radius', constants.average_angular_size.to('arcsec').value)))

    @property
    def coordinate_system(self):
        """Coordinate system used for x and y axes (ctype1/2)"""
        return {'x': self.meta.get('ctype1', 'HPLN-TAN'),
                'y': self.meta.get('ctype2', 'HPLT-TAN'),}

    @property
    def carrington_longitude(self):
        """Carrington longitude (crln_obs)"""
        return self.meta.get('crln_obs', 0.)

    @property
    def heliographic_latitude(self):
        """Heliographic latitude in degrees"""
        return self.meta.get('hglt_obs', self.meta.get('crlt_obs',
                                         self.meta.get('solar_b0', 0.)))

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
        return {'x': self.meta.get('cdelt1', 1.),
                'y': self.meta.get('cdelt2', 1.),}

    @property
    def units(self):
        """Image coordinate units along the x and y axes (cunit1/2)."""
        return {'x': self.meta.get('cunit1', 'arcsec'),
                'y': self.meta.get('cunit2', 'arcsec'),}

    #TODO: This needs to be WCS compliant!
    @property
    def rotation_angle(self):
        """The Rotation angle of each axis"""
        return {'x': self.meta.get('crota1', 0.),
                'y': self.meta.get('crota2', 0.),}

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
    
    def rotate(self, angle=None, rmatrix=None, scale=1.0, rotation_center=None, recenter=True,
               missing=0.0, interpolation='bicubic', interp_param=-0.5):
        """Returns a new rotated, rescaled and shifted map.

        Parameters
        ----------
        angle: float
           The angle to rotate the image by (radians). Specify angle or matrix.
        rmatrix: NxN
            Linear transformation rotation matrix. Specify angle or matrix.
        scale: float
           A scale factor for the image, default is no scaling
        rotation_center: tuple
           The point in the image to rotate around (Axis of rotation).
           Default: center of the array
        recenter: bool, or array-like
           Move the centroid (axis of rotation) to the center of the array
           or recenter coords.
           Default: True, recenter to the center of the array.
        missing: float
           The numerical value to fill any missing points after rotation.
           Default: 0.0
        interpolation: {'nearest' | 'bilinear' | 'spline' | 'bicubic'}
            Interpolation method to use in the transform.
            Spline uses the
            scipy.ndimage.interpolation.affline_transform routine.
            nearest, bilinear and bicubic all replicate the IDL rot() function.
            Default: 'bicubic'
        interp_par: Int or Float
            Optional parameter for controlling the interpolation.
            Spline interpolation requires an integer value between 1 and 5 for
            the degree of the spline fit.
            Default: 3
            BiCubic interpolation requires a flaot value between -1 and 0.
            Default: 0.5
            Other interpolation options ingore the argument.

        Returns
        -------
        New rotated, rescaled, translated map

        Notes
        -----
        Apart from interpolation='spline' all other options use a compiled
        C-API extension. If for some reason this is not compiled correctly this
        routine will fall back upon the scipy implementation of order = 3.
        For more infomation see:
        http://sunpy.readthedocs.org/en/latest/guide/troubleshooting.html#crotate-warning
        """
        assert angle is None or rmatrix is None
        #Interpolation parameter Sanity
        assert interpolation in ['nearest','spline','bilinear','bicubic']
        #Set defaults based on interpolation
        if interp_param is None:
            if interpolation is 'spline':
                interp_param = 3
            elif interpolation is 'bicubic':
                interp_param = 0.5
            else:
                interp_param = 0 #Default value for nearest or bilinear

        #Make sure recenter is a vector with shape (2,1)
        if not isinstance(recenter, bool):
            recenter = np.array(recenter).reshape(2,1)

        #Define Size and center of array
        center = (np.array(self.data.shape)-1)/2.0

        #If rotation_center is not set (None or False),
        #set rotation_center to the center of the image.
        if rotation_center is None:
            rotation_center = center
        else:
            #Else check rotation_center is a vector with shape (2,1)
            rotation_center = np.array(rotation_center).reshape(2,1)

        #recenter to the rotation_center if recenter is True
        if isinstance(recenter, bool):
            #if rentre is False then this will be (0,0)
            shift = np.array(rotation_center) - np.array(center)
        else:
            #recenter to recenter vector otherwise
            shift = np.array(recenter) - np.array(center)

        image = self.data.copy()

        if not angle is None:
            #Calulate the parameters for the affline_transform
            c = np.cos(angle)
            s = np.sin(angle)
            mati = np.array([[c, s],[-s, c]]) / scale   # res->orig
        if not rmatrix is None:
            mati = rmatrix / scale   # res->orig
        center = np.array([center]).transpose()  # the center of rotn
        shift = np.array([shift]).transpose()    # the shift
        kpos = center - np.dot(mati, (center + shift))
        # kpos and mati are the two transform constants, kpos is a 2x2 array
        rsmat, offs =  mati, np.squeeze((kpos[0,0], kpos[1,0]))

        if interpolation == 'spline':
            # This is the scipy call
            data = scipy.ndimage.interpolation.affine_transform(image, rsmat,
                           offset=offs, order=interp_param, mode='constant',
                           cval=missing)
        else:
            #Use C extension Package
            if not 'Crotate' in globals():
                warnings.warn("""The C extension sunpy.image.Crotate is not
installed, falling back to the interpolation='spline' of order=3""" ,Warning)
                data = scipy.ndimage.interpolation.affine_transform(image, rsmat,
                           offset=offs, order=3, mode='constant',
                           cval=missing)
            #Set up call parameters depending on interp type.
            if interpolation == 'nearest':
                interp_type = Crotate.NEAREST
            elif interpolation == 'bilinear':
                interp_type = Crotate.BILINEAR
            elif interpolation == 'bicubic':
                interp_type = Crotate.BICUBIC
            #Make call to extension
            data = Crotate.affine_transform(image,
                                      rsmat, offset=offs,
                                      kernel=interp_type, cubic=interp_param,
                                      mode='constant', cval=missing)

        #Return a new map
        #Copy Header
        new_map = deepcopy(self)

        # Create new map instance
        new_map.data = data
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

    def draw_grid(self, axes=None, grid_spacing=20, **kwargs):
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
        rsun = self.rsun_meters
        dsun = self.dsun

        b0 = self.heliographic_latitude
        l0 = self.heliographic_longitude
        units = [self.units['x'], self.units['y']]

        #Prep the plot kwargs
        plot_kw = {'color':'white',
                   'linestyle':'dotted',
                   'zorder':100}
        plot_kw.update(kwargs)

        #TODO: This function could be optimized. Does not need to convert the entire image
        # coordinates
        #lon_self, lat_self = wcs.convert_hpc_hg(rsun, dsun, angle_units = units[0], b0, l0, x, y)
        lon_self, lat_self = wcs.convert_hpc_hg(x, y, b0_deg=b0, l0_deg=l0, dsun_meters=dsun, angle_units='arcsec')
        # define the number of points for each latitude or longitude line
        num_points = 20

        #TODO: The following code is ugly. Fix it.
        lon_range = [lon_self.min(), lon_self.max()]
        lat_range = [lat_self.min(), lat_self.max()]
        if np.isfinite(lon_range[0]) == False:
            lon_range[0] = -90 + self.heliographic_longitude
        if np.isfinite(lon_range[1]) == False:
            lon_range[1] = 90 + self.heliographic_longitude
        if np.isfinite(lat_range[0]) == False:
            lat_range[0] = -90 + self.heliographic_latitude
        if np.isfinite(lat_range[1]) == False:
            lat_range[1] = 90 + self.heliographic_latitude

        hg_longitude_deg = np.linspace(lon_range[0], lon_range[1], num=num_points)
        hg_latitude_deg = np.arange(lat_range[0], lat_range[1]+grid_spacing, grid_spacing)

        # draw the latitude lines
        for lat in hg_latitude_deg:
            hg_latitude_deg_mesh, hg_longitude_deg_mesh = np.meshgrid(
                lat * np.ones(num_points), hg_longitude_deg)
            x, y = wcs.convert_hg_hpc(hg_longitude_deg_mesh, hg_latitude_deg_mesh, b0_deg=b0, l0_deg=l0,
                    dsun_meters=dsun, angle_units=units[0], occultation=False)

            axes.plot(x, y, **plot_kw)

        hg_longitude_deg = np.arange(lon_range[0], lon_range[1]+grid_spacing, grid_spacing)
        hg_latitude_deg = np.linspace(lat_range[0], lat_range[1], num=num_points)

        # draw the longitude lines
        for lon in hg_longitude_deg:
            hg_longitude_deg_mesh, hg_latitude_deg_mesh = np.meshgrid(
                lon * np.ones(num_points), hg_latitude_deg)
            x, y = wcs.convert_hg_hpc(hg_longitude_deg_mesh, hg_latitude_deg_mesh, b0_deg=b0, l0_deg=l0,
                    dsun_meters=dsun, angle_units=units[0], occultation=False)
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
    def peek(self, draw_limb=True, draw_grid=False, gamma=None,
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

        #Get current axes
        if not axes:
            axes = plt.gca()

        # Normal plot
        if annotate:
            axes.set_title("%s %s" % (self.name, parse_time(self.date).strftime("%Y-%m-%d %H:%M:%S.%f")))

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

            #make imshow kwargs a dict

        kwargs = {'origin':'lower',
                  'cmap':cmap,
                  'norm':self.norm,
                  'extent':extent,
                  'interpolation':'nearest'}
        kwargs.update(imshow_args)

        ret = axes.imshow(self.data, **kwargs)

        #Set current image (makes colorbar work)
        plt.sci(ret)
        return ret

    def _get_norm(self):
        """Default normalization method. Not yet implemented."""
        return None


class InvalidHeaderInformation(ValueError):
    """Exception to raise when an invalid header tag value is encountered for a
    FITS/JPEG 2000 file."""
    pass
