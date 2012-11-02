"""
Map is a generic Map class from which all other Map classes inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404,W0613
__authors__ = ["Keith Hughitt, Steven Christe"]
__email__ = "keith.hughitt@nasa.gov"

import os
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.interpolation
from matplotlib import patches
from matplotlib import colors
from matplotlib import cm
from copy import copy
try:
    import sunpy.image.Crotate as Crotate
except ImportError:
    print("C extension sunpy.image.Crotate cannot be found")
from sunpy.wcs import wcs as wcs
from sunpy.util.util import toggle_pylab
from sunpy.io import read_file, read_file_header
from sunpy.sun import constants
from sunpy.time import parse_time
from sunpy.util.util import to_signed
from sunpy.image.rescale import resample, reshape_image_to_4d_superpixel

from sunpy.util.cond_dispatch import ConditionalDispatch
from sunpy.util.create import Parent

"""
TODO
----
* Automatically include Map docstring when displaying help for subclasses?

Questions
---------
* Should we use Helioviewer or VSO's data model? (e.g. map.meas, map.wavelength
or something else?)
* Should 'center' be renamed to 'offset' and crpix1 & 2 be used for 'center'?
"""

class Map(np.ndarray, Parent):
    """
    Map(data, header)

    A spatially-aware data array based on the SolarSoft Map object

    Parameters
    ----------
    data : numpy.ndarray, list
        A 2d list or ndarray containing the map data
    header : dict
        A dictionary of the original image header tags

    Attributes
    ----------
    original_header : dict
        Dictionary representation of the original FITS header
    carrington_longitude : str
        Carrington longitude (crln_obs)
    center : dict
        X and Y coordinate of the center of the map in units.
        Usually represents the offset between the center of the Sun and the
        center of the map.
    cmap : matplotlib.colors.Colormap
        A Matplotlib colormap to be applied to the data
    coordinate_system : dict
        Coordinate system used for x and y axes (ctype1/2)
    date : datetime
        Image observation time
    detector : str
        Detector name
    dsun : float
        The observer distance from the Sun.
    exptime : float
        Exposure time of the image in seconds.
    heliographic_latitude : float
        Heliographic latitude in degrees
    heliographic_longitude : float
        Heliographic longitude in degrees
    instrument : str
        Instrument name
    measurement : str, int
        Measurement name. In some instances this is the wavelength of image.
    name: str
        Human-readable description of map-type
    nickname : str
        An abbreviated human-readable description of the map-type; part of
        the Helioviewer data model
    observatory : str
        Observatory name
    reference_coordinate : float
        Reference point WCS axes in data units (crval1/2) 
    reference_pixel : float
        Reference point axes in pixels (crpix1/2)
    rsun_arcseconds : float
        Radius of the sun in arcseconds
    rsun_meters : float
        Radius of the sun in meters
    scale : dict
        Image scale along the x and y axes in units/pixel (cdelt1/2).
    units : dict
        Image coordinate units along the x and y axes (cunit1/2).

    Methods
    -------
    std()
        Return the standard deviation of the map data
    mean()
        Return the mean of the map data
    min()
        Return the minimum value of the map data
    max()
        Return the maximum value of the map data
    resample(dimension, method)
        Returns a new map that has been resampled up or down
    superpixel(dimension, method)
        Returns a new map consisting of superpixels formed from the
        original data.
    save()
        Save the map to a fits file.
    submap(range_a, range_b, units)
        Returns a submap of the map with the specified range
    plot()
        Return a matplotlib plot figure object
    show()
        Display a matplotlib plot to the screen 
    get_header()
        Returns the original header from when the map was first created.

    Examples
    --------
    >>> aia = sunpy.Map(sunpy.AIA_171_IMAGE)
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
    >>> aia.show()
    >>> import matplotlib.cm as cm
    >>> import matplotlib.colors as colors
    >>> aia.show(cmap=cm.hot, norm=colors.Normalize(1, 2048))

    See Also
    --------
    numpy.ndarray Parent class for the Map object

    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    | http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    | http://docs.scipy.org/doc/numpy/reference/ufuncs.html
    | http://www.scipy.org/Subclasses

    """
    _create = ConditionalDispatch.from_existing(Parent._create)
    create = classmethod(_create.wrapper())

    def __new__(cls, data, header):
        """Creates a new Map instance"""
        if isinstance(data, np.ndarray):
            obj = data.view(cls)
        elif isinstance(data, list):
            obj = np.asarray(data).view(cls)
        else:
            raise TypeError('Invalid input')

        return obj

    def __init__(self, data, header):
        self._original_header = header
        
        # Set naxis1 and naxis2 if not specified
        if header.get('naxis1') is None:
            header['naxis1'] = self.shape[1]
        if header.get('naxis2') is None:
            header['naxis2'] = self.shape[0]

        # Parse header and set map attributes
        for attr, value in list(self.get_properties(header).items()):
            setattr(self, attr, value)
        
        # Validate properties
        self._validate()
        
    @classmethod
    def get_properties(cls, header):
        """Parses a map header and determines default properties."""
        return {
            "cmap": cm.gray,  # @UndefinedVariable
            "date": parse_time(header.get('date-obs', None)),
            "detector": header.get('detector', ''),
            "dsun": header.get('dsun_obs', constants.au),
            "exposure_time": header.get('exptime', 0.),
            "instrument": header.get('instrume', ''),
            "measurement": header.get('wavelnth', ''),
            "observatory": header.get('telescop', ''),
            "name": header.get('telescop', '') + " " + 
                    str(header.get('wavelnth', '')),
            "nickname": header.get('detector', ''),
            "rsun_meters": header.get('rsun_ref', constants.radius),
            "rsun_arcseconds": header.get('rsun_obs', header.get('solar_r',
                               header.get('radius',
                               constants.average_angular_size))),
            "coordinate_system": {
                'x': header.get('ctype1', 'HPLN-TAN'),
                'y': header.get('ctype2', 'HPLT-TAN')
            },
            "carrington_longitude": header.get('crln_obs', 0.),
            "heliographic_latitude": header.get('hglt_obs', 
                                     header.get('crlt_obs',
                                     header.get('solar_b0', 0.))),
            "heliographic_longitude": header.get('hgln_obs', 0.),
            "reference_coordinate": {
                'x': header.get('crval1', 0.),
                'y': header.get('crval2', 0.),
            },
            "reference_pixel": {
                'x': header.get('crpix1', (header.get('naxis1') + 1) / 2.),
                'y': header.get('crpix2', (header.get('naxis2') + 1) / 2.)
            },
            "scale": {
                'x': header.get('cdelt1', 1.),
                'y': header.get('cdelt2', 1.),
            },
            "units": {
                'x': header.get('cunit1', 'arcsec'),
                'y': header.get('cunit2', 'arcsec')
            },
            "rotation_angle": {
                'x': header.get('crota1', 0.),
                'y': header.get('crota2', 0.)
            }
        }

    def __array_finalize__(self, obj):
        """Finishes instantiation of the new map object"""
        if obj is None:
            return

        if hasattr(obj, '_original_header'):
            properties = ['_original_header', 'cmap', 'date', 'detector', 'dsun',
                          'exposure_time', 'instrument', 'measurement', 'name',
                          'observatory', 'rsun_arcseconds', 'rsun_meters',
                          'scale', 'units', 'reference_coordinate',
                          'reference_pixel', 'coordinate_system',
                          'heliographic_latitude', 'heliographic_longitude',
                          'carrington_longitude','rotation_angle']

            for attr in properties:
                setattr(self, attr, getattr(obj, attr))

    def __array_wrap__(self, out_arr, context=None):
        """Returns a wrapped instance of a Map object"""
        return np.ndarray.__array_wrap__(self, out_arr, context)

    def __getitem__(self, key):
        """Overriding indexing operation to ensure that header is updated.  Note
        that the indexing follows the ndarray row-column order, which is
        reversed from calling Map.submap()"""
        if isinstance(key, tuple):
            # Used when asking for a 2D sub-array
            # The header will be updated
            if type(key[1]) is slice:
                x_range = [key[1].start, key[1].stop]
            else:
                x_range = [key[1], key[1]+1]

            if type(key[0]) is slice:
                y_range = [key[0].start, key[0].stop]
            else:
                y_range = [key[0], key[0]+1]

            return self.submap(x_range, y_range, units="pixels")
        else:
            # Typically used by np.ndarray.__repr__() due to indexing with [-1]
            # The header will not be updated properly!
            return np.ndarray.__getitem__(self, key)

    def __add__(self, other):
        """Add two maps. Currently does not take into account the alignment
        between the two maps."""
        result = np.ndarray.__add__(self, other)

        return result

    def __repr__(self):
        if not hasattr(self, 'observatory'):
            return np.ndarray.__repr__(self)

        return (
"""SunPy Map
---------
Observatory:\t %s
Instrument:\t %s
Detector:\t %s
Measurement:\t %s
Obs Date:\t %s
dt:\t\t %f
Dimension:\t [%d, %d] 
[dx, dy] =\t [%f, %f]
 
""" % (self.observatory, self.instrument, self.detector, self.measurement,
       self.date.strftime("%Y-%m-%d %H:%M:%S"), self.exposure_time,
       self.shape[1], self.shape[0], self.scale['x'], self.scale['y']) 
     + np.ndarray.__repr__(self))

    def __sub__(self, other):
        """Subtract two maps. Currently does not take into account the
        alignment between the two maps.

        numpy dtype nums:
            1    int8
            2    uint8
            3    int16
            4    uint16
        """
        # if data is stored as unsigned, cast up (e.g. uint8 => int16)
        if self.dtype.kind == "u":
            self = self.astype(to_signed(self.dtype))
        if other.dtype.kind == "u":
            other = other.astype(to_signed(other.dtype))

        result = np.ndarray.__sub__(self, other)

        def norm():
            mean = result.mean()
            std = result.std()
            vmin = max(result.min(), mean - 6 * std)
            vmax = min(result.max(), mean + 6 * std)

            return colors.Normalize(vmin, vmax)

        result.norm = norm
        result.cmap = cm.gray  # @UndefinedVariable

        return result

    @property
    def xrange(self):
        """Return the X range of the image in arcsec from edge to edge."""
        xmin = self.center['x'] - self.shape[1] / 2 * self.scale['x']
        xmax = self.center['x'] + self.shape[1] / 2 * self.scale['x']
        return [xmin, xmax]

    @property
    def yrange(self):
        """Return the Y range of the image in arcsec from edge to edge."""
        ymin = self.center['y'] - self.shape[0] / 2 * self.scale['y']
        ymax = self.center['y'] + self.shape[0] / 2 * self.scale['y']
        return [ymin, ymax]
    
    @property
    def center(self):
        """Returns the offset between the center of the Sun and the center of 
        the map."""
        return {
            'x': wcs.get_center(self.shape[1], self.scale['x'], 
                                self.reference_pixel['x'], 
                                self.reference_coordinate['x']),
            'y': wcs.get_center(self.shape[0], self.scale['y'], 
                                self.reference_pixel['y'], 
                                self.reference_coordinate['y'])
        }

    def draw_limb(self, axes=None):
        """Draws a circle representing the solar limb 
        
            Parameters
            ----------
            axes: matplotlib.axes object or None
                Axes to plot limb on or None to use current axes.
        
            Returns
            -------
            matplotlib.axes object
        """
        
        if not axes:
            axes = plt.gca()
        
        circ = patches.Circle([0, 0], radius=self.rsun_arcseconds, fill=False,
                              color='white')
        axes.add_artist(circ)
        
        return axes

    def draw_grid(self, axes=None, grid_spacing=20):
        """Draws a grid over the surface of the Sun
        
        Parameters
        ----------
        axes: matplotlib.axes object or None
        Axes to plot limb on or None to use current axes.
        
        grid_spacing: float
            Spacing (in degrees) for logditude and latitude grid.
        
        Returns
        -------
        matplotlib.axes object
        """

        if not axes:
            axes = plt.gca()

        x, y = self.pixel_to_data()
        rsun = self.rsun_meters
        dsun = self.dsun
	
        b0 = self.heliographic_latitude
        l0 = self.heliographic_longitude
        units = [self.units.get('x'), self.units.get('y')]
	
	    #TODO: This function could be optimized. Does not need to convert the entire image
	    # coordinates
        lon_self, lat_self = wcs.convert_hpc_hg(rsun, dsun, units[0], units[1], b0, l0, x, y)
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
            x, y = wcs.convert_hg_hpc(self.rsun_meters,
                                      self.dsun, self.heliographic_latitude,
                                      self.heliographic_longitude,
                                      hg_longitude_deg_mesh,
                                      hg_latitude_deg_mesh, units='arcsec')
            axes.plot(x, y, color='white', linestyle='dotted')
            
        hg_longitude_deg = np.arange(lon_range[0], lon_range[1]+grid_spacing, grid_spacing)
        hg_latitude_deg = np.linspace(lat_range[0], lat_range[1], num=num_points)

        # draw the longitude lines
        for lon in hg_longitude_deg:
            hg_longitude_deg_mesh, hg_latitude_deg_mesh = np.meshgrid(
                lon * np.ones(num_points), hg_latitude_deg)
            x, y = wcs.convert_hg_hpc(self.rsun_meters,
                                      self.dsun, self.heliographic_latitude,
                                      self.heliographic_longitude,
                                      hg_longitude_deg_mesh,
                                      hg_latitude_deg_mesh, units='arcsec')
            axes.plot(x, y, color='white', linestyle='dotted')
            
        axes.set_ylim(self.yrange)
        axes.set_xlim(self.xrange)

        return axes

    def _validate(self):
        """Validates the meta-information associated with a Map.

        This function includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.map.sources package."""
        if (self.dsun <= 0 or self.dsun >= 40 * constants.au):
            raise InvalidHeaderInformation("Invalid value for DSUN")

    def std(self, *args, **kwargs):
        """overide np.ndarray.std()"""
        return np.array(self, copy=False, subok=False).std(*args, **kwargs)
    
    def mean(self, *args, **kwargs):
        """overide np.ndarray.mean()"""
        return np.array(self, copy=False, subok=False).mean(*args, **kwargs)
    
    def min(self, *args, **kwargs):
        """overide np.ndarray.min()"""
        return np.array(self, copy=False, subok=False).min(*args, **kwargs)
        
    def max(self, *args, **kwargs):
        """overide np.ndarray.max()"""
        return np.array(self, copy=False, subok=False).max(*args, **kwargs)

    def data_to_pixel(self, value, dim):
        """Convert pixel-center data coordinates to pixel values"""
        #TODO: This function should be renamed. It is confusing as data
        # coordinates are in something like arcsec but this function just changes how you
        # count pixels
        if dim not in ['x', 'y']:
            raise ValueError("Invalid dimension. Must be one of 'x' or 'y'.")

        size = self.shape[dim == 'x']  # 1 if dim == 'x', 0 if dim == 'y'.

        return (value - self.center[dim]) / self.scale[dim] + ((size - 1) / 2.)

    def pixel_to_data(self, x = None, y = None):
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

        scale = np.array([self.scale.get('x'), self.scale.get('y')])
        crpix = np.array([self.reference_pixel.get('x'), self.reference_pixel.get('y')])
        crval = np.array([self.reference_coordinate.get('x'), self.reference_coordinate.get('y')])
        coordinate_system = [self.coordinate_system.get('x'), self.coordinate_system.get('y')]
        x,y = wcs.convert_pixel_to_data(width, height, scale[0], scale[1], crpix[0], crpix[1], crval[0], crval[1], coordinate_system[0], x = x, y = y)

        return x, y

    def get_header(self, original=False):
        """Returns an updated MapHeader instance"""
        header = self._original_header.copy()
        
        # If requested, return original header as-is
        if original:
            return header
        
        # Bit-depth
        #
        #   8    Character or unsigned binary integer
        #  16    16-bit twos-complement binary integer
        #  32    32-bit twos-complement binary integer
        # -32    IEEE single precision floating point
        # -64    IEEE double precision floating point
        #
        if not header.has_key('bitpix'):
            bitdepth = 8 * self.dtype.itemsize
            
            if self.dtype.kind == "f":
                bitdepth = - bitdepth
                
            header['bitpix'] = bitdepth

        # naxis
        header['naxis'] = self.ndim
        header['naxis1'] = self.shape[1]
        header['naxis2'] = self.shape[0]
        
        # dsun
        if header.has_key('dsun_obs'):
            header['dsun_obs'] = self.dsun

        # rsun_obs
        if header.has_key('rsun_obs'):
            header['rsun_obs'] = self.rsun_arcseconds
        elif header.has_key('solar_r'):
            header['solar_r'] = self.rsun_arcseconds
        elif header.has_key('radius'):
            header['radius'] = self.rsun_arcseconds
            
        # cdelt
        header['cdelt1'] = self.scale['x']
        header['cdelt2'] = self.scale['y']

        # crpix
        header['crval1'] = self.reference_coordinate['x']
        header['crval2'] = self.reference_coordinate['y']
        
        # crval
        header['crpix1'] = self.reference_pixel['x']
        header['crpix2'] = self.reference_pixel['y']
        
        return header               

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
        data = resample(np.asarray(self).copy().T, dimensions,
                        method, center=True)

        # Update image scale and number of pixels
        header = self._original_header.copy()

        # Note that 'x' and 'y' correspond to 1 and 0 in self.shape,
        # respectively
        scale_factor_x = (float(self.shape[1]) / dimensions[0])
        scale_factor_y = (float(self.shape[0]) / dimensions[1])

        # Create new map instance
        new_map = self.__class__(data.T, header)

        # Update metadata
        new_map.scale['x'] *= scale_factor_x
        new_map.scale['y'] *= scale_factor_y
        new_map.reference_pixel['x'] = (dimensions[0] + 1) / 2.
        new_map.reference_pixel['y'] = (dimensions[1] + 1) / 2.
        new_map.reference_coordinate['x'] = self.center['x']
        new_map.reference_coordinate['y'] = self.center['y']

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
        reshaped = reshape_image_to_4d_superpixel(np.asarray(self).copy().T,
                                                  dimensions)
        if method == 'sum':
            data = reshaped.sum(axis=3).sum(axis=1)
        elif method == 'average':
            data = ((reshaped.sum(axis=3).sum(axis=1)) /
                    np.float32(dimensions[0] * dimensions[1]))
        
        
        #data = resample(np.asarray(self).copy().T, dimensions,
        #                method, center=True)

        # Update image scale and number of pixels
        header = self._original_header.copy()

        # Note that 'x' and 'y' correspond to 1 and 0 in self.shape,
        # respectively
        new_nx = self.shape[1] / dimensions[0]
        new_ny = self.shape[0] / dimensions[1]

        # Create new map instance
        new_map = self.__class__(data.T, header)

        # Update metadata
        new_map.scale['x'] = dimensions[0] * self.scale['x']
        new_map.scale['y'] = dimensions[1] * self.scale['y']
        new_map.reference_pixel['x'] = (new_nx + 1) / 2.
        new_map.reference_pixel['y'] = (new_ny + 1) / 2.
        new_map.reference_coordinate['x'] = self.center['x']
        new_map.reference_coordinate['y'] = self.center['y']

        return new_map
    
    def rotate(self, angle, scale=1.0, rotation_centre=None, recentre=True,
               missing=0.0, interpolation='bicubic', interp_param=-0.5):
        """Returns a new rotated, rescaled and shifted map.
        
        Parameters
        ---------
        angle: float
           The angle to rotate the image by (radians)        
        scale: float
           A scale factor for the image, default is no scaling
        rotation_centre: tuple
           The point in the image to rotate around (Axis of rotation).
           Default: Centre of the array
        recentre: bool, or array-like
           Move the centroid (axis of rotation) to the centre of the array
           or recentre coords. 
           Default: True, recentre to the centre of the array.
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
        """
        
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
        if not isinstance(recentre, bool):
            recentre = np.array(recentre).reshape(2,1)
                
        #Define Size and centre of array
        centre = (np.array(self.shape)-1)/2.0
        
        #If rotation_centre is not set (None or False),
        #set rotation_centre to the centre of the image.
        if rotation_centre is None:
            rotation_centre = centre 
        else:
            #Else check rotation_centre is a vector with shape (2,1)
            rotation_centre = np.array(rotation_centre).reshape(2,1)

        #Recentre to the rotation_centre if recentre is True
        if isinstance(recentre, bool):
            #if rentre is False then this will be (0,0)
            shift = np.array(rotation_centre) - np.array(centre) 
        else:
            #Recentre to recentre vector otherwise
            shift = np.array(recentre) - np.array(centre)
        
        image = np.asarray(self).copy()
    
        #Calulate the parameters for the affline_transform
        c = np.cos(angle)
        s = np.sin(angle)
        mati = np.array([[c, s],[-s, c]]) / scale   # res->orig
        centre = np.array([centre]).transpose()  # the centre of rotn
        shift = np.array([shift]).transpose()    # the shift
        kpos = centre - np.dot(mati, (centre + shift))  
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
                raise ValueError("You do not have the C extension sunpy.image.Crotate")
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
        header = self._original_header.copy()
        # Create new map instance
        new_map = self.__class__(data, header)
        
        return new_map
    
    def save(self, filepath):
        """Saves the SunPy Map object to a file.
        
        Currently SunPy can only save files in the FITS format. In the future
        support will be added for saving to other formats.
        
        Parameters
        ----------
        filepath : string
            Location to save file to.
        """
        pyfits_header = self.get_header().as_pyfits_header()
        hdu = pyfits.PrimaryHDU(self, header=pyfits_header)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(os.path.expanduser(filepath))        

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

        # Make a copy of the header with updated centering information
        header = self._original_header.copy()
        
        # Get ndarray representation of submap
        data = np.asarray(self)[y_pixels[0]:y_pixels[1],
                                x_pixels[0]:x_pixels[1]]
        
        # Instantiate new instance and update metadata
        new_map = self.__class__(data.copy(), header)
        new_map.reference_pixel['x'] = self.reference_pixel['x'] - x_pixels[0]
        new_map.reference_pixel['y'] = self.reference_pixel['y'] - y_pixels[0]

        return new_map

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
            axes.set_title("%s %s" % (self.name, self.date))
            
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
        
        cmap = copy(self.cmap)
        if gamma is not None:
            cmap.set_gamma(gamma)
            
            #make imshow kwargs a dict
        
        kwargs = {'origin':'lower',
                  'cmap':cmap,
                  'norm':self.norm(),
                  'extent':extent,
                  'interpolation':'nearest'}
        kwargs.update(imshow_args)
        
        ret = axes.imshow(self, **kwargs)
        
        #Set current image (makes colorbar work)
        plt.sci(ret)
        return ret
        
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
            axes = figure.add_subplot(111)

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

        plt.show()
        
        return figure
    
    def norm(self):
        """Default normalization method. Not yet implemented."""
        return None

    @classmethod
    def parse_file(cls, filepath):
        """Reads in a map file and returns a header and data array"""
        return read_file(filepath)

    @classmethod
    def read(cls, filepath):
        """Map class factory

        Attempts to determine the type of data associated with input and
        returns an instance of either the generic Map class or a subclass
        of Map such as AIAMap, EUVIMap, etc.

        Parameters
        ----------
        filepath : string
            Path to a valid FITS or JPEG 2000 file of a type supported by SunPy

        Returns
        -------
        out : Map
            Returns a Map instance for the particular type of data loaded.
        """
        data, header = cls.parse_file(filepath)

        if cls.__name__ is not "Map":
            return cls(data, header)

        for cls in Map.__subclasses__():
            if cls.is_datasource_for(header):
                return cls(data, header)
        
        return Map(data, header)

    @classmethod
    def read_header(cls, filepath):
        """Attempts to detect the datasource type and returns meta-information
        for that particular datasource."""
        header = read_file_header(filepath)

        for cls in Map.__subclasses__():
            if cls.is_datasource_for(header):
                properties = cls.get_properties(header)
                properties['header'] = header
                
                return properties

class InvalidHeaderInformation(ValueError):
    """Exception to raise when an invalid header tag value is encountered for a
    FITS/JPEG 2000 file."""
    pass

Map.create.im_func.__doc__ = Map._create.generate_docs()
