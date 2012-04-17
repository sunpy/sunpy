"""
BaseMap is a generic Map class from which all other Map classes inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404
__authors__ = ["Keith Hughitt, Steven Christe"]
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import colors
from matplotlib import cm
from copy import copy
from sunpy.wcs import wcs as wcs
from sunpy.util.util import toggle_pylab
from sunpy.io import read_file
from sunpy.sun import constants
from sunpy.time import parse_time
from sunpy.map.header import MapHeader

"""
TODO
----
* Automatically include BaseMap docstring when displaying help for subclasses?

Questions
---------
* Should we use Helioviewer or VSO's data model? (e.g. map.meas, map.wavelength
or something else?)
* Should 'center' be renamed to 'offset' and crpix1 & 2 be used for 'center'?
"""

class BaseMap(np.ndarray):
    """
    BaseMap(data, header)

    A spatially-aware data array based on the SolarSoft Map object

    Parameters
    ----------
    data : numpy.ndarray, list
        A 2d list or ndarray containing the map data
    header : dict
        A dictionary of the original image header tags

    Attributes
    ----------
    fits_header : dict
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
    >>> aia.fits_header.get('cunit1')
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
    def __new__(cls, data):  # pylint: disable=W0613
        """Creates a new BaseMap instance"""
        if isinstance(data, np.ndarray):
            obj = data.view(cls)
        elif isinstance(data, list):
            obj = np.asarray(data).view(cls)
        else:
            raise TypeError('Invalid input')

        return obj

    def __init__(self, header):
        self.fits_header = header

        # Set properties
        self.byte_scaled = self.dtype == np.uint8
        self.cmap = cm.gray  # @UndefinedVariable
        self.date = parse_time(header.get('date-obs'))
        self.detector = header.get('detector')
        self.dsun = header.get('dsun_obs', constants.au)
        self.exposure_time = header.get('exptime')
        self.instrument = header.get('instrume')
        self.measurement = header.get('wavelnth')
        self.observatory = header.get('telescop')
        self.name = header.get('telescop') + " " + str(header.get('wavelnth'))
        self.rsun_meters = header.get('RSUN_REF', constants.radius)
        self.rsun_arcseconds = header.get('rsun_obs', header.get('solar_r',
                               header.get('radius',
                               constants.average_angular_size)))
        
        self.coordinate_system = {
            'x': header.get('ctype1'),
            'y': header.get('ctype2')
        }
        
        self.carrington_longitude = header.get('crln_obs', 0)
        
        self.heliographic_latitude = header.get('hglt_obs', 
                                     header.get('crlt_obs',
                                     header.get('solar_b0', 0)))
        
        self.heliographic_longitude = header.get('hgln_obs', 0)

        self.reference_coordinate = {
            'x': header.get('crval1', 0),
            'y': header.get('crval2', 0),
        }

        self.reference_pixel = {
            'x': header.get('crpix1', (self.shape[0] + 1) / 2.),
            'y': header.get('crpix2', (self.shape[1] + 1) / 2.)
        }

        self.scale = {
            'x': header.get('cdelt1'),
            'y': header.get('cdelt2'),
        }

        self.units = {
            'x': header.get('cunit1', 'arcsec'),
            'y': header.get('cunit2', 'arcsec')
        }
        
        # Validate properties
        self._validate()

    def __array_finalize__(self, obj):
        """Finishes instantiation of the new map object"""
        if obj is None:
            return

        if hasattr(obj, 'header'):
            properties = ['header', 'cmap', 'date', 'detector', 'dsun',
                          'exposure_time', 'instrument', 'measurement', 'name',
                          'observatory', 'rsun_arcseconds', 'rsun_meters',
                          'scale', 'units', 'reference_coordinate', 'center',
                          'reference_pixel', 'coordinate_system',
                          'heliographic_latitude', 'heliographic_longitude',
                          'carrington_longitude']

            for attr in properties:
                setattr(self, attr, getattr(obj, attr))

    def __array_wrap__(self, out_arr, context=None):
        """Returns a wrapped instance of a Map object"""
        return np.ndarray.__array_wrap__(self, out_arr, context)

    def __getitem__(self, key):
        """Overiding indexing operation to ensure that header is updated"""
        if isinstance(key, tuple) and type(key[0]) is slice:
            x_range = [key[1].start, key[1].stop]
            y_range = [key[0].start, key[0].stop]

            return self.submap(y_range, x_range, units="pixels")
        else:
            return np.ndarray.__getitem__(self, key)

    def __add__(self, other):
        """Add two maps. Currently does not take into account the alignment
        between the two maps."""
        result = np.ndarray.__add__(self, other)

        return result

    def __repr__(self):
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
       self.shape[0], self.shape[1], self.scale['x'], self.scale['y']) 
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
            dtype = "int%d" % (int(self.dtype.name[4:]) * 2)
            self = self.astype(np.dtype(dtype))
        if other.dtype.kind == "u":
            dtype = "int%d" % (int(other.dtype.name[4:]) * 2)
            other = other.astype(np.dtype(dtype))

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
            'x': wcs.get_center(self.shape[0], self.scale['x'], 
                                self.reference_pixel['x'], 
                                self.reference_coordinate['x']),
            'y': wcs.get_center(self.shape[0], self.scale['y'], 
                                self.reference_pixel['y'], 
                                self.reference_coordinate['y'])
        }

    def _draw_limb(self, fig, axes):
        """Draws a circle representing the solar limb"""
        circ = patches.Circle([0, 0], radius=self.rsun_arcseconds, fill=False,
                              color='white')
        axes.add_artist(circ)
        return fig, axes

    def _draw_grid(self, fig, axes, grid_spacing=20):
        """Draws a grid over the surface of the Sun"""
        # define the number of points for each latitude or longitude line
        num_points = 20
        hg_longitude_deg = np.linspace(-90, 90, num=num_points)
        hg_latitude_deg = np.arange(-90, 90, grid_spacing)

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

        hg_longitude_deg = np.arange(-90, 90, grid_spacing)
        hg_latitude_deg = np.linspace(-90, 90, num=num_points)

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

        return fig, axes

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

    def data_to_pixel(self, value, dim):
        """Convert pixel-center data coordinates to pixel values"""
        if dim not in ['x', 'y']:
            raise ValueError("Invalid dimension. Must be one of 'x' or 'y'.")

        size = self.shape[dim == 'x']  # 1 if dim == 'x', 0 if dim == 'y'.

        return (value - self.center[dim]) / self.scale[dim] + ((size - 1) / 2.)

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
        from sunpy.image import resample

        # Note: because the underlying ndarray is transposed in sense when
        #   compared to the Map, the ndarray is transposed, resampled, then
        #   transposed back
        # Note: "center" defaults to True in this function because data
        #   coordinates in a Map are at pixel centers

        # Make a copy of the original data and perform resample
        data = resample(np.asarray(self).copy().T, dimensions,
                        method, center=True)

        # Update image scale and number of pixels
        header = self.fits_header.copy()

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
        new_map.reference_coordinate['y'] = self.center['x']

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
            #x_pixels = [self.data_to_pixel(elem, 'x') for elem in range_a]
            x_pixels = [np.ceil(self.data_to_pixel(range_a[0], 'x')),
                        np.floor(self.data_to_pixel(range_a[1], 'x')) + 1]
            #y_pixels = [self.data_to_pixel(elem, 'y') for elem in range_b]
            y_pixels = [np.ceil(self.data_to_pixel(range_b[0], 'y')),
                        np.floor(self.data_to_pixel(range_b[1], 'y')) + 1]
        elif units is "pixels":
            x_pixels = range_a
            y_pixels = range_b
        else:
            raise ValueError(
                "Invalid unit. Must be one of 'data' or 'pixels'")

        # Make a copy of the header with updated centering information
        header = self.fits_header.copy()
        
        # Get ndarray representation of submap
        data = np.asarray(self)[y_pixels[0]:y_pixels[1],
                                x_pixels[0]:x_pixels[1]]
        
        # Instantiate new instance and update metadata
        new_map = self.__class__(data.copy(), header)
        new_map.reference_pixel['x'] = self.reference_pixel['x'] - x_pixels[0]
        new_map.reference_pixel['y'] = self.reference_pixel['y'] - y_pixels[0]

        return new_map

    @toggle_pylab
    def plot(self, figure=None, overlays=None, draw_limb=True, gamma=None,
             draw_grid=False, colorbar=True, basic_plot=False, **matplot_args):
        """Plots the map object using matplotlib

        Parameters
        ----------
        overlays : list
            List of overlays to include in the plot
        draw_limb : bool
            Whether the solar limb should be plotted.
        draw_grid : bool
            Whether solar meridians and parallels
        grid_spacing : float
            Set the spacing between meridians and parallels for the grid
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
        if overlays is None:
            overlays = []
        if draw_limb:
            overlays = overlays + [self._draw_limb]
        # TODO: need to be able to pass the grid spacing to _draw_grid from the
        # plot command.
        if draw_grid:
            overlays = overlays + [self._draw_grid]

        # Create a figure and add title and axes
        if figure is None:
            figure = plt.figure(frameon=not basic_plot)

        # Basic plot
        if basic_plot:
            axes = plt.Axes(figure, [0., 0., 1., 1.])
            axes.set_axis_off()
            figure.add_axes(axes)
            
        # Normal plot
        else:
            axes = figure.add_subplot(111)
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

        # Matplotlib arguments
        params = {
            "cmap": self.cmap,
            "norm": self.norm()
        }
        params.update(matplot_args)

        if gamma is not None:
            params['cmap'] = copy(params['cmap'])
            params['cmap'].set_gamma(gamma)

        im = axes.imshow(self, origin='lower', extent=extent, **params)
        
        if colorbar and not basic_plot:
            figure.colorbar(im)

        for overlay in overlays:
            figure, axes = overlay(figure, axes)
        return figure

    def show(self, figure=None, overlays=None, draw_limb=False, gamma=1.0,
             draw_grid=False, colorbar=True, basic_plot=False, **matplot_args):
        """Displays map on screen. Arguments are same as plot()."""
        self.plot(figure, overlays, draw_limb, gamma, draw_grid, colorbar, 
                  basic_plot, **matplot_args).show()
    
    def norm(self):
        """Default normalization method"""
        return None

    @classmethod
    def parse_file(cls, filepath):
        """Reads in a map file and returns a header and data array"""
        data, dict_header = read_file(filepath)

        header = MapHeader(dict_header)

        return header, data

    @classmethod
    def read(cls, filepath):
        """Map class factory

        Attempts to determine the type of data associated with input and
        returns an instance of either the generic BaseMap class or a subclass
        of BaseMap such as AIAMap, EUVIMap, etc.

        Parameters
        ----------
        filepath : string
            Path to a valid FITS or JPEG 2000 file of a type supported by SunPy

        Returns
        -------
        out : Map
            Returns a Map instance for the particular type of data loaded.
        """
        header, data = cls.parse_file(filepath)

        if cls.__name__ is not "BaseMap":
            return cls(data, header)

        for cls in BaseMap.__subclasses__():
            if cls.is_datasource_for(header):
                return cls(data, header)
        raise UnrecognizedDataSouceError("File header not recognized by SunPy")

#    @classmethod
#    def detect_properties(cls, filepath):
#        """Attempts to detect the datasource type and returns meta-information
#        for that particular datasource."""
#        dict_header = read_header(filepath)
#
#        header = MapHeader(dict_header)
#
#        for cls in BaseMap.__subclasses__():
#            if cls.is_datasource_for(header):
#                return cls.get_properties(header)
#


class UnrecognizedDataSouceError(ValueError):
    """Exception to raise when an unknown datasource is encountered"""
    pass


class InvalidHeaderInformation(ValueError):
    """Exception to raise when an invalid header tag value is encountered for a
    FITS/JPEG 2000 file."""
    pass
