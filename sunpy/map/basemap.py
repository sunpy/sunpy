"""
BaseMap is a generic Map class from which all other Map classes inherit from.
"""
#pylint: disable=E1101,E1121
__authors__ = ["Keith Hughitt, Steven Christe"]
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors
import matplotlib.cm as cm
from datetime import datetime
from sunpy.sun import sun

"""
Questions
---------
1. Which is better? center['x'] & center['y'] or center[0] and center[1], or?
2. map.wavelength, map.meas or? (use hv/vso/etc conventions?)
3. Are self.r_sun and radius below different? (rsun or rsun_obs for AIA?)
4. Should default cmap and normalization settings be chosen for each image?

Requests
--------
1. Would be nice to be able to easily extract the data as a ndarray from the map.
2. Need to provide a way to get a sub map.
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
    header : dict
        A dictionary representation of the image header
    date : datetime
        Image observation time
    det : str
        Detector name
    inst : str
        Instrument name
    meas : str, int
        Measurement name. For AIA this is the wavelength of image
    obs : str
        Observatory name
    r_sun : float
        Radius of the sun
    name : str
        Nickname for the image type (e.g. "AIA 171")
    center : dict
        X and Y coordinate for the center of the sun in units
    scale: dict
        Image scale along the x and y axes in units/pixel
    units: dict
        Image coordinate units along the x and y axes

    Examples
    --------
    >>> aia = sunpy.Map(sunpy.AIA_171_IMAGE)
    >>> aia.T
    Map([[ 0.3125,  1.    , -1.1875, ..., -0.625 ,  0.5625,  0.5   ],
    [-0.0625,  0.1875,  0.375 , ...,  0.0625,  0.0625, -0.125 ],
    [-0.125 , -0.8125, -0.5   , ..., -0.3125,  0.5625,  0.4375],
    ..., 
    [ 0.625 ,  0.625 , -0.125 , ...,  0.125 , -0.0625,  0.6875],
    [-0.625 , -0.625 , -0.625 , ...,  0.125 , -0.0625,  0.6875],
    [ 0.    ,  0.    , -1.1875, ...,  0.125 ,  0.    ,  0.6875]])
    >>> aia.header.get('cunit1')
    'arcsec'
    >>> aia.plot()
    >>> import matplotlib.cm as cm
    >>> import matplotlib.colors as colors
    >>> aia.plot(cmap=cm.hot, norm=colors.Normalize(1, 2048))
    
    See Also:
    ---------
    numpy.ndarray Parent class for the Map object
    
    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    | http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    | http://docs.scipy.org/doc/numpy/reference/ufuncs.html
    | http://www.scipy.org/Subclasses

    """
    def __new__(cls, data, header=None): #pylint: disable=W0613
        """Creates a new BaseMap instance"""        
        if isinstance(data, np.ndarray):
            obj = data.view(cls)
        elif isinstance(data, list):
            obj = np.asarray(data).view(cls)
        
        return obj
    
    def __init__(self, data, header=None):
        """BaseMap constructor"""
        if header:
            self.header = header
            self.date = None
            self.name = None
            self.cmap = None
            self.norm = None
            #self.exptime = None
            
            # Set object attributes dynamically
            for attr, value in list(self.get_properties(header).items()):
                setattr(self, attr, value)

            self.center = {
                "x": header.get('cdelt1')*data.shape[0]/2 + header.get('crval1') - header.get('crpix1')*header.get('cdelt1'),
                "y": header.get('cdelt2')*data.shape[1]/2 + header.get('crval2') - header.get('crpix2')*header.get('cdelt2')
            }
            self.scale = {
                "x": header.get('cdelt1'),
                "y": header.get('cdelt2')
            }
            self.units = {
                "x": header.get('cunit1'), 
                "y": header.get('cunit2')
            }
                
    def __add__(self, other):
        """Add two maps. Currently does not take into account the alignment  
        between the two maps."""
        result = np.ndarray.__add__(self, other)
        
        return result
    
    def __sub__(self, other):
        """Add two maps. Currently does not take into account the alignment 
        between the two maps."""
        result = np.ndarray.__sub__(self, other)

        minmax = np.array([abs(result.min()), abs(result.max())]).max()
        result.norm = colors.Normalize(-minmax, minmax, True)
        
        result.cmap = cm.gray #@UndefinedVariable 
        
        return result
    
    @classmethod
    def get_properties(cls):
        """Returns default map properties""" 
        return {
            'cmap': cm.gray,  #@UndefinedVariable
            'norm': colors.Normalize(5, 1024, True),
            'date': datetime.today(),
            'det': "None",
            'inst': "None",
            'meas': "None",
            'obs': "None",
            'name': "Default Map",
            'r_sun': None
        }
        
    def get_xrange(self):
        """Return the X range of the image in arcsec."""        
        xmin = self.center['x'] - self.shape[0]/2*self.scale['x']
        xmax = self.center['x'] + self.shape[0]/2*self.scale['x']
        return [xmin,xmax]
        
    def get_yrange(self):
        """Return the Y range of the image in arcsec."""
        ymin = self.center['y'] - self.shape[1]/2*self.scale['y']
        ymax = self.center['y'] + self.shape[1]/2*self.scale['y']
        return [ymin,ymax]
    
    def _transform_pixel_to_coord(self, pixel, axis):    
        """Given a pixel number return the coordinate value of that pixel."""
        pixels = np.array(pixel)
        if axis == 'x': n = self.shape[0]
        if axis == 'y': n = self.shape[1]
        return self.scale[axis] * pixels + self.center[axis] - n / 2.0 * self.scale[axis]
        
    def _transform_coord_to_pixel(self, coord, axis):
        """Given a data coordinate return the pixel value (not necessarily an integer)."""
        coords = np.array(coord)
        if axis == 'x': n = self.shape[0]
        if axis == 'y': n = self.shape[1]
        return 1 / self.scale[axis] * (coords + n / 2.0 * self.scale[axis] - self.center[axis])

    def submap(self, x_range, y_range):        
        #convert data coordinates to pixel coordinates        
        xrange_pixelcoord = self._transform_coord_to_pixel(x_range, 'x')
        yrange_pixelcoord = self._transform_coord_to_pixel(y_range, 'y')

        xrange_pixelcoord = xrange_pixelcoord.astype('int')
        yrange_pixelcoord = yrange_pixelcoord.astype('int')
        
        dpixel = [0.5*(xrange_pixelcoord[1] - xrange_pixelcoord[0]), 0.5*(yrange_pixelcoord[1] - yrange_pixelcoord[0])]
        
        xcenter_datacoord = self._transform_pixel_to_coord(xrange_pixelcoord[0] + dpixel[0], 'x')
        ycenter_datacoord = self._transform_pixel_to_coord(yrange_pixelcoord[0] + dpixel[1], 'y')
        
        self.center['x'] = xcenter_datacoord
        self.center['y'] = ycenter_datacoord
        
        self = self[xrange_pixelcoord[0]:xrange_pixelcoord[1], yrange_pixelcoord[0]:yrange_pixelcoord[1]]
        return self
        
    def plot(self, draw_limb=True, **matplot_args):
        """Plots the map object using matplotlib
        
        Parameters
        ----------
        draw_limb : bool
            Whether a circle should be drawn around the solar limb.
        **matplot_args : dict
            Matplotlib Any additional im_show arguments that should be used
            when plotting the image.
        """
        # Create a figure and add title and axes
        fig = plt.figure()
        
        axes = fig.add_subplot(111)
        axes.set_title("%s %s" % (self.name, self.date))
        axes.set_xlabel('X-postion [' + self.units['x'] + ']')
        axes.set_ylabel('Y-postion [' + self.units['y'] + ']')
        
        # Draw circle at solar limb
        if draw_limb:
            circ = patches.Circle([0, 0], radius=sun.radius(self.date), 
                fill=False, color='white')
            axes.add_artist(circ)

        # Determine extent
        extent = self.get_xrange() + self.get_yrange()
        
        # Matplotlib arguments
        params = {
            "cmap": self.cmap,
            "norm": self.norm
        }
        params.update(matplot_args)
            
        plt.imshow(self, origin='lower', extent=extent, **params)
        plt.colorbar()
        plt.show()

    def __array_finalize__(self, obj):
        """Finishes instantiation of the new map object"""
        if obj is None:
            return

        if hasattr(obj, 'header'):
            self.header = obj.header

            # preserve object properties
            properties = self.get_properties(obj.header)
            for attr, value in list(properties.items()):
                setattr(self, attr, getattr(obj, attr, value))
                
            self.center = obj.center
            self.scale = obj.scale
            self.units = obj.units
        
    def __array_wrap__(self, out_arr, context=None):
        """Returns a wrapped instance of a Map object"""
        return np.ndarray.__array_wrap__(self, out_arr, context)

class UnrecognizedDataSouceError(ValueError):
    """Exception to raise when an unknown datasource is encountered"""
    pass
