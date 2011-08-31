"""
BaseMap is a generic Map class from which all other Map classes inherit from.
"""
#pylint: disable=E1101,E1121
__authors__ = ["Keith Hughitt, Steven Christe"]
__email__ = "keith.hughitt@nasa.gov"

import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors
import matplotlib.cm as cm
from datetime import datetime
from sunpy.solwcs import solwcs as wcs

"""
Questions
---------
1. Which is better? center['x'] & center['y'] or center[0] and center[1], or?
2. map.wavelength, map.meas or? (use hv/vso/etc conventions?)
3. Are self.r_sun and radius below different? (rsun or rsun_obs for AIA?)
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
    
    def __init__(self, data, header=None): #pylint: disable=W0613
        """BaseMap constructor"""
        if header:
            self.header = header
            self.date = None
            self.name = None
            self.cmap = None
            self.norm = None
            
            # Set object attributes dynamically
            for attr, value in list(self.get_properties(header).items()):
                setattr(self, attr, value)

            self.center = {
                "x": wcs.get_center(header, axis='x'),
                "y": wcs.get_center(header, axis='y')
            }
            self.scale = {
                "x": header.get('cdelt1'),
                "y": header.get('cdelt2')
            }
            self.units = {
                "x": wcs.get_units(header, axis='x'), 
                "y": wcs.get_units(header, axis='y')
            }
            self.rsun = wcs.solar_limb(header)
            
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
                
    def __add__(self, other):
        """Add two maps. Currently does not take into account the alignment  
        between the two maps."""
        result = np.ndarray.__add__(self, other)
        
        return result
    
    def __getitem__(self, key):
        """Overiding indexing operation to ensure that header is updated"""
        if isinstance(key, tuple):
            x_range = [key[1].start, key[1].stop]
            y_range = [key[0].start, key[0].stop]

            return self.submap(x_range, y_range, units="pixels")
        else:
            return np.ndarray.__getitem__(self, key)
    
    def __sub__(self, other):
        """Subtract two maps. Currently does not take into account the alignment 
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
        xmin = self.center['x'] - self.shape[1] / 2 * self.scale['x']
        xmax = self.center['x'] + self.shape[1] / 2 * self.scale['x']
        return [xmin, xmax]
        
    def get_yrange(self):
        """Return the Y range of the image in arcsec."""
        ymin = self.center['y'] - self.shape[0] / 2 * self.scale['y']
        ymax = self.center['y'] + self.shape[0] / 2 * self.scale['y']
        return [ymin, ymax]

    def submap(self, x_range, y_range, units="arcseconds"):
        """Returns a submap of the map with the specified range
        
        Keith [08/19/2011]
         * Slicing in numpy expects [y, x]. Should we break convention here? 
        """
        height = self.shape[0]
        width = self.shape[1]
        
        # Arcseconds => Pixels
        #  x_px = (x / cdelt1) + (width / 2)
        #
        if units is "arcseconds":
            x_pixels = (np.array(x_range) / self.scale['x']) + (width / 2)
            y_pixels = (np.array(y_range) / self.scale['y']) + (height / 2)

        elif units is "pixels":
            x_pixels = x_range
            y_pixels = y_range

        # Make a copy of the header with updated centering information        
        header = copy.deepcopy(self.header)
        header['crpix1'] = header['crpix1'] - x_pixels[0]
        header['crpix2'] = header['crpix2'] - y_pixels[0]
        header['naxis1'] = x_pixels[1] - x_pixels[0]
        header['naxis2'] = y_pixels[1] - y_pixels[0]
        
        # Get ndarray representation of submap
        data = np.asarray(self)[y_pixels[0]:y_pixels[1], 
                                x_pixels[0]:x_pixels[1]]

        return self.__class__(data, header)
   
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
            circ = patches.Circle([0, 0], radius=self.rsun, 
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

class UnrecognizedDataSouceError(ValueError):
    """Exception to raise when an unknown datasource is encountered"""
    pass
