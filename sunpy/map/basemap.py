from __future__ import absolute_import

"""
BaseMap is a generic Map class from which all other Map classes inherit from.
"""

#pylint: disable=E1101,E1121
__authors__ = ["Keith Hughitt, Steven Christe"]
__email__ = "keith.hughitt@nasa.gov"

from copy import copy

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
1. map.wavelength, map.meas or? (use hv/vso/etc conventions?)
2. Are self.r_sun and radius below different? (rsun or rsun_obs for AIA?)
"""


def draw_limb(map_, fig, axes):
    circ = patches.Circle([0, 0], radius=map_.rsun, 
        fill=False, color='white')
    axes.add_artist(circ)
    return fig, axes
_draw_limb = draw_limb


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
    rsun : float
        Radius of the sun
    exptime: float
        Exposure time of the image in seconds.
    name : str
        Nickname for the image type (e.g. "AIA 171")
    center : dict
        X and Y coordinate of the center of the map in units. Usually represents the offset
        between the center of the Sun and the center of the map.
    scale: dict
        Image scale along the x and y axes in units/pixel
    units: dict
        Image coordinate units along the x and y axes

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
    >>> aia.header.get('cunit1')
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
    def __new__(cls, data, header=None): #pylint: disable=W0613
        """Creates a new BaseMap instance"""        
        if isinstance(data, np.ndarray):
            obj = data.view(cls)
        elif isinstance(data, list):
            obj = np.asarray(data).view(cls)
        else:
            raise TypeError('Invalid data')
        
        return obj
    
    def __init__(self, data, header=None): #pylint: disable=W0613
        """BaseMap constructor"""
        if header:
            self.header = header
            self.date = None
            self.name = None
            self.cmap = None
            self.norm = None
            self.exptime = None
            
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
            self.rsun = wcs.get_solar_limb(header)
            
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
        if isinstance(key, tuple) and type(key[0]) is slice:
            x_range = [key[1].start, key[1].stop]
            y_range = [key[0].start, key[0].stop]

            return self.submap(y_range, x_range, units="pixels")
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
            'date': datetime.today(),
            'det': "None",
            'inst': "None",
            'meas': "None",
            'obs': "None",
            'exptime': "None", 
            'name': "Default Map",
            'r_sun': None
        }
    
    @property
    def xrange(self):
        """Return the X range of the image in arcsec."""        
        xmin = self.center['x'] - self.shape[1] / 2 * self.scale['x']
        xmax = self.center['x'] + self.shape[1] / 2 * self.scale['x']
        return [xmin, xmax]
    
    @property
    def yrange(self):
        """Return the Y range of the image in arcsec."""
        ymin = self.center['y'] - self.shape[0] / 2 * self.scale['y']
        ymax = self.center['y'] + self.shape[0] / 2 * self.scale['y']
        return [ymin, ymax]
    
    def arcsecs_to_pixels(self, value, dim):
        if dim not in ['x', 'y']:
            raise ValueError("Invalid dimension. Must be one of 'x' or 'y'.")
        size = self.shape[dim == 'y'] # 0 if dim == 'x', 1 if dim == 'y'.
        
        return value / self.scale[dim] + (size / 2)

    def submap(self, range_a, range_b, units="arcseconds"):
        """Returns a submap of the map with the specified range
        
        Parameters
        ----------
        range_a : list
            The range of data to select across either the x axis (if
            units='arcseconds') or the y axis (if units='pixels').
        range_b : list
            The range of data to select across either the y axis (if
            units='arcseconds') or the x axis (if units='pixels').
        units : {'arcseconds' | 'pixels'}, optional
            The units for which the submap region has been specified.
            
        Returns
        -------
        out : BaseMap
            A new map instance is returned representing to specified sub-region.
        
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
        # Arcseconds => Pixels
        #  x_px = (x / cdelt1) + (width / 2)
        #
        if units is "arcseconds":
            x_pixels = [self.arcsecs_to_pixels(elem, 'x') for elem in range_a]
            y_pixels = [self.arcsecs_to_pixels(elem, 'y') for elem in range_b]
        elif units is "pixels":
            x_pixels = range_b
            y_pixels = range_a
        else:
            raise ValueError(
                "Invalid unit. Must be one of 'arcseconds' or 'pixels'")

        # Make a copy of the header with updated centering information        
        header = self.header.copy()
        header['crpix1'] = header['crpix1'] - x_pixels[0]
        header['crpix2'] = header['crpix2'] - y_pixels[0]
        header['naxis1'] = x_pixels[1] - x_pixels[0]
        header['naxis2'] = y_pixels[1] - y_pixels[0]
        
        # Get ndarray representation of submap
        data = np.asarray(self)[y_pixels[0]:y_pixels[1], 
                                x_pixels[0]:x_pixels[1]]

        return self.__class__(data, header)
   
    def plot(self, overlays=[], draw_limb=False, gamma=None, **matplot_args):
        """Plots the map object using matplotlib
        
        Parameters
        ----------
        draw_limb : bool
            Whether the solar limb should be plotted.
        gamma : float
            Gamma value to use for the color map
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting the image.
        """
        if draw_limb:
            overlays = overlays + [_draw_limb]
        # Create a figure and add title and axes
        fig = plt.figure()
        
        axes = fig.add_subplot(111)
        axes.set_title("%s %s" % (self.name, self.date))
        axes.set_xlabel('X-position [' + self.units['x'] + ']')
        axes.set_ylabel('Y-position [' + self.units['y'] + ']')

        # Determine extent
        extent = self.xrange + self.yrange

        # Matplotlib arguments
        params = {
            "cmap": self.cmap,
            "norm": self.norm
        }
        params.update(matplot_args)
        if gamma is not None:
            params['cmap'] = copy(params['cmap'])
            params['cmap'].set_gamma(gamma) 
        im = axes.imshow(self, origin='lower', extent=extent, **params)
        fig.colorbar(im)
        
        for overlay in overlays:
            fig, axes = overlay(self, fig, axes)
        return fig
    
    def show(self, overlays=[], draw_limb=False, gamma=1.0, **matplot_args):
        """Displays map on screen. Arguments are same as plot()."""
        self.plot(overlays, draw_limb, gamma, **matplot_args).show()


class UnrecognizedDataSouceError(ValueError):
    """Exception to raise when an unknown datasource is encountered"""
    pass
