"""A Generic Python Map Object

Authors: `Keith Hughitt <keith.hughitt@nasa.gov>`_ and `Steven Christe 
<steven.d.christe@nasa.gov>`_

Questions
---------
1. Which is better? centerX & centerY or center[0] and center[1], or?
2. map.wavelength, map.meas or? (use hv/vso/etc conventions?)
3. Are self.r_sun and radius below different? (rsun or rsun_obs for AIA?)
4. Should default cmap and normalization settings be chosen for each image?

References
----------
| http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
| http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
| http://www.scipy.org/Subclasses
"""
__author__ = "Keith Hughitt and Steven Christe"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sunpy import Sun

class BaseMap(np.ndarray):
    """
    BaseMap(data, header)
    
    A spatially-aware data array based on the SolarSoft Map object

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
    centerX : float
        X-coordinate for the center of the sun in arcseconds
    centerY : float
        Y-coordinate for the center of the sun in arcseconds
    scaleX : float
        Image scale along the x-axis in arcseconds/pixel
    scaleY : float
        Image scale along the y-axis in arcseconds/pixel

    Parameters
    ----------
    data : numpy.ndarray, list
        A 2d list or ndarray containing the map data
    header : dict
        A dictionary of the original image header tags
        
    See Also:
    ---------
    numpy.ndarray Parent class for the Map object
        
    Examples
    --------
    >>> map = sunpy.Map('sunpy/dev/sample-data/AIA20110319_105400_0171.fits')
    >>> map.T
    Map([[ 0.3125,  1.    , -1.1875, ..., -0.625 ,  0.5625,  0.5   ],
    [-0.0625,  0.1875,  0.375 , ...,  0.0625,  0.0625, -0.125 ],
    [-0.125 , -0.8125, -0.5   , ..., -0.3125,  0.5625,  0.4375],
    ..., 
    [ 0.625 ,  0.625 , -0.125 , ...,  0.125 , -0.0625,  0.6875],
    [-0.625 , -0.625 , -0.625 , ...,  0.125 , -0.0625,  0.6875],
    [ 0.    ,  0.    , -1.1875, ...,  0.125 ,  0.    ,  0.6875]])
    >>> map.header['cunit1']
    'arcsec'
    >>> map.plot()
    >>> map.plot(cmap=cm.hot, norm=colors.Normalize(1, 2048))
    """
    def __new__(cls, data, header=None):
        """Creates a new BaseMap instance"""        
        if isinstance(data, np.ndarray):
            obj = data.view(cls)
        elif isinstance(data, list):
            obj = np.asarray(data).view(cls)

        if header:
            obj.header = header
            obj.centerX = header['crpix1']
            obj.centerY = header['crpix2']
            obj.scaleX = header['cdelt1']
            obj.scaleY = header['cdelt2']

        return obj
            
    def __array_finalize__(self, obj):
        """Finishes instantiation of the new map object"""
        if obj is None: return
        
    def plot(self, cmap=None, norm=None, draw_limb=True):
        """Plots the map object using matplotlib
        
        Parameters
        ----------
        cmap : matplotlib.colors.Colormap
            Colormap to apply to the image. Defaults to a adaptive logarithmic
            grayscale colormap.
        norm : matplotlib.colors.Normalize
            Normalization method to use on the data when plotting
        draw_limb : bool
            Whether a circle should be drawn around the solar limb
        """
        # Create a figure and add title and axes
        fig = plt.figure()
        
        ax = fig.add_subplot(111)
        ax.set_title("%s %s" % (self.name, self.date))
        ax.set_xlabel('X-postion (arcseconds)')
        ax.set_ylabel('Y-postion (arcseconds)')
        
        # Draw circle at solar limb
        if draw_limb:
            circ = patches.Circle([0, 0], radius=Sun.radius(self.date), 
                fill=False, color='white')
            ax.add_artist(circ)

        # Determine extent
        xmin = -(self.centerX - 1) * self.scaleX
        xmax = (self.centerX - 1) * self.scaleX
        ymin = -(self.centerY - 1) * self.scaleY
        ymax = (self.centerY - 1) * self.scaleY
        
        extent = [xmin, xmax, ymin, ymax]
        
        # Use default cmap and norm values if none set
        if cmap is None:
            cmap = self.cmap
        if norm is None:
            norm = self.norm
            
        plt.imshow(self, cmap=cmap, origin='lower', extent=extent, norm=norm)
        plt.colorbar()
        plt.show()      

