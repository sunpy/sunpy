"""Definition for a Python Map Object

| Author: Keith Hughitt <keith.hughitt@nasa.gov>
| Author: Steven Christe <steven.d.christe@nasa.gov>

Questions
---------
    1. what is better? centerX & centerY or center[0] and center[1], or?
    2. map.wavelength, map.meas or? (use hv/vso/etc conventions?)
    3. Are self.r_sun and radius below different?
    4. Should default cmap and normalization settings be chosen for each image?

References
----------
| http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
| http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
| http://www.scipy.org/Subclasses
"""

import numpy as np
import pyfits
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.patches as patches
import Sun
from datetime import datetime

class Map(np.ndarray):
    """Map(input)
    
        A spatially-aware data array based on the SolarSoft Map object
    
        Attributes
        ----------
        header : dict
            A dictionary representation of the image header
        date : datetime
            Image observation time
        det: str
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
        input_ : filepath, data array
            The data source used to create the map object
            
        See Also:
        ---------
        numpy.ndarray Parent class for the Map object
            
        Examples
        --------
        >>>map = sunpy.Map('sunpy/dev/sample-data/AIA20110319_105400_0171.fits')
        >>>map.T
        Map([[ 0.3125,  1.    , -1.1875, ..., -0.625 ,  0.5625,  0.5   ],
        [-0.0625,  0.1875,  0.375 , ...,  0.0625,  0.0625, -0.125 ],
        [-0.125 , -0.8125, -0.5   , ..., -0.3125,  0.5625,  0.4375],
        ..., 
        [ 0.625 ,  0.625 , -0.125 , ...,  0.125 , -0.0625,  0.6875],
        [-0.625 , -0.625 , -0.625 , ...,  0.125 , -0.0625,  0.6875],
        [ 0.    ,  0.    , -1.1875, ...,  0.125 ,  0.    ,  0.6875]])
        >>>map.header['cunit1']
        'arcsec'
        >>>map.plot()
        >>>map.plot(cmap=cm.hot, norm=colors.Normalize(1, 2048))
    """
    def __new__(cls, input_):
        """Creates a new Map instance"""
        if isinstance(input_, str):
            try:
                fits = pyfits.open(input_)
            except IOError:
                sys.exit("Unable to read the file %s" % input_)

            #np.ndarray.__new__(self, fits[0].data)
            obj = np.asarray(fits[0].data).view(cls)
            
            obj.header = fits[0].header
            
            obj.centerX = obj.header['crpix1']
            obj.centerY = obj.header['crpix2']
            obj.scaleX = obj.header['cdelt1']
            obj.scaleY = obj.header['cdelt2']
            
            norm_header = parse_header(obj.header)
            
            obj.cmap = norm_header['cmap']
            obj.norm = norm_header['norm']
            obj.date = norm_header['date']
            obj.det = norm_header['det']
            obj.inst = norm_header['inst']
            obj.meas = norm_header['meas']
            obj.obs = norm_header['obs']
            obj.name = norm_header['name']
            obj.r_sun = norm_header['r_sun']

        elif isinstance(input_, list):
            obj = np.asarray(input_).view(cls)
        elif isinstance(input_, np.ndarray):
            obj = input_
            
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
        ax.set_title("%s %s" % (self.inst, self.date))
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

#
# External functions - These should be moved to a separate file or files
# once an appropriate location is determined
#         
def parse_header(header):
    """Parses a FITS, etc image header
    
    Attempts to detect the type of data (e.g. AIA) based on values in the 
    file header and returns a mapped dictionary of of some important values.

    Parameters
    ----------
    header : dict
        A dictionary container the header keywords from the file being read in
    """
    # AIA
    if header['telescop'] == 'SDO/AIA':
        return get_norm_header_tags(header, "aia")
        
    # HMI
    elif header['telescop'] == 'SDO/HMI':
        return get_norm_header_tags(header, "hmi")

def get_norm_header_tags(header, type_):
    """Returns a normalized dictionary of header values
    
    A normalized mapping of important header values is created and returned.
    Not all of the header values are used, but instead only those that are
    required for the Map class to function are included. Note that some values
    may be cast to new types in the process.
    
    Parameters
    ----------
    header : dict
        A dictionary container the header keywords from the file being read in
    type_ : str
        A short string describing the type of data being mapped
        
    TODO
    ----
    Rename method to reflect the fact that is now does more than just map
    header values (i.e. cmap, norm)
    
    Returns
    -------
    out: dict
        A new mapped dictionary of useful header values
    """
    date_fmt1 = "%Y-%m-%dT%H:%M:%S.%f"
    
    if type_ == "aia":
        return {
            "cmap": cm.gray,
            "norm": colors.Normalize(5, 1024, True),
            "date": datetime.strptime(header['date-obs'], date_fmt1),
            "det": "AIA",
            "inst": "AIA",
            "meas": header['wavelnth'],
            "obs": "SDO",
            "name": "AIA %s" % header['wavelnth'],
            "r_sun": header['r_sun']
        }
    elif type_ == "hmi":
        return {
            "cmap": None,
            "norm": None,
            "date": datetime.strptime(header['date-obs'], date_fmt1),
            "det": "HMI",
            "inst": "HMI",
            "meas": header['content'].lower(),
            "obs": "SDO",
            "name": "HMI %s" % header['content'].lower(),
            "r_sun": header['r_sun']
        }
