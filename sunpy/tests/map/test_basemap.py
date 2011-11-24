"""
BaseMap tests
"""
from __future__ import absolute_import

#pylint: disable=C0103,R0904,W0201,W0232
import sunpy
import pyfits
import numpy as np

class TestBaseMap:
    """Tests the BaseMap class"""
    def setup_class(self):
        self.file = sunpy.AIA_171_IMAGE
        self.map = sunpy.Map(self.file)
        self.fits = pyfits.open(self.file)
        self.fits.verify('silentfix')

    def teardown_class(self):
        self.map = None
        self.fits = None
        
    def test_data_to_pixel(self):
        """Make sure conversion from data units to pixels is accurate"""
        # Convert center coords (not FITS starts from 1,1)
        assert self.map.data_to_pixel(0, 'x') == self.map.header['crpix1'] - 1
        assert self.map.data_to_pixel(0, 'y') == self.map.header['crpix2'] - 1
        
        # TODO: also check data value that should map to pixel 1023, 1023
        
    def test_data_range(self):
        """Make sure xrange and yrange work"""
        xmax = self.map.header['cdelt1'] * self.map.header['naxis1'] / 2
        ymax = self.map.header['cdelt2'] * self.map.header['naxis2'] / 2
        
        assert (np.array([-1, 1]) * xmax == self.map.xrange).all()
        assert (np.array([-1, 1]) * ymax == self.map.yrange).all()
        
    def test_submap(self):
        """Check data and header information for a submap"""
        width = self.map.shape[1]
        height = self.map.shape[0]

        # Create a submap of the top-right quadrant of the image
        submap = self.map[height/2:height, width/2:width]
        
        # Expected offset for center
        offset = {
            "x": self.map.header.get('crpix1') - width / 2,
            "y": self.map.header.get('crpix2') - height / 2,
        }
        
        # Check to see if submap header was updated properly
        assert submap.header.get('crpix1') == offset['x'] 
        assert submap.header.get('crpix1') == offset['y']
        assert submap.header.get('naxis1') == width / 2
        assert submap.header.get('naxis2') == height / 2
        
        # Check data
        assert (np.asarray(self.map)[height/2:height, 
                                     width/2:width] == submap).all()
        
    def test_fits_data_comparison(self):
        """Make sure the data is the same in pyfits and SunPy"""
        assert (self.map == self.fits[0].data).all()

    def test_fits_header_comparison(self):
        """Make sure the header is the same in pyfits and SunPy"""

        # Access fits data once to apply scaling-related changes and update
        # header information in fits[0].header
        self.fits[0].data #pylint: disable=W0104
        
        assert dict(self.map.header) == dict(self.fits[0].header)
