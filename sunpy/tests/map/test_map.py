"""
Map tests
"""
from __future__ import absolute_import

#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103
import sunpy
import pyfits
import numpy as np

class TestMap:
    """Tests the Map class"""
    def setup_class(self):
        self.file = sunpy.AIA_171_IMAGE
        self.map = sunpy.make_map(self.file)
        self.fits = pyfits.open(self.file)
        self.fits.verify('silentfix')

    def teardown_class(self):
        self.map = None
        self.fits = None
        
    def test_data_to_pixel(self):
        """Make sure conversion from data units to pixels is accurate"""
        # Check conversion of reference pixel
        # Note: FITS pixels starts from 1,1
        assert self.map.data_to_pixel(self.map._original_header['crval1'], 'x') == self.map._original_header['crpix1'] - 1
        assert self.map.data_to_pixel(self.map._original_header['crval2'], 'y') == self.map._original_header['crpix2'] - 1
        
        # Check conversion of map center
        assert self.map.data_to_pixel(self.map.center['x'], 'x') == (self.map._original_header['naxis1'] - 1) / 2.
        assert self.map.data_to_pixel(self.map.center['y'], 'y') == (self.map._original_header['naxis2'] - 1) / 2.
        
        # Check conversion of map edges
        # Note: data coords are at pixel centers, so edges are 0.5 pixels wider
        assert self.map.data_to_pixel(self.map.xrange[0], 'x') == 0. - 0.5
        assert self.map.data_to_pixel(self.map.yrange[0], 'y') == 0. - 0.5
        assert self.map.data_to_pixel(self.map.xrange[1], 'x') == (self.map._original_header['naxis1'] - 1) + 0.5
        assert self.map.data_to_pixel(self.map.yrange[1], 'y') == (self.map._original_header['naxis2'] - 1) + 0.5
    
    def test_data_range(self):
        """Make sure xrange and yrange work"""
        assert self.map.xrange[1] - self.map.xrange[0] == self.map._original_header['cdelt1'] * self.map._original_header['naxis1']
        assert self.map.yrange[1] - self.map.yrange[0] == self.map._original_header['cdelt2'] * self.map._original_header['naxis2']
        
        assert np.average(self.map.xrange) == self.map.center['x']
        assert np.average(self.map.yrange) == self.map.center['y']
        
    def test_submap(self):
        """Check data and header information for a submap"""
        width = self.map.shape[1]
        height = self.map.shape[0]

        # Create a submap of the top-right quadrant of the image
        submap = self.map[height/2:height, width/2:width]
        
        # Expected offset for center
        offset = {
            "x": self.map._original_header.get('crpix1') - width / 2,
            "y": self.map._original_header.get('crpix2') - height / 2,
        }
        
        # Check to see if submap properties were updated properly
        assert submap.reference_pixel['x'] == offset['x'] 
        assert submap.reference_pixel['y'] == offset['y']
        assert submap.shape[0] == width / 2
        assert submap.shape[1] == height / 2
        
        # Check to see if header was updated
        submap_header = submap.get_header()
        assert submap_header.get('naxis1') == width / 2
        assert submap_header.get('naxis2') == height / 2
        
        # Check data
        assert (np.asarray(self.map)[height/2:height, 
                                     width/2:width] == submap).all()
        
    def test_fits_data_comparison(self):
        """Make sure the data is the same in pyfits and SunPy"""
        assert (self.map == self.fits[0].data).all()

    def test__original_header_comparison(self):
        """Make sure the header is the same in pyfits and SunPy"""

        # Access fits data once to apply scaling-related changes and update
        # header information in fits[0].header
        self.fits[0].data #pylint: disable=W0104
        
        assert dict(self.map._original_header) == dict(self.fits[0].header)
