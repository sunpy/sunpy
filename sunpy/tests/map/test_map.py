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
        
        # include full comment
        fits_comment = self.fits[0].header.get_comment()
        
        # PyFITS 2.x
        if isinstance(fits_comment[0], basestring):
            comments = [val for val in fits_comment]       
        else:
            # PyFITS 3.x
            comments = [card.value for card in fits_comment]
        comment = "".join(comments).strip()
        
        # touch data to apply scaling up front
        self.fits[0].data        
        
        self.fits[0].header.update('COMMENT', comment)

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
    
    def test_sub(self):
        map_ = sunpy.Map(np.array([[1, 1], [2, 2]], dtype=np.uint8), {})
        minus = map_ - map_
        assert minus.dtype == np.int16
        assert (minus == 0).all()
        assert np.array_equal(
            map_ - 2 * map_,
            np.array([[-1, -1], [-2, -2]], dtype=np.int16)
        )
    
    def test_original_header_comparison(self):
        """Make sure the header is the same in pyfits and SunPy.
        
        PyFITS makes a number of changes to the data and header when reading
        it in including applying scaling and removing the comment from the
        header cards to handle it separately.
        
        The manipulations in the setup_class method and here attempt to
        level the playing field some so that the rest of the things that
        should be the same can be tested.
        
        // Keith, July 2012
        """

        # Access fits data once to apply scaling-related changes and update
        # header information in fits[0].header
        #self.fits[0].data #pylint: disable=W0104

        fits_header = dict(self.fits[0].header)
        map_header = dict(self.map._original_header)
        
        # Ignore fields modified by PyFITS
        for key in ['COMMENT', 'BZERO', 'BSCALE', 'BITPIX']:
            if key in fits_header:
                del fits_header[key]
            if key in map_header:
                del map_header[key]
                
        # Remove empty field (newline?) that is added when data is accesed for first time
        if '' in fits_header:
            fits_header.pop('')
        
        for k,v in map_header.items():
            if v != fits_header[k]:
                print k
        
        assert map_header == fits_header

    # TODO: Add tests for other resample methods (neighbour, nearest, spline)
    def test_linear_resample_dimensions(self):
        """Check that resampled map has expected dimensions."""

        new_dimensions = (100, 200)
        resampled_map = self.map.resample(new_dimensions)
        assert resampled_map.shape[1] == new_dimensions[0]
        assert resampled_map.shape[0] == new_dimensions[1]
