"""
Map tests
"""
from __future__ import absolute_import

import sunpy
import sunpy.map
import sunpy.data.test

from astropy.io import fits
import numpy as np

from itertools import izip

import pytest

filepath = sunpy.data.test.rootdir


# Try different dimensions to ensure that the resample method works
# correctly in cases where the dimensions of the original map are
# are exactly divisible by those of the output map as well as the cases
# in which they aren't.
resample_params = [
    ('linear', (100, 200)),
    ('neighbor', (128, 256)),
    ('nearest', (512, 128)),
    ('spline', (200, 200)),
]

class TestGenericMap:
    """Tests the Map class"""
    def setup_class(self):
        self.file = sunpy.AIA_171_IMAGE
        self.map = sunpy.map.Map(self.file)
        self.fits = fits.open(self.file)
        self.fits.verify('silentfix')

        # include full comment
        fits_comment = self.fits[0].header['COMMENT']

        # PyFITS 2.x
        if isinstance(fits_comment[0], basestring):
            comments = [val for val in fits_comment]
        else:
            # PyFITS 3.x
            comments = [card.value for card in fits_comment]
        comment = "".join(comments).strip()

        # touch data to apply scaling up front
        self.fits[0].data

        self.fits[0].header['COMMENT'] = comment

    def teardown_class(self):
        self.map = None
        self.fits = None

    def test_data_to_pixel(self):
        """Make sure conversion from data units to pixels is accurate"""
        # Check conversion of reference pixel
        # Note: FITS pixels starts from 1,1
        assert self.map.data_to_pixel(self.map.meta['crval1'], 'x') == self.map.meta['crpix1'] - 1
        assert self.map.data_to_pixel(self.map.meta['crval2'], 'y') == self.map.meta['crpix2'] - 1

        # Check conversion of map center
        assert self.map.data_to_pixel(self.map.center['x'], 'x') == (self.map.meta['naxis1'] - 1) / 2.
        assert self.map.data_to_pixel(self.map.center['y'], 'y') == (self.map.meta['naxis2'] - 1) / 2.

        # Check conversion of map edges
        # Note: data coords are at pixel centers, so edges are 0.5 pixels wider
        assert self.map.data_to_pixel(self.map.xrange[0], 'x') == 0. - 0.5
        assert self.map.data_to_pixel(self.map.yrange[0], 'y') == 0. - 0.5
        assert self.map.data_to_pixel(self.map.xrange[1], 'x') == (self.map.meta['naxis1'] - 1) + 0.5
        assert self.map.data_to_pixel(self.map.yrange[1], 'y') == (self.map.meta['naxis2'] - 1) + 0.5

    def test_data_range(self):
        """Make sure xrange and yrange work"""
        assert self.map.xrange[1] - self.map.xrange[0] == self.map.meta['cdelt1'] * self.map.meta['naxis1']
        assert self.map.yrange[1] - self.map.yrange[0] == self.map.meta['cdelt2'] * self.map.meta['naxis2']

        assert np.average(self.map.xrange) == self.map.center['x']
        assert np.average(self.map.yrange) == self.map.center['y']

    def test_submap(self):
        """Check data and header information for a submap"""
        width = self.map.shape[1]
        height = self.map.shape[0]

        # Create a submap of the top-right quadrant of the image
        submap = self.map.submap([height/2.,height], [width/2.,width],
                                 units='pixels')

        # Expected offset for center
        offset = {
            "x": self.map.meta['crpix1'] - width / 2.,
            "y": self.map.meta['crpix2'] - height / 2.,
        }

        # Check to see if submap properties were updated properly
        assert submap.reference_pixel['x'] == offset['x']
        assert submap.reference_pixel['y'] == offset['y']
        assert submap.shape[0] == width / 2.
        assert submap.shape[1] == height / 2.

        # Check to see if header was updated
        assert submap.meta['naxis1'] == width / 2.
        assert submap.meta['naxis2'] == height / 2.

        # Check data
        assert (self.map.data[height/2:height,
                                     width/2:width] == submap.data).all()

    def test_fits_data_comparison(self):
        """Make sure the data is the same in pyfits and SunPy"""
        assert (self.map.data == self.fits[0].data).all()

#TODO: What the really?
#    def test_original_header_comparison(self):
#        """Make sure the header is the same in pyfits and SunPy.
#
#        PyFITS makes a number of changes to the data and header when reading
#        it in including applying scaling and removing the comment from the
#        header cards to handle it separately.
#
#        The manipulations in the setup_class method and here attempt to
#        level the playing field some so that the rest of the things that
#        should be the same can be tested.
#
#        // Keith, July 2012
#        """
#
#        # Access fits data once to apply scaling-related changes and update
#        # header information in fits[0].header
#        #self.fits[0].data #pylint: disable=W0104
#
#        fits_header = dict(self.fits[0].header)
#        map_header = self.map.meta
#
#        # Ignore fields modified by PyFITS
#        for key in ['COMMENT', 'BZERO', 'BSCALE', 'BITPIX']:
#            if key in fits_header:
#                del fits_header[key]
#            if key in map_header:
#                del map_header[key]
#
#        # Remove empty field (newline?) that is added when data is accesed for first time
#        if '' in fits_header:
#            fits_header.pop('')
#
#        for k,v in map_header.items():
#            if v != fits_header[k]:
#                print k
#
#        assert map_header == fits_header

    @pytest.mark.parametrize('sample_method,new_dimensions', resample_params)
    def test_resample_dimensions(self, sample_method, new_dimensions):
        """Check that resampled map has expected dimensions."""
        resampled_map = self.map.resample(new_dimensions, method=sample_method)
        assert resampled_map.shape[1] == new_dimensions[0]
        assert resampled_map.shape[0] == new_dimensions[1]

    @pytest.mark.parametrize('sample_method,new_dimensions', resample_params)
    def test_resample_metadata(self, sample_method, new_dimensions):
        """
        Check that the resampled map has correctly adjusted metadata.
        """
        resampled_map = self.map.resample(new_dimensions, method=sample_method)
        assert float(resampled_map.meta['cdelt1']) / self.map.meta['cdelt1'] \
            == float(self.map.shape[1]) / resampled_map.shape[1]
        assert float(resampled_map.meta['cdelt2']) / self.map.meta['cdelt2'] \
            == float(self.map.shape[0]) / resampled_map.shape[0]
        assert resampled_map.meta['crpix1'] == (resampled_map.shape[1] + 1) / 2.
        assert resampled_map.meta['crpix2'] == (resampled_map.shape[0] + 1) / 2.
        assert resampled_map.meta['crval1'] == self.map.center['x']
        assert resampled_map.meta['crval2'] == self.map.center['y']
        for key in self.map.meta:
            if key not in ('cdelt1', 'cdelt2', 'crpix1', 'crpix2',
                           'crval1', 'crval2'):
                assert resampled_map.meta[key] == self.map.meta[key]
            

    def test_superpixel(self):

        dimensions = (2, 2)
        superpixel_map_sum = self.map.superpixel(dimensions)
        assert superpixel_map_sum.shape[0] == self.map.shape[0]/dimensions[1]
        assert superpixel_map_sum.shape[1] == self.map.shape[1]/dimensions[0]
        assert superpixel_map_sum.data[0][0] == self.map.data[0][0] + self.map.data[0][1] + self.map.data[1][0] + self.map.data[1][1]

        dimensions = (2, 2)
        superpixel_map_avg = self.map.superpixel(dimensions, 'average')
        assert superpixel_map_avg.shape[0] == self.map.shape[0]/dimensions[1]
        assert superpixel_map_avg.shape[1] == self.map.shape[1]/dimensions[0]
        assert superpixel_map_avg.data[0][0] == (self.map.data[0][0] + self.map.data[0][1] + self.map.data[1][0] + self.map.data[1][1])/4.0


    def calc_new_matrix(self, angle):
        angle *= -1  # Counter-clockwise rotation
        #Calulate the parameters for the affine_transform
        c = np.cos(np.deg2rad(angle))
        s = np.sin(np.deg2rad(angle))
        return np.matrix([[c, -s], [s, c]])

    def test_rotate(self):

        rotated_map_1 = self.map.rotate(20)
        rotated_map_2 = rotated_map_1.rotate(20)
        assert rotated_map_2.center == rotated_map_1.center == self.map.center
        assert rotated_map_2.shape == rotated_map_1.shape == self.map.shape
        np.testing.assert_allclose(rotated_map_1.rotation_matrix,
                                   np.dot(self.map.rotation_matrix,
                                          self.calc_new_matrix(20).T))
        np.testing.assert_allclose(rotated_map_2.rotation_matrix,
                                   np.dot(self.map.rotation_matrix,
                                          self.calc_new_matrix(40).T))
        # Rotation of a square map by non-integral multiple of 90 degrees cuts off the corners
        # and assigns the value of 0 to corner pixels. This results in reduction
        # of the mean and an increase in standard deviation.
        assert rotated_map_2.mean() < rotated_map_1.mean() < self.map.mean()
        assert rotated_map_2.std() > rotated_map_1.std() > self.map.std()

        rotated_map_3 = self.map.rotate(0, scale=1.5)
        assert rotated_map_3.mean() > self.map.mean()

        # Mean and std should be equal when angle of rotation is integral multiple
        # of 90 degrees for a square map
        rotated_map_4 = self.map.rotate(90, scale=1.5)
        rotated_map_5 = self.map.rotate(180, scale=1.5)
        assert int(rotated_map_3.mean()) == int(rotated_map_4.mean()) == int(rotated_map_5.mean())
        assert int(rotated_map_3.std()) == int(rotated_map_4.std()) == int(rotated_map_5.std())

    def test_rotate_recenter(self):
        # Check recentering
        image_center = np.array((200, 100))
        rotated_map_6 = self.map.rotate(20, image_center=image_center, recenter=True)
        shift = image_center - np.array(self.map.data.shape)/2. + 0.5
        np.testing.assert_allclose(rotated_map_6.reference_pixel.values(),
                                   np.array(self.map.reference_pixel.values()) + shift[::-1])