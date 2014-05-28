from __future__ import absolute_import

import sunpy
from sunpy.instr.aia import aiaprep
from sunpy.image.tests.test_rotate import alltests as test_rotation
import numpy as np
import matplotlib.pyplot as plt
import skimage.data as images
import datetime as dt

# Define test image first so it's accessable to all functions.
data = images.camera()
# An entirely made up test header.
header = {'cdelt1': 0.6,
          'cdelt2': 0.6,
          'crpix1': data.shape[1]/2.0 - 0.5,
          'crpix2': data.shape[0]/2.0 - 0.5,
          'cunit1': 'arcsec',
          'cunit2': 'arcsec',
          'lvl_num': 1.0,
          'naxis': 2,
          'naxis1': data.shape[1],
          'naxis2': data.shape[0],
          'date-obs': sunpy.time.parse_time(dt.datetime.now())}
original = sunpy.map.Map(data, header)


def test_aiaprep():
    # Test that the affine transformation is working first
    test_rotation()
    
    # Test that header info for the map has been correctly updated
    # Check all of these for Map attributes and .meta values?
    # Also may be worth checking they stay the same when saved, I'm sure I've had issues with that before.
    print 'Testing header values ...',
    prep_map = aiaprep(original)
    # Check crpix values
    assert prep_map.meta['crpix1'] == prep_map.shape[1]/2.0 - 0.5
    assert prep_map.meta['crpix2'] == prep_map.shape[1]/2.0 - 0.5
    # Check cdelt values
    assert prep_map.meta['cdelt1'] == 0.6
    assert prep_map.meta['cdelt2'] == 0.6
    # Check crota values
    assert prep_map.meta['crota1'] == 0.0
    assert prep_map.meta['crota2'] == 0.0
    # Check level number
    assert prep_map.meta['lvl_num'] == 1.5
    print 'Passed'


if __name__ == "__main__":
    try:
        test_aiaprep()
    except AssertionError:
        print 'Failed'
        raise
    else:
        print 'All tests completed successfully'