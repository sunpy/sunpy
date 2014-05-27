from __future__ import absolute_import

import sunpy
from sunpy.instr.aia import aiaprep
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

# Tolerance for tests
rtol = 1.0e-6

def plot_results(expect, result, diff):
    fig = plt.figure()
    
    fig.add_subplot(1, 3, 1)
    expect.plot()
    plt.colorbar(orientation='horizontal')
    
    fig.add_subplot(1, 3, 2)
    result.plot()
    plt.colorbar(orientation='horizontal')
    
    fig.add_subplot(1, 3, 3)
    diff.plot()
    plt.colorbar(orientation='horizontal')


def compare_results(expect, result, testmessg):
    exp = expect[1:-1, 1:-1]
    res = result[1:-1, 1:-1]
    #print exp.mean(), res.mean(), abs(exp.mean() - res.mean())
    print 'Testing', testmessg, '...',
    assert abs(exp.mean() - res.mean()) <= rtol*exp.mean()
    assert np.allclose(exp, res, rtol=rtol)
    print 'Passed'


def test_aiaprep_all(scale=0.5):
    print '\n==== Combined tests ===='
    # Check a shifted, rotated and scaled shape against expected outcome
    # Rotate original image
    #rot_data = np.rot90(original.data.copy(), -1)
    # Shift rotated image
    shift_data = np.zeros(original.shape)
    shift_data[:-20, 100:] = original.data.copy()[20:, :-100]
    """# Scale rotated image
    newim = tf.rescale(rot_data.copy(), scale, order=4, mode='constant', cval=original.min()) * rot_data.max()
    w = original.shape[0]/2.0
    new_w = newim.shape[0]/2.0
    scale_data = np.zeros(original.shape)
    scale_data[w-new_w:w+new_w, w-new_w:w+new_w] = newim
    input = sunpy.map.GenericMap(shift_data, original.meta.copy())"""
    input = sunpy.map.GenericMap(shift_data.copy(), original.meta.copy())
    # Adjust header values so aiaprep() will reproduce expect_data
    #input.meta['crota2'] = 90.0
    #input.meta['crpix1'] -= 100#*scale
    #input.meta['crpix2'] += 20#*scale
    input.meta['cdelt1'] = 0.6#*scale
    input.meta['cdelt2'] = 0.6#*scale
    print input.reference_pixel
    prep_map = aiaprep(original)#input)
    print prep_map.reference_pixel
    diff_map = sunpy.map.GenericMap(abs(original.data-prep_map.data), prep_map.meta)
    plot_results(original, prep_map, diff_map)
    compare_results(original.data, prep_map.data, 'combined aiaprep() functionality')
    plt.close()

    """# Check a prepped and de-prepped shape against original image
    original.meta['crota2'] = 10.0
    prep_map = aiaprep(original)
    prep_map.meta['crota2'] = -10.0
    prep_map.meta['crpix1'], original.meta['crpix2'] = data.shape[1]/2.0 - 0.5, data.shape[0]/2.0 - 0.5 # Not entirely sure about this either
    prep_map.meta['cdelt1'], original.meta['cdelt2'] = 1.2, 1.2
    prep_map = aiaprep(prep_map)
    diff_map = sunpy.map.GenericMap(original.data-prep_map.data, prep_map.meta)
    plot_results(original, prep_map, diff_map)
    assert np.allclose(expected.data, prep_map.data, rtol=rtol)
    assert abs(expected.data.mean() - prep_map.data.mean()) <= rtol
    plt.close()"""

    # Test that header info for the map has been correctly updated
    # Check all of these for Map attributes and .meta values?
    # Also may be worth checking they stay the same when saved, I'm sure I've had issues with that before.
    print 'Testing header values ...',
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


try:
    test_aiaprep_rotation()
    test_aiaprep_shift()
    test_aiaprep_scale()
    test_aiaprep_all()
except AssertionError:
    print 'Failed'
    plt.show()
    raise