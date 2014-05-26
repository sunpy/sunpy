from __future__ import absolute_import

import sunpy
from sunpy.instr.aia import aiaprep
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
from skimage import transform as tf
import skimage.data as images
import datetime as dt

# Define test image first so it's accessable to all functions.
data = images.camera()#checkerboard()
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
"""
fig = plt.figure()
original.plot(cmap=cm.binary)
plt.show()

header = original.meta.keys()
header.sort()
for item in header:
    print item, original.meta[item]
"""

# Tolerance for tests - high for now just so they run. I'll make it much smaller
# when actually testing
rtol = 1.0e-5#0.001

def plot_results(expect, result, diff):
    fig = plt.figure()
    
    fig.add_subplot(1, 3, 1)
    expect.plot()#cmap=cm.binary)
    plt.colorbar(orientation='horizontal')
    
    fig.add_subplot(1, 3, 2)
    result.plot()#cmap=cm.binary)
    plt.colorbar(orientation='horizontal')
    
    fig.add_subplot(1, 3, 3)
    diff.plot()#cmap=cm.binary)
    plt.colorbar(orientation='horizontal')


def compare_results(expect, result, testmessg):
    exp = expect[1:-1, 1:-1]
    res = result[1:-1, 1:-1]
    #print exp.mean(), res.mean(), abs(exp.mean() - res.mean())
    print 'Testing', testmessg, '...',
    assert abs(exp.mean() - res.mean()) <= rtol*exp.mean()
    assert np.allclose(exp, res, rtol=rtol)
    print 'Passed'


def test_aiaprep_rotation():    
    print '\n==== Testing rotation ===='
    # Test 90 degree rotation against expected outcome
    original.meta['crota2'] = 90.0
    expected = sunpy.map.GenericMap(np.rot90(original.data.copy()), original.meta.copy())
    rot_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(abs(expected.data-rot_map.data), rot_map.meta)
    plot_results(expected, rot_map, diff_map)
    compare_results(expected.data, rot_map.data, '90 degree rotation')
    plt.close()

    # Test 90 degree rotation against -270 degree rotation
    original.meta['crota2'] = -270.0
    expected = sunpy.map.GenericMap(np.rot90(original.data.copy()), original.meta.copy())
    rot_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(abs(expected.data-rot_map.data), rot_map.meta)
    plot_results(expected, rot_map, diff_map)
    compare_results(expected.data, rot_map.data, '90 degree rotation against -270 degree rotation')
    plt.close()

    # Test -90 degree rotation against 270 degree rotation
    original.meta['crota2'] = 90.0
    expected = sunpy.map.GenericMap(np.rot90(original.data.copy(), -3), original.meta.copy())
    rot_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(abs(expected.data-rot_map.data), rot_map.meta)
    plot_results(expected, rot_map, diff_map)
    compare_results(expected.data, rot_map.data, '-90 degree rotation against 270 degree rotation')
    plt.close()

    # Check 360 degree rotation against original image
    original.meta['crota2'] = 360.0
    rot_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(abs(original.data-rot_map.data), rot_map.meta)
    plot_results(original, rot_map, diff_map)
    compare_results(original.data, rot_map.data, '360 degree rotation')
    plt.close()
    
    # Check incremental 360 degree rotation against original image

    # This doesn't work but I have no idea why just now.
    """# Check rotated and derotated image against original
    original.meta['crota2'] = 10.0
    rot_map = aiaprep(original)
    rot_map.meta['crota2'] = -10.0
    rot_map = aiaprep(rot_map)
    diff_map = sunpy.map.GenericMap(abs(original.data-rot_map.data), rot_map.meta)
    plot_results(original, rot_map, diff_map)
    compare_results(original.data, rot_map.data, 'rotation and derotation')
    plt.close()"""


def test_aiaprep_shift():
    print '\n==== Testing translation ===='
    # Check a shifted shape against expected outcome
    original.meta['crota2'] = 0.0
    expect_data = np.zeros(original.shape)
    expect_data[:-20, 100:] = original.data[20:, :-100]
    expected = sunpy.map.GenericMap(expect_data, original.meta.copy())
    original.meta['crpix1'] -= 100
    original.meta['crpix2'] += 20
    shift_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(abs(expected.data-shift_map.data), shift_map.meta)
    plot_results(expected, shift_map, diff_map)
    compare_results(expected.data, shift_map.data, 'translation')
    plt.close()

    # This doesn't work either
    """# Check shifted and unshifted shape against original image
    shift_map = aiaprep(original)
    shift_map.meta['crpix1'], shift_map.meta['crpix2'] = data.shape[1]/2.0 - 0.5, data.shape[0]/2.0 - 0.5
    shift_map = aiaprep(shift_map)
    diff_map = sunpy.map.GenericMap(abs(original.data-shift_map.data), shift_map.meta)
    plot_results(original, shift_map, diff_map)
    compare_results(original.data, shift_map.data, 'translation and detranslation')
    plt.close()"""


def test_aiaprep_scale(scale=0.5):
    print '\n==== Testing scaling ===='
    # Check a scaled image against the expected outcome
    original.meta['crota2'] = 0.0
    original.meta['crpix1'] = data.shape[1]/2.0 - 0.5
    original.meta['crpix2'] = data.shape[0]/2.0 - 0.5
    expect_data = tf.rescale(original.data.copy(), scale, order=4, mode='constant', cval=original.min()) * original.max()
    w = original.shape[0]/2.0
    new_w = expect_data.shape[0]/2.0
    newim = expect_data.copy()
    expect_data = np.zeros(original.shape)
    expect_data[w-new_w:w+new_w, w-new_w:w+new_w] = newim
    expected = sunpy.map.GenericMap(expect_data, original.meta.copy())
    original.meta['cdelt1'], original.meta['cdelt2'] = 0.6*scale, 0.6*scale
    scale_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(abs(expected.data-scale_map.data), scale_map.meta)
    plot_results(expected, scale_map, diff_map)
    compare_results(expected.data, scale_map.data, 'scaling')
    plt.close()
    
    # I doubt this will work based on the other similar tests
    """# Check a scaled and descaled image against the original
    scale_map = aiaprep(original)
    scale_map.meta['cdelt1'], scale_map.meta['cdelt2'] = 1.2, 1.2 # I think this should be the same as changing the target scale in aiaprep()
    scale_map = aiaprep(scale_map)
    diff_map = sunpy.map.GenericMap(abs(original.data-scale_map.data), scale_map.meta)
    plot_results(original, scale_map, diff_map)
    compare_results(original.data, scale_map.data, 'scaling and descaling')
    plt.close()"""


def test_aiaprep_all(scale=0.5):
    print '\n==== Combined tests ===='
    # Check a shifted, rotated and scaled shape against expected outcome
    # Rotate original image
    rot_data = np.rot90(original.data.copy())
    # Scale rotated image
    newim = tf.rescale(rot_data.copy(), scale, order=4, mode='constant', cval=original.min()) * original.max()
    w = original.shape[0]/2.0
    new_w = newim.shape[0]/2.0
    scale_data = np.zeros(original.shape)
    scale_data[w-new_w:w+new_w, w-new_w:w+new_w] = newim
    # Shift rotated image
    shift_data = np.zeros(original.shape)
    shift_data[:-20, 100:] = scale_data[20:, :-100]
    expected = sunpy.map.GenericMap(shift_data, original.meta.copy())
    # Adjust header values so aiaprep() will reproduce expect_data
    original.meta['crota2'] = 90.0
    original.meta['crpix1'] -= 100*scale
    original.meta['crpix2'] += 20*scale
    original.meta['cdelt1'] = 0.6*scale
    original.meta['cdelt2'] = 0.6*scale
    prep_map = aiaprep(original)
    diff_map = sunpy.map.GenericMap(abs(expected.data-prep_map.data), prep_map.meta)
    plot_results(expected, prep_map, diff_map)
    compare_results(expected.data, prep_map.data, 'combined aiaprep() things')
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


try:
    #test_aiaprep_rotation()
    #test_aiaprep_shift()
    #test_aiaprep_scale()
    test_aiaprep_all()
except AssertionError:
    print 'Failed'
    plt.show()
    raise