from __future__ import absolute_import

from sunpy.image.rotate import affine_transform as aff
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from skimage import transform as tf
import skimage.data as images

# Define test image first so it's accessable to all functions.
original = images.camera()

# Tolerance for tests
rtol = 1.0e-6

def plot_results(expect, result, diff):
    """
    Function to plot the results to be shown in the event that the test fails.
    """
    fig = plt.figure()
    
    fig.add_subplot(1, 3, 1)
    plt.imshow(expect, cmap=cm.gray)
    plt.colorbar(orientation='horizontal')
    plt.title('Expected image')
    
    fig.add_subplot(1, 3, 2)
    plt.imshow(result, cmap=cm.gray)
    plt.colorbar(orientation='horizontal')
    plt.title('Output produced')
    
    fig.add_subplot(1, 3, 3)
    plt.imshow(diff, cmap=cm.gray)
    plt.colorbar(orientation='horizontal')
    plt.title('Difference in images')


def compare_results(expect, result, testmessg):
    """
    Function to check that the obtained results are what was expected, to
    within the relative tolerance defined above.
    """
    # Outermost pixels can contain artefacts which will be ignored.
    exp = expect[1:-1, 1:-1]
    res = result[1:-1, 1:-1]
    print 'Testing', testmessg, '...',
    assert abs(exp.mean() - res.mean()) <= rtol*exp.mean()
    assert np.allclose(exp, res, rtol=rtol)
    print 'Passed'


def test_rotation():    
    print '\n==== Testing rotation ===='
    # Rotation center for all rotation tests.
    rotation_center = np.array(original.shape)/2.0 - 0.5
    
    # Test 90 degree rotation against expected outcome
    angle = np.radians(-90.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    expected = np.rot90(original)
    rot = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rotation_center)
    diff = abs(expected-rot)
    plot_results(expected, rot, diff)
    compare_results(expected, rot, '90 degree rotation')
    plt.close()

    # Test 90 degree rotation against -270 degree rotation
    angle = np.radians(270.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    expected = np.rot90(original)
    rot = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rotation_center)
    diff = abs(expected-rot)
    plot_results(expected, rot, diff)
    compare_results(expected, rot, '90 degree rotation against -270 degree rotation')
    plt.close()

    # Test -90 degree rotation against 270 degree rotation
    angle = np.radians(-90.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    expected = np.rot90(original, -3)
    rot = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rotation_center)
    diff = abs(expected-rot)
    plot_results(expected, rot, diff)
    compare_results(expected, rot, '-90 degree rotation against 270 degree rotation')
    plt.close()

    # Check 360 degree rotation against original image
    angle = np.radians(-360.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    rot = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rotation_center)
    diff = abs(original-rot)
    plot_results(original, rot, diff)
    compare_results(original, rot, '360 degree rotation')
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


def test_shift():
    print '\n==== Testing translation ===='
    # Rotation center for all translation tests.
    rotation_center = np.array(original.shape)/2.0 - 0.5
    # No rotation for all translation tests.
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])

    # Check a shifted shape against expected outcome
    expected = np.zeros(original.shape)
    expected[:-20, 100:] = original[20:, :-100]
    rcen = rotation_center + np.array([20, -100])
    shift = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rcen)
    diff = abs(expected-shift)
    plot_results(expected, shift, diff)
    compare_results(expected, shift, 'translation')
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


def test_scale(scale_factor=0.5):
    print '\n==== Testing scaling ===='
    # Rotation center for all scaling tests.
    rotation_center = np.array(original.shape)/2.0 - 0.5
    # No rotation for all scaling tests.
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])
    
    # Check a scaled image against the expected outcome
    newim = tf.rescale(original, scale_factor, order=4, mode='constant') * original.max()
    # Old width and new centre of image
    w = original.shape[0]/2.0
    new_c = (newim.shape[0]/2.0)
    expected = np.zeros(original.shape)
    if scale_factor > 1:
        expected = newim[new_c-w:new_c+w, new_c-w:new_c+w]
    else:
        expected[w-new_c:w+new_c, w-new_c:w+new_c] = newim
    scale = aff(original, rmatrix=rmatrix, scale=scale_factor, recenter=True, rotation_center=rotation_center)
    diff = abs(expected-scale)
    plot_results(expected, scale, diff)
    compare_results(expected, scale, 'scaling')
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


def test_all(scale_factor=0.5):
    print '\n==== Combined tests ===='
    angle = np.radians(-90.0)
    rotation_center = np.array(original.shape)/2.0 - 0.5
    # Check a shifted, rotated and scaled shape against expected outcome
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    scale = tf.rescale(original, scale_factor, order=4, mode='constant') * original.max()
    new = np.zeros(original.shape)
    # Old width and new centre of image
    w = np.array(scale.shape)
    new_c = (np.array(scale.shape)/2.0 - 0.5)
    im_min = rotation_center - new_c
    if scale_factor > 1:
        new = scale[new_c-w:new_c+w, new_c-w:new_c+w]
    else:
        new[im_min[0]:w[0]+im_min[0], im_min[1]:w[1]+im_min[1]] = scale
    shift = np.zeros(new.shape)
    disp = np.array([20, -100])
    shift[:-disp[0], -disp[1]:] = new[disp[0]:, :disp[1]]
    rcen = rotation_center + disp
    rot = np.rot90(shift)
    expected = rot
    rotscaleshift = aff(original, rmatrix=rmatrix, scale=scale_factor, recenter=True, rotation_center=rcen)
    diff = abs(expected-rotscaleshift)
    plot_results(expected, rotscaleshift, diff)
    compare_results(expected, rotscaleshift, 'combined rotation, scaling and translation')
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


try:
    test_rotation()
    test_shift()
    test_scale()
    test_all()
except AssertionError:
    print 'Failed'
    plt.show()
    raise