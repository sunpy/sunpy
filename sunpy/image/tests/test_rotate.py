from __future__ import absolute_import

from sunpy.image.rotate import affine_transform as aff
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import patches
from skimage import transform as tf
import skimage.data as images

# Define test image first so it's accessable to all functions.
original = images.camera()

# Tolerance for tests
rtol = 1.0e-5

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
    print 'Testing', testmessg, '...'#,
    print exp.mean(), res.mean()
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
    
    # TODO: Check incremental 360 degree rotation against original image

    # Check rotated and derotated image against original
    angle = np.radians(-90.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    rot = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rotation_center, missing=original.mean()/original.max())
    angle = np.radians(90.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    derot = aff(rot/rot.max(), rmatrix=rmatrix, recenter=True, rotation_center=rotation_center) * rot.max()
    diff = abs(original-derot)
    plot_results(original, derot, diff)
    compare_results(original, derot, 'rotation and derotation')
    plt.close()


def test_shift():
    print '\n==== Testing translation ===='
    # Rotation center for all translation tests.
    rotation_center = np.array(original.shape)/2.0 - 0.5
    # No rotation for all translation tests.
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])

    # Check a shifted shape against expected outcome
    expected = np.zeros(original.shape)
    expected[:-20, 100:] = original[20:, :-100]
    rcen = rotation_center + np.array([100, -20])
    shift = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rcen)
    diff = abs(expected-shift)
    plot_results(expected, shift, diff)
    compare_results(expected, shift, 'translation')
    plt.close()

    # Check shifted and unshifted shape against original image
    rcen = rotation_center + np.array([100, -20])
    shift = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rcen)
    rcen = rotation_center - np.array([100, -20])
    unshift = aff(shift/shift.max(), rmatrix=rmatrix, recenter=True, rotation_center=rcen) * shift.max()
    diff = abs(original-unshift)
    plot_results(original, unshift, diff)
    # Need to ignore the portion of the image cut off by the first shift
    compare_results(original[20:,:-100], unshift[20:,:-100], 'translation and inverse translation')
    plt.close()


def test_scale(scale_factor=0.5):
    print '\n==== Testing scaling ===='
    # Rotation center for all scaling tests.
    rotation_center = np.array(original.shape)/2.0 - 0.5
    # No rotation for all scaling tests.
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])
    
    # Check a scaled image against the expected outcome
    newim = tf.rescale(original, scale_factor, order=4, mode='constant') * original.max()
    # Old width and new centre of image
    w = original.shape[0]/2.0 - 0.5
    new_c = (newim.shape[0]/2.0) -0.5
    expected = np.zeros(original.shape)
    upper = w+new_c+1
    if scale_factor > 1:
        lower = new_c-w
        expected = newim[lower:upper, lower:upper]
    else:
        lower = w-new_c
        expected[lower:upper, lower:upper] = newim
    scale = aff(original, rmatrix=rmatrix, scale=scale_factor, recenter=True, rotation_center=rotation_center)
    diff = abs(expected-scale)
    plot_results(expected, scale, diff)
    compare_results(expected, scale, 'scaling')
    plt.close()
    
    # Check a scaled and descaled image against the original
    print scale_factor, 1.0/scale_factor
    scale = aff(original, rmatrix=rmatrix, scale=scale_factor, recenter=True, rotation_center=rotation_center)
    descale = aff(scale/scale.max(), rmatrix=rmatrix, scale=1.0/scale_factor, recenter=True, rotation_center=rotation_center) * scale.max()
    diff = abs(original-descale)
    if scale_factor > 1:
        # Need to ignore the parts of the image cropped by the initial scaling
        scl_w = (original.shape[0]/2.0 - 0.5)/scale_factor
        lower = w-scl_w
        upper = w+scl_w+1
        plot_results(original[lower:upper, lower:upper], descale[lower:upper, lower:upper], diff[lower:upper, lower:upper])
        compare_results(original[lower:upper, lower:upper], descale[lower:upper, lower:upper], 'scaling and inverse scaling')
    else:
        plot_results(original, descale, diff)
        compare_results(original, descale, 'scaling and inverse scaling')
    plt.close()


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
    #expected[:-20, 100:] = original[20:, :-100]
    #rcen = rotation_center + np.array([100, -20])
    disp = np.dot(rmatrix, np.array([100, -20]))
    shift[:disp[1], disp[0]:] = new[-disp[1]:, :-disp[0]]
    rcen = rotation_center + disp
    rot = np.rot90(shift)
    expected = rot
    rotscaleshift = aff(original, rmatrix=rmatrix, scale=scale_factor, recenter=True, rotation_center=rcen)
    diff = abs(expected-rotscaleshift)
    plot_results(expected, rotscaleshift, diff)
    compare_results(expected, rotscaleshift, 'combined rotation, scaling and translation')
    plt.close()

    # Check a prepped and de-prepped shape against original image
    angle = np.radians(-20.0)
    rotation_center = np.array(original.shape)/2.0 - 0.5
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    rcen = rotation_center + np.array([20, -100])
    transformed = aff(original, rmatrix=rmatrix, scale=scale_factor, recenter=True, rotation_center=rcen)
    angle = np.radians(20.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    inverse = aff(transformed, rmatrix=rmatrix, scale=scale_factor, recenter=True, rotation_center=rotation_center)
    diff = abs(original-inverse)
    plot_results(original, inverse, diff)
    compare_results(original, inverse, 'combined rotation, scaling and translation')
    plt.close()


def alltests():
    try:
        test_rotation()
        test_shift()
        test_scale(1.0)
        test_all(1.0)
    except AssertionError:
        print 'Failed'
        plt.show()
        raise

if __name__ == "__main__":
    alltests()