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
rtol = 1.0e-5

def plot_results(expect, result):#, diff):
    """
    Function to plot the results to be shown in the event that the test fails.
    """
    diff = abs(expect - result)

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


def compare_results(expect, result):
    """
    Function to check that the obtained results are what was expected, to
    within the relative tolerance defined above.
    """
    # Outermost pixels can contain artefacts which will be ignored.
    exp = expect[1:-1, 1:-1]
    res = result[1:-1, 1:-1]
    assert abs(exp.mean() - res.mean()) <= rtol*exp.mean()
    assert np.allclose(exp, res, rtol=rtol)


def test_rotation():
    # Rotation center for all rotation tests.
    rotation_center = np.array(original.shape)/2.0 - 0.5
    
    # Test 90 degree rotation against expected outcome
    angle = np.radians(-90.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    expected = np.rot90(original)
    rot = aff(original, rmatrix=rmatrix)#, recenter=True, rotation_center=rotation_center)
    plot_results(expected, rot)
    compare_results(expected, rot)
    plt.close()

    # Test 90 degree rotation against -270 degree rotation
    angle = np.radians(270.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    expected = np.rot90(original)
    rot = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rotation_center)
    plot_results(expected, rot)
    compare_results(expected, rot)#, '90 degree rotation against -270 degree rotation')
    plt.close()

    # Test -90 degree rotation against 270 degree rotation
    angle = np.radians(90.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    expected = np.rot90(original, 3)
    rot = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rotation_center)
    plot_results(expected, rot)
    compare_results(expected, rot)#, '-90 degree rotation against 270 degree rotation')
    plt.close()

    # Check 360 degree rotation against original image
    angle = np.radians(-360.0)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    rot = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rotation_center)
    plot_results(original, rot)
    compare_results(original, rot)#, '360 degree rotation')
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
    plot_results(original, derot)
    compare_results(original, derot)#, 'rotation and derotation')
    plt.close()


def test_shift():
    # Rotation center for all translation tests.
    rotation_center = np.array(original.shape)/2.0 - 0.5
    # No rotation for all translation tests.
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])

    # Check a shifted shape against expected outcome
    expected = np.zeros(original.shape)
    expected[:-20, 100:] = original[20:, :-100]
    rcen = rotation_center + np.array([100, -20])
    shift = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rcen)
    plot_results(expected, shift)
    compare_results(expected, shift)#, 'translation')
    plt.close()

    # Check shifted and unshifted shape against original image
    rcen = rotation_center + np.array([100, -20])
    shift = aff(original, rmatrix=rmatrix, recenter=True, rotation_center=rcen)
    rcen = rotation_center - np.array([100, -20])
    unshift = aff(shift/shift.max(), rmatrix=rmatrix, recenter=True, rotation_center=rcen) * shift.max()
    plot_results(original, unshift)
    # Need to ignore the portion of the image cut off by the first shift
    compare_results(original[20:,:-100], unshift[20:,:-100])#, 'translation and inverse translation')
    plt.close()


def test_scale(scale_factor=0.5):
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
    plot_results(expected, scale)
    compare_results(expected, scale)#, 'scaling')
    plt.close()


def test_all(scale_factor=0.5):
    angle = np.radians(-90.0)
    rotation_center = np.array(original.shape)/2.0 - 0.5
    # Check a shifted, rotated and scaled shape against expected outcome
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    scale = tf.rescale(original, scale_factor, order=4, mode='constant') * original.max()
    new = np.zeros(original.shape)
    # Old width and new centre of image
    w = np.array(original.shape[0])/2.0 - 0.5
    new_c = (np.array(scale.shape[0])/2.0 - 0.5)
    upper = w+new_c+1
    if scale_factor > 1:
        lower = new_c-w
        new = scale[lower:upper, lower:upper]
    else:
        lower = w-new_c
        new[lower:upper, lower:upper] = scale
    disp = np.array([100, -20])
    rcen = rotation_center + disp
    rot = np.rot90(new)
    shift = np.zeros(rot.shape)
    shift[:disp[1], disp[0]:] = rot[-disp[1]:, :-disp[0]]
    expected = shift
    rotscaleshift = aff(original, rmatrix=rmatrix, scale=scale_factor, recenter=True, rotation_center=rcen)
    w = np.array(expected.shape[0])/2.0 - 0.5
    new_c = (np.array(rotscaleshift.shape[0])/2.0 - 0.5)
    upper = w+new_c+1
    if scale_factor > 1:
        lower = new_c-w
        expected = rotscaleshift[lower:upper, lower:upper]
    else:
        lower = w-new_c
        expected[lower:upper, lower:upper] = rotscaleshift
    plot_results(expected, rotscaleshift)
    compare_results(expected, rotscaleshift)#, 'combined rotation, scaling and translation')
    plt.close()

    # Check a rotated/shifted and restored image against original
    angle = np.radians(-90.0)
    rotation_center = np.array(original.shape)/2.0 - 0.5
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    rcen = rotation_center + np.array([20, -100])
    transformed = aff(original, rmatrix=rmatrix, scale=1.0, recenter=True, rotation_center=rcen)
    angle = np.radians(90.0)
    rcen = rotation_center - np.dot(rmatrix, np.array([20, -100]))
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, s], [-s, c]])
    inverse = aff(transformed/transformed.max(), rmatrix=rmatrix, scale=1.0, recenter=True, rotation_center=rcen) * transformed.max()
    plot_results(original, inverse)
    # Need to ignore the portion of the image cut off by the first shift
    # (which isn't the portion you'd expect, because of the rotation)
    compare_results(original[:-20,:-100], inverse[:-20,:-100])#, 'combined rotation, scaling and translation')
    plt.close()


def alltests():
    try:
        test_rotation()
        test_shift()
        test_scale()
        test_all()
    except AssertionError:
        plt.show()
        raise

if __name__ == "__main__":
    alltests()