from __future__ import absolute_import

from sunpy.image.transform import affine_transform
import numpy as np
from skimage import transform as tf
import skimage.data as images
import pytest

# Define test image first so it's accessable to all functions.
original = images.camera().astype('float')

# Tolerance for tests
rtol = 1.0e-15

def compare_results(expect, result, allclose=True):
    """
    Function to check that the obtained results are what was expected, to
    within the relative tolerance defined above.
    """
    # Outermost pixels can contain artefacts which will be ignored.
    exp = expect[1:-1, 1:-1]
    res = result[1:-1, 1:-1]
    t1 = abs(exp.mean() - res.mean()) <= rtol*exp.mean()

    #Don't do the allclose test for scipy as the bicubic algorthm has edge effects
    if allclose:
        t2 = np.allclose(exp, res, rtol=rtol)  #TODO: Develop a better way of testing this
    else:
        t2 = True

    return t1 and t2

@pytest.mark.parametrize("angle, k", [(90.0, 1), (-90.0, -1), (-270.0, 1),
                                      (-90.0, 3), (360.0, 0), (-360.0, 0)])
def test_rotation(angle, k):
    # Test rotation against expected outcome
    angle = np.radians(angle)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, -s], [s, c]])
    expected = np.rot90(original, k=k)

    #Run the tests at order 4 as it produces more accurate 90 deg rotations
    rot = affine_transform(original, order=4, rmatrix=rmatrix)
    assert compare_results(expected, rot)

    # TODO: Check incremental 360 degree rotation against original image

    # Check derotated image against original
    derot_matrix = np.array([[c, s], [-s, c]])
    derot = affine_transform(rot, order=4, rmatrix=derot_matrix)
    assert compare_results(original, derot)

@pytest.mark.parametrize("angle, k", [(90.0, 1), (-90.0, -1), (-270.0, 1),
                                      (-90.0, 3), (360.0, 0), (-360.0, 0)])
def test_scipy_rotation(angle, k):
    # Test rotation against expected outcome
    angle = np.radians(angle)
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, -s], [s, c]])
    expected = np.rot90(original, k=k)
    rot = affine_transform(original, rmatrix=rmatrix, use_scipy=True)
    assert compare_results(expected, rot, allclose=False)

    # TODO: Check incremental 360 degree rotation against original image

    # Check derotated image against original
    derot_matrix = np.array([[c, s], [-s, c]])
    derot = affine_transform(rot, rmatrix=derot_matrix, use_scipy=True)
    assert compare_results(original, derot, allclose=False)

dx_values, dy_values = range(-100, 101, 100)*3, range(-100, 101, 100)*3
dy_values.sort()
@pytest.mark.parametrize("dx, dy", zip(dx_values, dy_values))
def test_shift(dx, dy):
    # Rotation center for all translation tests.
    image_center = np.array(original.shape)/2.0 - 0.5

    # No rotation for all translation tests.
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])

    # Check a shifted shape against expected outcome
    expected = np.roll(np.roll(original, dx, axis=1), dy, axis=0)
    rcen = image_center + np.array([dx, dy])
    shift = affine_transform(original, rmatrix=rmatrix, recenter=True, image_center=rcen)
    ymin, ymax = max([0, dy]), min([original.shape[1], original.shape[1]+dy])
    xmin, xmax = max([0, dx]), min([original.shape[0], original.shape[0]+dx])
    compare_results(expected[ymin:ymax, xmin:xmax], shift[ymin:ymax, xmin:xmax])

    # Check shifted and unshifted shape against original image
    rcen = image_center - np.array([dx, dy])
    unshift = affine_transform(shift, rmatrix=rmatrix, recenter=True, image_center=rcen)
    # Need to ignore the portion of the image cut off by the first shift
    ymin, ymax = max([0, -dy]), min([original.shape[1], original.shape[1]-dy])
    xmin, xmax = max([0, -dx]), min([original.shape[0], original.shape[0]-dx])
    compare_results(original[ymin:ymax, xmin:xmax], unshift[ymin:ymax, xmin:xmax])


@pytest.mark.parametrize("scale_factor", [0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
def test_scale(scale_factor):
    # No rotation for all scaling tests.
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])

    # Check a scaled image against the expected outcome
    newim = tf.rescale(original/original.max(), scale_factor, order=4,
                       mode='constant') * original.max()
    # Old width and new center of image
    w = original.shape[0]/2.0 - 0.5
    new_c = (newim.shape[0]/2.0) - 0.5
    expected = np.zeros(original.shape)
    upper = w+new_c+1
    if scale_factor > 1:
        lower = new_c-w
        expected = newim[lower:upper, lower:upper]
    else:
        lower = w-new_c
        expected[lower:upper, lower:upper] = newim
    scale = affine_transform(original, rmatrix=rmatrix, scale=scale_factor)
    compare_results(expected, scale)


@pytest.mark.parametrize("angle, dx, dy, scale_factor", [(90, -100, 50, 0.25),
                                                          (-90, 50, -100, 0.75),
                                                          (180, 100, 50, 1.5)])
def test_all(angle, dx, dy, scale_factor):
    """
    Tests to make sure that combinations of scaling, shifting and rotation
    produce the expected output.
    """
    k = int(angle/90)
    angle = np.radians(angle)
    image_center = np.array(original.shape)/2.0 - 0.5

    # Check a shifted, rotated and scaled shape against expected outcome
    c = np.cos(angle); s = np.sin(angle)
    rmatrix = np.array([[c, -s], [s, c]])
    scale = tf.rescale(original/original.max(), scale_factor, order=4,
                       mode='constant') * original.max()
    new = np.zeros(original.shape)

    # Old width and new center of image
    w = np.array(original.shape[0])/2.0 - 0.5
    new_c = (np.array(scale.shape[0])/2.0 - 0.5)
    upper = w+new_c+1
    if scale_factor > 1:
        lower = new_c-w
        new = scale[lower:upper, lower:upper]
    else:
        lower = w-new_c
        new[lower:upper, lower:upper] = scale
    disp = np.array([dx, dy])
    rcen = image_center + disp
    rot = np.rot90(new, k=k)
    shift = np.roll(np.roll(rot, dx, axis=1), dy, axis=0)
    expected = shift
    rotscaleshift = affine_transform(original, rmatrix=rmatrix, scale=scale_factor,
                        recenter=True, image_center=rcen)
    w = np.array(expected.shape[0])/2.0 - 0.5
    new_c = (np.array(rotscaleshift.shape[0])/2.0 - 0.5)
    upper = w+new_c+1
    if scale_factor > 1:
        lower = new_c-w
        expected = rotscaleshift[lower:upper, lower:upper]
    else:
        lower = w-new_c
        expected[lower:upper, lower:upper] = rotscaleshift
    compare_results(expected, rotscaleshift)

    # Check a rotated/shifted and restored image against original
    transformed = affine_transform(original, rmatrix=rmatrix, scale=1.0, recenter=True,
                      image_center=rcen)
    rcen = image_center - np.dot(rmatrix, np.array([dx, dy]))
    dx, dy = np.dot(rmatrix, disp)
    rmatrix = np.array([[c, s], [-s, c]])
    inverse = affine_transform(transformed, rmatrix=rmatrix, scale=1.0, recenter=True,
                  image_center=rcen)

    # Need to ignore the portion of the image cut off by the first shift
    # (which isn't the portion you'd expect, because of the rotation)
    ymin, ymax = max([0, -dy]), min([original.shape[1], original.shape[1]-dy])
    xmin, xmax = max([0, -dx]), min([original.shape[0], original.shape[0]-dx])
    compare_results(original[ymin:ymax, xmin:xmax], inverse[ymin:ymax, xmin:xmax])
