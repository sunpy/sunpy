import numpy as np
import pytest
import skimage.data as images
from skimage import transform as tf

from sunpy.image.transform import affine_transform
from sunpy.util import SunpyUserWarning

# Tolerance for tests
RTOL = 1.0e-10


@pytest.fixture
def original():
    # Test image
    return images.camera().astype('float')


@pytest.fixture
def identity():
    return np.array([[1, 0], [0, 1]])


def compare_results(expect, result, allclose=True):
    """
    Function to check that the obtained results are what was expected, to
    within the relative tolerance defined above.
    """
    # Outermost pixels can contain artefacts which will be ignored.
    exp = expect[1:-1, 1:-1]
    res = result[1:-1, 1:-1]
    t1 = abs(exp.mean() - res.mean()) <= RTOL*exp.mean()

    # Don't do the allclose test for scipy as the bicubic algorithm has edge effects
    # TODO: Develop a way of testing this for scipy
    if not allclose:
        return t1
    else:
        notclose = ~np.isclose(exp, res, rtol=RTOL)
        t2 = not np.any(notclose)

        # Print out every mismatch
        if not t2:
            mismatches = np.stack([*notclose.nonzero(), exp[notclose], res[notclose]]).T
            for row in mismatches:
                print(f"i={int(row[0]+1)}, j={int(row[1]+1)}: expected={row[2]}, result={row[3]}, "
                      f"adiff={row[2]-row[3]}, rdiff={(row[2]-row[3])/row[2]}")

    return t1 and t2


@pytest.mark.parametrize("angle, k", [(90.0, 1), (-90.0, -1), (-270.0, 1),
                                      (-90.0, 3), (360.0, 0), (-360.0, 0)])
def test_rotation(original, angle, k):
    # Test rotation against expected outcome
    angle = np.radians(angle)
    c = np.round(np.cos(angle))
    s = np.round(np.sin(angle))
    rmatrix = np.array([[c, -s], [s, c]])
    expected = np.rot90(original, k=k)

    error = False

    # Run the tests at order 4 as it produces more accurate 90 deg rotations
    print("Order 4 rotated")
    rot = affine_transform(original, order=4, rmatrix=rmatrix)
    if not compare_results(expected, rot):
        print(rot.flags)
        print("*** Running again ***")
        compare_results(expected, rot)
        error = True

    # TODO: Check incremental 360 degree rotation against original image

    print("Order 4 derotated")
    # Check derotated image against original
    derot_matrix = np.array([[c, s], [-s, c]])
    derot = affine_transform(rot, order=4, rmatrix=derot_matrix)
    if not compare_results(original, derot):
        print(derot.flags)
        print("*** Running again ***")
        compare_results(original, derot)
        error = True

    # Repeat everything for order=3

    print("Order 3 rotated")
    rot = affine_transform(original, order=3, rmatrix=rmatrix)
    if not compare_results(expected, rot):
        print(rot.flags)
        print("*** Running again ***")
        compare_results(expected, rot)
        error = True

    # TODO: Check incremental 360 degree rotation against original image

    # Check derotated image against original
    print("Order 3 derotated")
    derot_matrix = np.array([[c, s], [-s, c]])
    derot = affine_transform(rot, order=3, rmatrix=derot_matrix)
    if not compare_results(original, derot):
        print(derot.flags)
        print("*** Running again ***")
        compare_results(original, derot)
        error = True

    assert not error


@pytest.mark.parametrize("angle, k", [(90.0, 1), (-90.0, -1), (-270.0, 1),
                                      (-90.0, 3), (360.0, 0), (-360.0, 0)])
def test_scipy_rotation(original, angle, k):
    # Test rotation against expected outcome
    angle = np.radians(angle)
    c = np.round(np.cos(angle))
    s = np.round(np.sin(angle))
    rmatrix = np.array([[c, -s], [s, c]])
    expected = np.rot90(original, k=k)
    rot = affine_transform(original, rmatrix=rmatrix, use_scipy=True)
    assert compare_results(expected, rot, allclose=False)

    # TODO: Check incremental 360 degree rotation against original image

    # Check derotated image against original
    derot_matrix = np.array([[c, s], [-s, c]])
    derot = affine_transform(rot, rmatrix=derot_matrix, use_scipy=True)
    assert compare_results(original, derot, allclose=False)


dx_values, dy_values = list(range(-100, 101, 100))*3, list(range(-100, 101, 100))*3
dy_values.sort()
@pytest.mark.parametrize("dx, dy", list(zip(dx_values, dy_values)))
def test_shift(original, dx, dy):
    # Rotation center for all translation tests.
    image_center = np.array(original.shape)/2.0 - 0.5

    # No rotation for all translation tests.
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])

    # Check a shifted shape against expected outcome
    expected = np.roll(np.roll(original, dx, axis=1), dy, axis=0)
    rcen = image_center - np.array([dx, dy])
    shift = affine_transform(original, rmatrix=rmatrix, recenter=True, image_center=rcen)
    ymin, ymax = max([0, dy]), min([original.shape[1], original.shape[1]+dy])
    xmin, xmax = max([0, dx]), min([original.shape[0], original.shape[0]+dx])
    assert compare_results(expected[ymin:ymax, xmin:xmax], shift[ymin:ymax, xmin:xmax])

    # Check shifted and unshifted shape against original image
    rcen = image_center + np.array([dx, dy])
    unshift = affine_transform(shift, rmatrix=rmatrix, recenter=True, image_center=rcen)
    # Need to ignore the portion of the image cut off by the first shift
    ymin, ymax = max([0, -dy]), min([original.shape[1], original.shape[1]-dy])
    xmin, xmax = max([0, -dx]), min([original.shape[0], original.shape[0]-dx])
    assert compare_results(original[ymin:ymax, xmin:xmax], unshift[ymin:ymax, xmin:xmax])


@pytest.mark.parametrize("scale_factor", [0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
def test_scale(original, scale_factor):
    # No rotation for all scaling tests.
    rmatrix = np.array([[1.0, 0.0], [0.0, 1.0]])

    # Check a scaled image against the expected outcome
    newim = tf.rescale(original / original.max(), scale_factor, order=4,
                       mode='constant', multichannel=False, anti_aliasing=False) * original.max()
    # Old width and new center of image
    w = original.shape[0] / 2.0 - 0.5
    new_c = (newim.shape[0] / 2.0) - 0.5
    expected = np.zeros(original.shape)
    upper = int(w + new_c + 1)
    if scale_factor > 1:
        lower = int(new_c - w)
        expected = newim[lower:upper, lower:upper]
    else:
        lower = int(w - new_c)
        expected[lower:upper, lower:upper] = newim
    scale = affine_transform(original, rmatrix=rmatrix, scale=scale_factor, order=4)
    assert compare_results(expected, scale)


@pytest.mark.parametrize("angle, dx, dy, scale_factor", [(90, -100, 40, 0.25),
                                                         (-90, 40, -80, 0.75),
                                                         (180, 20, 50, 1.5)])
def test_all(original, angle, dx, dy, scale_factor):
    """
    Tests to make sure that combinations of scaling, shifting and rotation
    produce the expected output.
    """
    k = int(angle / 90)
    angle = np.radians(angle)
    image_center = np.array(original.shape) / 2.0 - 0.5

    # Check a shifted, rotated and scaled shape against expected outcome
    c = np.round(np.cos(angle))
    s = np.round(np.sin(angle))
    rmatrix = np.array([[c, -s], [s, c]])
    scale = tf.rescale(original / original.max(), scale_factor, order=4,
                       mode='constant', multichannel=False, anti_aliasing=False) * original.max()
    new = np.zeros(original.shape)

    disp = np.array([dx, dy])
    dxs, dys = np.asarray(disp * scale_factor, dtype=int)
    # Old width and new center of image
    w = np.array(original.shape[0])/2.0 - 0.5
    new_c = (np.array(scale.shape[0])/2.0 - 0.5)
    upper = int(w+new_c+1)
    if scale_factor > 1:
        lower = int(new_c-w)
        new = scale[lower-dys:upper-dys, lower-dxs:upper-dxs]
    else:
        lower = int(w-new_c)
        new[lower+dys:upper+dys, lower+dxs:upper+dxs] = scale
    rcen = image_center - disp
    expected = np.rot90(new, k=k)

    rotscaleshift = affine_transform(original, rmatrix=rmatrix, scale=scale_factor, order=4,
                                     recenter=True, image_center=rcen)
    assert compare_results(expected, rotscaleshift)

    # Check a rotated/shifted and restored image against original
    transformed = affine_transform(original, rmatrix=rmatrix, scale=1.0, order=4, recenter=True,
                                   image_center=rcen)
    inv_rcen = image_center + np.dot(rmatrix.T, np.array([dx, dy]))
    inverse = affine_transform(transformed, rmatrix=rmatrix.T, scale=1.0, order=4, recenter=True,
                               image_center=inv_rcen)

    # Need to ignore the portion of the image cut off by the first shift
    ymin, ymax = max([0, -dy]), min([original.shape[1], original.shape[1]-dy])
    xmin, xmax = max([0, -dx]), min([original.shape[0], original.shape[0]-dx])
    assert compare_results(original[ymin:ymax, xmin:xmax], inverse[ymin:ymax, xmin:xmax])


def test_flat(identity):
    # Test that a flat array can be rotated using scikit-image
    in_arr = np.array([[100]], dtype=np.float64)
    out_arr = affine_transform(in_arr, rmatrix=identity)
    assert np.allclose(in_arr, out_arr, rtol=RTOL)


def test_nan_skimage_low(identity):
    # Test non-replacement of NaN values for scikit-image rotation with order <= 3
    in_arr = np.array([[np.nan]])
    out_arr = affine_transform(in_arr, rmatrix=identity, order=3)
    assert np.all(np.isnan(out_arr))


def test_nan_skimage_high(identity):
    # Test replacement of NaN values for scikit-image rotation with order >=4
    in_arr = np.array([[np.nan]])
    with pytest.warns(SunpyUserWarning, match='Setting NaNs to 0 for higher-order scikit-image rotation.'):
        out_arr = affine_transform(in_arr, rmatrix=identity, order=4)
    assert not np.all(np.isnan(out_arr))


def test_nan_scipy(identity):
    # Test replacement of NaN values for scipy rotation
    in_arr = np.array([[np.nan]])
    with pytest.warns(SunpyUserWarning, match='Setting NaNs to 0 for SciPy rotation.'):
        out_arr = affine_transform(in_arr, rmatrix=identity, use_scipy=True)
    assert not np.all(np.isnan(out_arr))


def test_int(identity):
    # Test casting of integer array to float array
    in_arr = np.array([[100]], dtype=int)
    with pytest.warns(SunpyUserWarning, match='Input data has been cast to float64.'):
        out_arr = affine_transform(in_arr, rmatrix=identity)
    assert np.issubdtype(out_arr.dtype, np.floating)


def test_lots_of_rotations(original):
    # Test a 90-degree rotation many times
    expected = np.rot90(original, k=1)
    matrix = np.array([[0, -1, 511], [1, 0, 0], [0, 0, 1]])
    match = np.zeros(100, bool)
    for i in range(len(match)):
        result = tf.warp(original, matrix, order=4)
        match[i] = compare_results(expected, result)
    assert np.sum(~match) == 0


def test_broken_apart(original):
    from skimage.transform import ProjectiveTransform, warp_coords
    from scipy.ndimage import map_coordinates
    matrix = np.array([[0, -1, 511], [1, 0, 0], [0, 0, 1]])
    coord_map = ProjectiveTransform(matrix=matrix)
    expected_coords = warp_coords(coord_map, original.shape)
    expected_warped = map_coordinates(original, expected_coords, order=4)

    print("Running warp_coords() 100 times")

    match1 = np.zeros(100, bool)
    for i in range(len(match1)):
        result_coords = warp_coords(coord_map, original.shape)
        match1[i] = compare_results(expected_coords[0, : :], result_coords[0, : :])
        match1[i] &= compare_results(expected_coords[1, : :], result_coords[1, : :])

    print(np.sum(~match1))

    print("Running map_coordinates() 100 times")

    match2 = np.zeros(100, bool)
    for i in range(len(match2)):
        result_warped = map_coordinates(original, expected_coords, order=4)
        match2[i] = compare_results(expected_warped, result_warped)

    print(np.sum(~match2))

    print("Running ProjectiveTransform 100 times")

    match3 = np.zeros(100, bool)
    base = np.indices(original.shape, dtype=np.float64).reshape(2, -1).T
    expected = coord_map(base)
    for i in range(len(match3)):
        result = coord_map(base)
        match3[i] = np.allclose(expected, result)

    print(np.sum(~match3))

    print("Running matrix multiplication 100 times")

    match4 = np.zeros(100, bool)
    x, y = np.transpose(base)
    src = np.vstack((x, y, np.ones_like(x)))
    expected = src.T @ matrix.T
    for i in range(len(match4)):
        result = src.T @ matrix.T
        match4[i] = np.allclose(expected, result)

    print(np.sum(~match4))

    assert (np.sum(~match1), np.sum(~match2), np.sum(~match3), np.sum(~match4)) == (0, 0, 0, 0)


def test_minimal_example():
    base = np.indices((512, 512), dtype=np.float64).reshape(2, -1).T
    x, y = np.transpose(base)
    src = np.vstack((x, y, np.ones_like(x)))

    matrix = np.array([[0, -1, 511], [1, 0, 0], [0, 0, 1]])
    expected = src.T @ matrix.T

    mismatches = np.zeros(100, int)
    for i in range(len(mismatches)):
        result = src.T @ matrix.T
        mismatches[i] = (~np.isclose(result, expected)).sum()
        if mismatches[i] != 0:
            print(f"{mismatches[i]} mismatching elements in this multiplication")

    assert np.sum(mismatches != 0) == 0
