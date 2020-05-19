#
# Testing functions for a mapsequence solar derotation functionality.
#
import os
from copy import deepcopy

import numpy as np
import pytest
from numpy.testing import assert_allclose

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose

import sunpy.data.test
import sunpy.map
from sunpy.physics.solar_rotation import calculate_solar_rotate_shift, mapsequence_solar_derotate


@pytest.fixture
def aia171_test_map():
    testpath = sunpy.data.test.rootdir
    return sunpy.map.Map(os.path.join(testpath, 'aia_171_level1.fits'))


@pytest.fixture
def aia171_test_submap(aia171_test_map):
    return aia171_test_map.submap(SkyCoord(((0, 0), (400, 500))*u.arcsec,
                                           frame=aia171_test_map.coordinate_frame))


@pytest.fixture
def aia171_test_mapsequence(aia171_test_submap):
    m2header = deepcopy(aia171_test_submap.meta)
    m2header['date-obs'] = '2011-02-15T01:00:00.34'
    m2 = sunpy.map.Map((aia171_test_submap.data, m2header))
    m3header = deepcopy(aia171_test_submap.meta)
    m3header['date-obs'] = '2011-02-15T02:00:00.34'
    m3 = sunpy.map.Map((aia171_test_submap.data, m3header))
    return sunpy.map.Map([aia171_test_submap, m2, m3], sequence=True)


# Known displacements for these mapsequence layers when the layer index is set to 0
@pytest.fixture
def known_displacements_layer_index0():
    return {'x': np.asarray([1.136868e-13, -1.048550e+01, -2.099555e+01]),
            'y': np.asarray([1.136868e-13, 2.685361e-01, 5.230608e-01])}


# Known displacements for these mapsequence layers when the layer index is set to 1
@pytest.fixture
def known_displacements_layer_index1():
    return {'x': np.asarray([1.045963e+01, 8.526513e-14, -1.048549e+01]),
            'y': np.asarray([-2.825134e-01, 1.136868e-13, 2.685360e-01])}


def test_calculate_solar_rotate_shift(aia171_test_mapsequence, known_displacements_layer_index0, known_displacements_layer_index1):
    # Test that the default works
    test_output = calculate_solar_rotate_shift(aia171_test_mapsequence)
    assert_allclose(test_output['x'].to('arcsec').value,
                    known_displacements_layer_index0['x'], rtol=5e-2, atol=1e-5)
    assert_allclose(test_output['y'].to('arcsec').value,
                    known_displacements_layer_index0['y'], rtol=5e-2, atol=1e-5)

    # Test that the rotation relative to a nonzero layer_index works
    test_output = calculate_solar_rotate_shift(aia171_test_mapsequence, layer_index=1)
    print(test_output['x'].to('arcsec').value)
    print(test_output['y'].to('arcsec').value)
    assert_allclose(test_output['x'].to('arcsec').value,
                    known_displacements_layer_index1['x'], rtol=5e-2, atol=1e-5)
    assert_allclose(test_output['y'].to('arcsec').value,
                    known_displacements_layer_index1['y'], rtol=5e-2, atol=1e-5)


def test_mapsequence_solar_derotate(aia171_test_mapsequence, aia171_test_submap):
    # Test that a mapsequence is returned when the clipping is False.
    tmc = mapsequence_solar_derotate(aia171_test_mapsequence, clip=False)
    assert(isinstance(tmc, sunpy.map.MapSequence))

    # Test that all entries have the same shape when clipping is False
    for m in tmc:
        assert(m.data.shape == aia171_test_submap.data.shape)

    # Test that a mapsequence is returned on default clipping (clipping is True)
    tmc = mapsequence_solar_derotate(aia171_test_mapsequence)
    assert(isinstance(tmc, sunpy.map.MapSequence))

    # Test that the shape of data is correct when clipped
    clipped_shape = (26, 20)
    for m in tmc:
        assert(m.data.shape == clipped_shape)

    # Test that the returned reference pixels are correctly displaced.
    layer_index = 0
    derotated = mapsequence_solar_derotate(
        aia171_test_mapsequence, clip=True, layer_index=layer_index)
    tshift = calculate_solar_rotate_shift(aia171_test_mapsequence, layer_index=layer_index)
    derotated_reference_pixel_at_layer_index = derotated[layer_index].reference_pixel
    for i, m_derotated in enumerate(derotated):
        for i_s, s in enumerate(['x', 'y']):
            diff_in_rotated_reference_pixel = derotated[i].reference_pixel[i_s] - \
                derotated_reference_pixel_at_layer_index[i_s]
            diff_arcsec = tshift[s][i] - tshift[s][layer_index]
            diff_pixel = diff_arcsec / m.scale[0]
            assert_quantity_allclose(diff_in_rotated_reference_pixel, diff_pixel, rtol=5e-2)
