# Author: Jack Ireland
#
# Testing functions for a mapcube solar derotation functionality.
#
import os
from copy import deepcopy

import pytest
import numpy as np
from numpy.testing import assert_allclose

import astropy.units as u

import sunpy.data.test
from sunpy import map
from sunpy.physics.solar_rotation import calculate_solar_rotate_shift, mapcube_solar_derotate



@pytest.fixture
def aia171_test_map():
    testpath = sunpy.data.test.rootdir
    return sunpy.map.Map(os.path.join(testpath, 'aia_171_level1.fits'))


@pytest.fixture
def aia171_test_submap(aia171_test_map):
    return aia171_test_map.submap((0, 400)*u.arcsec, (0, 500)*u.arcsec)


@pytest.fixture
def aia171_test_mapcube(aia171_test_submap):
    m2header = deepcopy(aia171_test_submap.meta)
    m2header['date-obs'] = '2011-02-15T01:00:00.34'
    m2 = map.Map((aia171_test_submap.data, m2header))
    m3header = deepcopy(aia171_test_submap.meta)
    m3header['date-obs'] = '2011-02-15T02:00:00.34'
    m3 = map.Map((aia171_test_submap.data, m3header))
    return map.Map([aia171_test_submap, m2, m3], cube=True)


# Known displacements for these mapcube layers when the layer index is set to 0
@pytest.fixture
def known_displacements_layer_index0():
    return {'x': np.asarray([1.222134e-12, -9.142634e+00, -1.830314e+01]),
            'y': np.asarray([1.506351e-12, 1.955316e-01, 3.805245e-01])}


# Known displacements for these mapcube layers when the layer index is set to 1
@pytest.fixture
def known_displacements_layer_index1():
    return {'x': np.asarray([9.123862e+00, 6.821210e-13, -9.142596e+00]),
            'y': np.asarray([-2.061722e-01, 9.663381e-13, 1.956522e-01])}


def test_calculate_solar_rotate_shift(aia171_test_mapcube, known_displacements_layer_index0, known_displacements_layer_index1):
    # Test that the default works
    test_output = calculate_solar_rotate_shift(aia171_test_mapcube)
    assert_allclose(test_output['x'].to('arcsec').value, known_displacements_layer_index0['x'], rtol=5e-2, atol=0)
    assert_allclose(test_output['y'].to('arcsec').value, known_displacements_layer_index0['y'], rtol=5e-2, atol=0)

    # Test that the rotation relative to a nonzero layer_index works
    test_output = calculate_solar_rotate_shift(aia171_test_mapcube, layer_index=1)
    assert_allclose(test_output['x'].to('arcsec').value, known_displacements_layer_index1['x'], rtol=5e-2, atol=0)
    assert_allclose(test_output['y'].to('arcsec').value, known_displacements_layer_index1['y'], rtol=5e-2, atol=0)


def test_mapcube_solar_derotate(aia171_test_mapcube, aia171_test_submap):
    # Test that a mapcube is returned when the clipping is False
    tmc = mapcube_solar_derotate(aia171_test_mapcube, clip=False)
    assert(isinstance(tmc, map.MapCube))

    # Test that all entries have the same shape - nothing clipped
    for m in tmc:
        assert(m.data.shape == aia171_test_submap.data.shape)

    # Test that the returned reference pixels are correctly displaced.
    tmc = mapcube_solar_derotate(aia171_test_mapcube, clip=True)
    tshift = calculate_solar_rotate_shift(aia171_test_mapcube, layer_index=1)
    for im, m in enumerate(tmc):
        for i_s, s in enumerate(['x', 'y']):
            assert_allclose(m.reference_pixel[i_s],
                            aia171_test_submap.reference_pixel[i_s] +
                            tshift[s][im] / m.scale[i_s] -
                            tshift[s][0] / m.scale[i_s],
                            rtol=5e-2, atol=0)

    # Test that a mapcube is returned on default clipping (clipping is True)
    tmc = mapcube_solar_derotate(aia171_test_mapcube)
    assert(isinstance(tmc, map.MapCube))

    # Test that the shape of data is correct when clipped
    clipped_shape = (25, 19)
    for m in tmc:
        assert(m.data.shape == clipped_shape)
