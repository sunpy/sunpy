# Author: Jack Ireland
#
# Testing functions for a mapcube solar derotation functionality.
#
import numpy as np
from astropy import units as u
from numpy.testing import assert_allclose
from sunpy import AIA_171_IMAGE
from sunpy import map
from sunpy.physics.transforms.solar_rotation import calculate_solar_rotate_shift, mapcube_solar_derotate
from copy import deepcopy

# fake a mapcube
m = map.Map(AIA_171_IMAGE)
m1 = m.submap((0, 400), (0, 500))
m2header = deepcopy(m1.meta)
m2header['date-obs'] = '2011-03-19T11:54:00.34'
m2 = map.Map((m1.data, m2header))
m3header = deepcopy(m1.meta)
m3header['date-obs'] = '2011-03-19T12:54:00.34'
m3 = map.Map((m1.data, m3header))

mc = map.Map([m1, m2, m3], cube=True)

# Known displacements for these mapcube layers when the layer index is set to 0
known_displacements_layer_index0 = {'x': np.asarray([-4.51905180e-12, -9.06805914e+00, -1.81541844e+01]),
                                    'y': np.asarray([-5.57065505e-12, 2.56807726e-01, 5.02761067e-01])}

# Known displacements for these mapcube layers when the layer index is set to 1
known_displacements_layer_index1 = {'x': np.asarray([9.04899290e+00, 2.67164069e-12, -9.06791780e+00]),
                                    'y': np.asarray([-2.67659143e-01, 3.21165317e-12, 2.56829202e-01])}


def test_calculate_solar_rotate_shift():
    # Test that the default works
    test_output = calculate_solar_rotate_shift(mc)
    assert_allclose(test_output['x'].to('arcsec'), known_displacements_layer_index0['x'] * u.arcsec, rtol=5e-2, atol=0)
    assert_allclose(test_output['y'].to('arcsec'), known_displacements_layer_index0['y'] * u.arcsec, rtol=5e-2, atol=0)

    # Test that the rotation relative to a nonzero layer_index works
    test_output = calculate_solar_rotate_shift(mc, layer_index=1)
    assert_allclose(test_output['x'].to('arcsec'), known_displacements_layer_index1['x'] * u.arcsec, rtol=5e-2, atol=0)
    assert_allclose(test_output['y'].to('arcsec'), known_displacements_layer_index1['y'] * u.arcsec, rtol=5e-2, atol=0)


def test_mapcube_solar_derotate():
    # Test that a mapcube is returned when the clipping is False
    tmc = mapcube_solar_derotate(mc, clip=False)
    assert(isinstance(tmc, map.MapCube))

    # Test that all entries have the same shape - nothing clipped
    for m in tmc:
        assert(m.data.shape == m1.data.shape)

    # Test that the returned centers are correctly displaced.
    tshift = calculate_solar_rotate_shift(mc)
    for im, m in enumerate(tmc):
        for s in ['x', 'y']:
            assert_allclose(m.center[s], m1.center[s] -
                            tshift[s][im].to('arcsec').value, rtol=5e-2, atol=0)

    # Test that a mapcube is returned on default clipping (clipping is True)
    tmc = mapcube_solar_derotate(mc)
    assert(isinstance(tmc, map.MapCube))

    # Test that the shape of data is correct when clipped
    clipped_shape = (206, 159)
    for m in tmc:
        assert(m.data.shape == clipped_shape)






