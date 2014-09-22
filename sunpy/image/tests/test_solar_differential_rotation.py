# Author: Jack Ireland
#
# Testing functions for a mapcube solar derotation functionality.
#
import numpy as np
from astropy import units as u
from numpy.testing import assert_allclose
from scipy.ndimage.interpolation import shift
from sunpy import AIA_171_IMAGE
from sunpy import map
from sunpy.image.solar_differential_rotation import mapcube_solar_derotate
from copy import deepcopy

def test_mapcube_solar_derotate():
    # fake a mapcube
    m = map.Map(AIA_171_IMAGE)
    m1 = m.submap((0, 400), (0, 500))
    m2header = deepcopy(m1.meta)
    m2header['date-obs'] = '2011-03-19T11:54:00.34'
    m2 = map.Map((m1.data, m2header))
    mc = map.MapCube([m1, m2])

    # Known displacements for these mapcube layers
    known_displacements = {'x': np.asarray([ -4.51905180e-12,  -9.06805914e+00]),
                           'y': np.asarray([ -5.57065505e-12,   2.56807726e-01])}

    # Test that a mapcube is returned
    test_output = mapcube_solar_derotate(mc)
    # Assert
    assert(isinstance(test_output, map.MapCube))

    # Test the return of the only the displacements.
    test_output = mapcube_solar_derotate(mc, return_displacements_only=True)
    # Assert
    assert(isinstance(test_output, dict))
    assert_allclose(test_output['x'].value, known_displacements['x'], rtol=5e-2, atol=0)
    assert_allclose(test_output['y'].value, known_displacements['y'], rtol=5e-2, atol=0 )
    assert(isinstance(test_output['x'], u.Quantity) and test_output['x'].unit == 'arcsec')
    assert(isinstance(test_output['y'], u.Quantity) and test_output['y'].unit == 'arcsec')

    # Test returning using the "with_displacements" option
    test_output = mapcube_solar_derotate(mc, with_displacements=True)
    # Assert
    assert(isinstance(test_output[0], map.MapCube))
    assert(isinstance(test_output[1], dict))
    assert_allclose(test_output[1]['x'], known_displacements['x'], rtol=5e-2, atol=0)
    assert_allclose(test_output[1]['y'], known_displacements['y'], rtol=5e-2, atol=0 )
    assert(isinstance(test_output[1]['x'], u.Quantity) and test_output[1]['x'].unit == 'arcsec')
    assert(isinstance(test_output[1]['y'], u.Quantity) and test_output[1]['y'].unit == 'arcsec')
