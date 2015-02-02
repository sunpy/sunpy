from __future__ import absolute_import

import pytest

import numpy as np
from astropy import units as u
from astropy.coordinates import Longitude

from sunpy.physics.transforms.differential_rotation import diff_rot
from sunpy.tests.helpers import assert_quantity_allclose

# Please note the numbers in these tests are not checked for physical
# accuracy, only that they are the values the function was outputting upon
# implementation.

def test_single():
    rot = diff_rot(10, 30 * u.deg)
    assert rot == 136.8216 * u.deg

def test_array():
    rot = diff_rot(10, np.linspace(-70, 70, 2) * u.deg)
    assert_quantity_allclose(rot, Longitude(np.array([110.2725,  110.2725]) * u.deg))

def test_synodic():
    rot = diff_rot(10, 30 * u.deg, rot_type='howard', frame_time='synodic')
    assert rot == 126.9656 * u.deg

def test_sidereal():
    rot = diff_rot(10, 30 * u.deg, rot_type='howard', frame_time='sidereal')
    assert rot == 136.8216 * u.deg

def test_howard():
    rot = diff_rot(10, 30 * u.deg, rot_type='howard')
    assert rot == 136.8216 * u.deg

def test_allen():
    rot = diff_rot(10, 30 * u.deg, rot_type='allen')
    assert rot == 136.9 * u.deg

def test_snodgrass():
    rot = diff_rot(10, 30 * u.deg, rot_type='snodgrass')
    assert rot == 135.4232 * u.deg

def test_fail():
    with pytest.raises(ValueError):
        rot = diff_rot(10, 30 * u.deg, rot_type='garbage')
