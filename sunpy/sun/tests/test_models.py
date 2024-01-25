import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import Longitude, SkyCoord
from astropy.tests.helper import assert_quantity_allclose

import sunpy.map
from sunpy.sun.models import diff_rot


@pytest.fixture
def all_off_disk_map(aia171_test_map):
    return aia171_test_map.submap((1, 1)*u.pix, top_right=(11, 12)*u.pix)


@pytest.fixture
def all_on_disk_map(aia171_test_map):
    return aia171_test_map.submap((30, 60)*u.pix, top_right=(50, 85)*u.pix)


@pytest.fixture
def straddles_limb_map(aia171_test_map):
    return aia171_test_map.submap((64, 80)*u.pix, top_right=(120, 127)*u.pix)


@pytest.fixture
def aia171_test_map_with_mask(aia171_test_map):
    shape = aia171_test_map.data.shape
    mask = np.zeros_like(aia171_test_map.data, dtype=bool)
    mask[0:shape[0]//2, 0:shape[1]//2] = True
    return sunpy.map.Map(np.ma.array(aia171_test_map.data, mask=mask), aia171_test_map.meta)


@pytest.fixture
def aia171_test_submap(aia171_test_map):
    bl = SkyCoord(-512 * u.arcsec, 100 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    ur = SkyCoord(-100 * u.arcsec, 400 * u.arcsec, frame=aia171_test_map.coordinate_frame)
    return aia171_test_map.submap(bl, top_right=ur)


@pytest.fixture
def seconds_per_day():
    return 24 * 60 * 60.0 * u.s




def test_single(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg)
    assert_quantity_allclose(rot, 136.8216 * u.deg, rtol=1e-3)


def test_array(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, np.linspace(-70, 70, 2) * u.deg)
    assert_quantity_allclose(rot, Longitude(np.array([110.2725, 110.2725]) * u.deg), rtol=1e-3)


def test_synodic(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='howard', frame_time='synodic')
    assert_quantity_allclose(rot, 126.9656 * u.deg, rtol=1e-3)


def test_sidereal(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='howard', frame_time='sidereal')
    assert_quantity_allclose(rot, 136.8216 * u.deg, rtol=1e-3)


def test_howard(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='howard')
    assert_quantity_allclose(rot, 136.8216 * u.deg, rtol=1e-3)


def test_allen(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='allen')
    assert_quantity_allclose(rot, 136.9 * u.deg, rtol=1e-3)


def test_snodgrass(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='snodgrass')
    assert_quantity_allclose(rot, 135.4232 * u.deg, rtol=1e-3)


def test_rigid(seconds_per_day):
    rot = diff_rot(10 * seconds_per_day, [0, 30, 60] * u.deg, rot_type='rigid')
    assert_quantity_allclose(rot, [141.844 * u.deg] * 3, rtol=1e-3)


def test_fail(seconds_per_day):
    with pytest.raises(ValueError):
        diff_rot(10 * seconds_per_day, 30 * u.deg, rot_type='garbage')
