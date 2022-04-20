import numpy as np
import pytest

import astropy.units as u

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import asdf_entry_points
from .helpers import roundtrip_object


def assert_roundtrip_map(old):
    new = roundtrip_object(old)
    np.testing.assert_allclose(old.data, new.data)

    # Test the meta by force!
    for ok, ov in old.meta.items():
        assert ok in new.meta
        assert new.meta[ok] == ov

    assert u.allclose(old.shifted_value, new.shifted_value)
    if old.mask is not None and new.mask is not None:
        np.testing.assert_allclose(old.mask, new.mask)
    assert old.unit == new.unit


@pytest.fixture
def aia171_test_map():
    aia_path = get_test_filepath("aia_171_level1.fits")
    return sunpy.map.Map(aia_path)


@asdf_entry_points
def test_genericmap_basic(aia171_test_map, tmpdir):
    assert_roundtrip_map(aia171_test_map)


@asdf_entry_points
def test_genericmap_mask(aia171_test_map, tmpdir):

    mask = np.zeros_like(aia171_test_map.data)
    mask[10, 10] = 1

    aia171_test_map.mask = mask
    aia171_test_map._unit = u.m

    assert_roundtrip_map(aia171_test_map)
