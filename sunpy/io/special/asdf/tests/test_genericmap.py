import numpy as np
import pytest

import astropy.units as u
from asdf.testing.helpers import roundtrip_object

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import asdf_entry_points


@pytest.fixture
def aia171_test_map():
    aia_path = get_test_filepath("aia_171_level1.fits")
    return sunpy.map.Map(aia_path)


@asdf_entry_points
def test_genericmap_basic(aia171_test_map, tmpdir):
    new_map = roundtrip_object(aia171_test_map)


@asdf_entry_points
def test_genericmap_mask(aia171_test_map, tmpdir):

    mask = np.zeros_like(aia171_test_map.data)
    mask[10, 10] = 1

    aia171_test_map.mask = mask
    aia171_test_map._unit = u.m

    new_map = roundtrip_object(aia171_test_map)
