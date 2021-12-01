"""
isort:skip_file.
"""
import numpy as np
import pytest
import sys

import astropy.units as u
from asdf.tests.helpers import assert_roundtrip_tree

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import asdf_entry_points
from sunpy.io.special.asdf.extension import SunpyExtension

# TODO: Delete after a major pytest release
if sys.version_info > (3, 9):
    pytest.skip("pytest + asdf do not play well", allow_module_level=True)


@pytest.fixture
def aia171_test_map():
    aia_path = get_test_filepath("aia_171_level1.fits")
    return sunpy.map.Map(aia_path)


@asdf_entry_points
def test_genericmap_basic(aia171_test_map, tmpdir):

    tree = {'smap': aia171_test_map}

    assert_roundtrip_tree(tree, tmpdir, extensions=SunpyExtension())


@asdf_entry_points
def test_genericmap_mask(aia171_test_map, tmpdir):

    mask = np.zeros_like(aia171_test_map.data)
    mask[10, 10] = 1

    aia171_test_map.mask = mask
    aia171_test_map._unit = u.m

    tree = {'smap': aia171_test_map}

    assert_roundtrip_tree(tree, tmpdir, extensions=SunpyExtension())
