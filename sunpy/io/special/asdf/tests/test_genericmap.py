import numpy as np
import pytest
from packaging.version import Version

import asdf
import astropy.units as u
from asdf.testing.helpers import roundtrip_object

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import asdf_entry_points


def assert_roundtrip_map(old):
    new = roundtrip_object(old)
    np.testing.assert_allclose(old.data, new.data)

    # Test the meta by force!
    for ok, ov in old.meta.items():
        assert ok in new.meta
        assert new.meta[ok] == ov

    if old.mask is not None and new.mask is not None:
        np.testing.assert_allclose(old.mask, new.mask)

    assert old.unit == new.unit


def asdf_open_memory_mapping_kwarg(memmap: bool) -> dict:
    if Version(asdf.__version__) >= Version("3.1.0"):
        return {"memmap": memmap}
    else :
        return {"copy_arrays": not memmap}


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


@asdf_entry_points
def test_load_100_file_with_shift():
    fname = get_test_filepath("aiamap_shift_genericmap_1.0.0.asdf")
    with asdf.open(fname, **asdf_open_memory_mapping_kwarg(memmap=False)) as af:
        aiamap = af['object']
        assert isinstance(aiamap, sunpy.map.sources.AIAMap)
        assert "crval1" in aiamap.meta.modified_items
        crval1 = aiamap.meta.modified_items["crval1"]
        assert crval1.current - crval1.original == 10


@asdf_entry_points
def test_load_100_file_with_no_shift():
    fname = get_test_filepath("aiamap_genericmap_1.0.0.asdf")
    with asdf.open(fname, **asdf_open_memory_mapping_kwarg(memmap=False)) as af:
        aiamap = af['object']
        assert isinstance(aiamap, sunpy.map.sources.AIAMap)
        assert "crval1" not in aiamap.meta.modified_items
        assert "crval2" not in aiamap.meta.modified_items
