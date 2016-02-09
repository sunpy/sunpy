# -*- coding: utf-8 -*-
"""
Test mapcube functionality
"""
from __future__ import absolute_import

import numpy as np
import astropy.units as u
import sunpy
import sunpy.map
from sunpy.map.header import MapMeta
import pytest
import os
import sunpy.data.test


@pytest.fixture
def aia_map():
    """
    Load SunPy's test AIA image.
    """
    testpath = sunpy.data.test.rootdir
    aia_file = os.path.join(testpath, "aia_171_level1.fits")
    return sunpy.map.Map(aia_file)


@pytest.fixture
def masked_aia_map(aia_map):
    """
    Put a simple mask in the test AIA image.  A rectangular (not square) block
    of True values are included to test that operations on the mask respect how
    the mask is stored.
    """
    aia_map_data = aia_map.data
    aia_map_mask = np.zeros_like(aia_map_data)
    aia_map_mask[0:2, 0:3] = True
    return sunpy.map.Map(np.ma.masked_array(aia_map_data, mask=aia_map_mask),
                         aia_map.meta)


@pytest.fixture
def mapcube_all_the_same(aia_map):
    """ Simple `sunpy.map.mapcube` for testing."""
    return sunpy.map.Map([aia_map, aia_map], cube=True)


@pytest.fixture
def mapcube_all_the_same_all_have_masks(masked_aia_map):
    """ Simple `sunpy.map.mapcube` for testing, in which all the maps have
    masks."""
    return sunpy.map.Map([masked_aia_map, masked_aia_map], cube=True)


@pytest.fixture
def mapcube_all_the_same_some_have_masks(aia_map, masked_aia_map):
    """ Simple `sunpy.map.mapcube` for testing, in which at least some of the
    maps have masks."""
    return sunpy.map.Map([masked_aia_map, masked_aia_map, aia_map], cube=True)


@pytest.fixture()
def mapcube_different(aia_map):
    """ Mapcube allows that the size of the image data in each map be
    different.  This mapcube contains such maps."""
    return sunpy.map.Map([aia_map, aia_map.superpixel((4, 4)*u.pix)], cube=True)


def test_all_maps_same_shape(mapcube_all_the_same, mapcube_different):
    """Make sure that Mapcube knows if all the maps have the same shape"""
    assert mapcube_all_the_same.all_maps_same_shape()
    assert not mapcube_different.all_maps_same_shape()


def test_at_least_one_map_has_mask(mapcube_all_the_same,
                                   mapcube_all_the_same_all_have_masks,
                                   mapcube_all_the_same_some_have_masks
                                   ):
    """ Test that we can detect the presence of at least one masked map."""
    assert not mapcube_all_the_same.at_least_one_map_has_mask()
    assert mapcube_all_the_same_all_have_masks.at_least_one_map_has_mask()
    assert mapcube_all_the_same_some_have_masks.at_least_one_map_has_mask()


def test_as_array(mapcube_all_the_same,
                  mapcube_different,
                  mapcube_all_the_same_all_have_masks,
                  mapcube_all_the_same_some_have_masks):
    """Make sure the data in the mapcube returns correctly, when all the
    maps have the same shape.  When they don't have the same shape, make
    sure an error is raised."""
    # Should raise a ValueError if the mapcube has differently shaped maps in
    # it.
    with pytest.raises(ValueError):
        mapcube_different.as_array()

    # Test the case when none of the maps have a mask
    returned_array = mapcube_all_the_same.as_array()
    assert isinstance(returned_array, np.ndarray)
    assert returned_array.ndim == 3
    assert len(returned_array.shape) == 3
    assert returned_array.shape[0] == 128
    assert returned_array.shape[1] == 128
    assert returned_array.shape[2] == 2
    assert np.ma.getmask(returned_array) is np.ma.nomask

    # Test the case when all the maps have masks
    returned_array = mapcube_all_the_same_all_have_masks.as_array()
    assert isinstance(returned_array, np.ma.masked_array)
    data = np.ma.getdata(returned_array)
    assert data.ndim == 3
    assert len(data.shape) == 3
    assert data.shape[0] == 128
    assert data.shape[1] == 128
    assert data.shape[2] == 2
    mask = np.ma.getmask(returned_array)
    assert mask.ndim == 3
    assert len(mask.shape) == 3
    assert mask.shape[0] == 128
    assert mask.shape[1] == 128
    assert mask.shape[2] == 2
    assert mask.dtype == bool

    # Test the case when some of the maps have masks
    returned_array = mapcube_all_the_same_some_have_masks.as_array()
    assert isinstance(returned_array, np.ma.masked_array)
    data = np.ma.getdata(returned_array)
    assert data.ndim == 3
    assert len(data.shape) == 3
    assert data.shape[0] == 128
    assert data.shape[1] == 128
    assert data.shape[2] == 3
    mask = np.ma.getmask(mapcube_all_the_same_some_have_masks.as_array())
    assert mask.ndim == 3
    assert len(mask.shape) == 3
    assert mask.shape[0] == 128
    assert mask.shape[1] == 128
    assert mask.shape[2] == 3
    assert np.all(mask[0:2, 0:3, 0])
    assert np.all(mask[0:2, 0:3, 1])
    assert np.all(np.logical_not(mask[0:2, 0:3, 2]))


def test_all_meta(mapcube_all_the_same):
    """Tests that the correct number of map meta objects are returned, and
    that they are all map meta objects."""
    meta = mapcube_all_the_same.all_meta()
    assert len(meta) == 2
    assert np.all(np.asarray([isinstance(h, MapMeta) for h in meta]))
    assert np.all(np.asarray([meta[i] == mapcube_all_the_same[i].meta for i in range(0, len(meta))]))
