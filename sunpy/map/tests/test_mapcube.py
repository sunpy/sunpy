# -*- coding: utf-8 -*-
"""
Test mapcube functionality
"""
from __future__ import absolute_import

import numpy as np
import sunpy
import sunpy.map
from sunpy.map.header import MapMeta
import pytest

# Test maps and mapcubes
aia_map = sunpy.map.Map(sunpy.AIA_171_IMAGE)
mapcube_all_the_same = sunpy.map.Map([aia_map.superpixel((2, 1)),
                                      aia_map.superpixel((2, 1))], cube=True)
mapcube_different = sunpy.map.Map([aia_map,
                                   aia_map.superpixel((4, 4))], cube=True)


def test_all_maps_same_shape():
    """Make sure that Mapcube knows if all the maps have the same shape"""
    assert mapcube_all_the_same.all_maps_same_shape()
    assert not mapcube_different.all_maps_same_shape()


def test_as_array():
    """Make sure the data in the mapcube returns correctly, when all the
    maps have the same shape.  When they don't have the same shape, make
    sure an error is raised."""
    returned_array = mapcube_all_the_same.as_array()
    assert isinstance(returned_array, np.ndarray)
    assert returned_array.ndim == 3
    assert len(returned_array.shape) == 3
    assert returned_array.shape[0] == 1024
    assert returned_array.shape[1] == 512
    assert returned_array.shape[2] == 2
    # Should raise a ValueError if the mapcube has differently shaped maps in
    # it.
    with pytest.raises(ValueError):
        mapcube_different.as_array()


def test_all_meta():
    """Tests that the correct number of map meta objects are returned, and
    that they are all map meta objects."""
    meta = mapcube_all_the_same.all_meta()
    assert len(meta) == 2
    assert np.all(np.asarray([isinstance(h, MapMeta) for h in meta]))
    assert np.all(np.asarray([meta[i] == mapcube_all_the_same[i].meta for i in range(0, len(meta))]))
