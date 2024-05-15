"""
Test mapsequence functionality
"""

import numpy as np
import pytest

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose
from astropy.visualization import ImageNormalize

import sunpy
import sunpy.data.test
import sunpy.map
from sunpy.tests.helpers import figure_test
from sunpy.util.metadata import MetaDict


@pytest.fixture
def mapsequence_all_the_same(aia171_test_map):
    """ Simple `sunpy.map.mapsequence` for testing."""
    return sunpy.map.Map([aia171_test_map, aia171_test_map], sequence=True)


@pytest.fixture
def mapsequence_different_maps(aia171_test_map, eit_test_map):
    """
    Simple `sunpy.map.mapsequence` for testing, in which there are different maps
    """
    return sunpy.map.Map([aia171_test_map, eit_test_map], sequence=True)


@pytest.fixture
def mapsequence_all_the_same_all_have_masks(aia171_test_map_with_mask):
    """ Simple `sunpy.map.mapsequence` for testing, in which all the maps have
    masks."""
    return sunpy.map.Map([aia171_test_map_with_mask, aia171_test_map_with_mask], sequence=True)


@pytest.fixture
def mapsequence_all_the_same_some_have_masks(aia171_test_map, aia171_test_map_with_mask):
    """ Simple `sunpy.map.mapsequence` for testing, in which at least some of the
    maps have masks."""
    return sunpy.map.Map([aia171_test_map_with_mask, aia171_test_map_with_mask, aia171_test_map], sequence=True)


@pytest.fixture()
def mapsequence_different(aia171_test_map):
    """ Mapsequence allows that the size of the image data in each map be
    different.  This mapsequence contains such maps."""
    return sunpy.map.Map([aia171_test_map, aia171_test_map.superpixel((4, 4) * u.pix)], sequence=True)


def test_all_maps_same_shape(mapsequence_all_the_same, mapsequence_different):
    """Make sure that Mapsequence knows if all the maps have the same shape"""
    assert mapsequence_all_the_same.all_maps_same_shape()
    assert not mapsequence_different.all_maps_same_shape()


def test_at_least_one_map_has_mask(mapsequence_all_the_same,
                                   mapsequence_all_the_same_all_have_masks,
                                   mapsequence_all_the_same_some_have_masks
                                   ):
    """ Test that we can detect the presence of at least one masked map."""
    assert not mapsequence_all_the_same.at_least_one_map_has_mask()
    assert mapsequence_all_the_same_all_have_masks.at_least_one_map_has_mask()
    assert mapsequence_all_the_same_some_have_masks.at_least_one_map_has_mask()


def test_as_array(mapsequence_all_the_same,
                  mapsequence_different,
                  mapsequence_all_the_same_all_have_masks,
                  mapsequence_all_the_same_some_have_masks):
    """Make sure the data in the mapsequence returns correctly, when all the
    maps have the same shape.  When they don't have the same shape, make
    sure an error is raised."""
    # Should raise a ValueError if the mapsequence has differently shaped maps in
    # it.
    with pytest.raises(ValueError):
        mapsequence_different.as_array()

    # Test the case when none of the maps have a mask
    returned_array = mapsequence_all_the_same.as_array()
    assert isinstance(returned_array, np.ndarray)
    assert returned_array.ndim == 3
    assert len(returned_array.shape) == 3
    assert returned_array.shape[0] == 128
    assert returned_array.shape[1] == 128
    assert returned_array.shape[2] == 2
    assert np.ma.getmask(returned_array) is np.ma.nomask

    # Test the case when all the maps have masks
    returned_array = mapsequence_all_the_same_all_have_masks.as_array()
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
    returned_array = mapsequence_all_the_same_some_have_masks.as_array()
    assert isinstance(returned_array, np.ma.masked_array)
    data = np.ma.getdata(returned_array)
    assert data.ndim == 3
    assert len(data.shape) == 3
    assert data.shape[0] == 128
    assert data.shape[1] == 128
    assert data.shape[2] == 3
    mask = np.ma.getmask(mapsequence_all_the_same_some_have_masks.as_array())
    assert mask.ndim == 3
    assert len(mask.shape) == 3
    assert mask.shape[0] == 128
    assert mask.shape[1] == 128
    assert mask.shape[2] == 3
    assert np.all(mask[0:2, 0:3, 0])
    assert np.all(mask[0:2, 0:3, 1])
    assert np.all(np.logical_not(mask[0:2, 0:3, 2]))


def test_all_meta(mapsequence_all_the_same):
    """Tests that the correct number of map meta objects are returned, and
    that they are all map meta objects."""
    meta = mapsequence_all_the_same.all_meta()
    assert len(meta) == 2
    assert np.all(np.asarray([isinstance(h, MetaDict) for h in meta]))
    assert np.all(np.asarray(
        [meta[i] == mapsequence_all_the_same[i].meta for i in range(0, len(meta))]))


def test_repr(mapsequence_all_the_same, mapsequence_different_maps):
    """
    Tests that overridden __repr__ functionality works as expected. Test
    for mapsequence of same maps as well that of different maps.
    """
    # Test the case of MapSequence having same maps
    expected_out = 'MapSequence of 2 elements, with maps from AIAMap'
    obtained_out = repr(mapsequence_all_the_same)
    assert obtained_out.startswith(object.__repr__(mapsequence_all_the_same))
    assert len(mapsequence_all_the_same) == 2
    assert expected_out in obtained_out

    # Test the case of MapSequence having different maps
    expected_out1 = 'MapSequence of 2 elements, with maps from AIAMap, EITMap'
    expected_out2 = 'MapSequence of 2 elements, with maps from EITMap, AIAMap'
    obtained_out = repr(mapsequence_different_maps)
    assert obtained_out.startswith(object.__repr__(mapsequence_different_maps))
    assert len(mapsequence_different_maps) == 2
    assert expected_out1 in obtained_out or expected_out2 in obtained_out


def test_derotate():
    with pytest.raises(NotImplementedError):
        sunpy.map.MapSequence(derotate=True)


def test_repr_html(mapsequence_all_the_same):
    html_string = mapsequence_all_the_same._repr_html_()
    for m in mapsequence_all_the_same.maps:
        assert m._repr_html_() in html_string


def test_quicklook(mocker, mapsequence_all_the_same):
    mockwbopen = mocker.patch('webbrowser.open_new_tab')
    mapsequence_all_the_same.quicklook()
    # Check that the mock web browser was opened with a file URL
    mockwbopen.assert_called_once()
    file_url = mockwbopen.call_args[0][0]
    assert file_url.startswith('file://')
    # Open the file specified in the URL and confirm that it contains the HTML
    with open(file_url[7:]) as f:
        html_string = f.read()
    for m in mapsequence_all_the_same.maps:
        assert m._repr_html_() in html_string


@figure_test
def test_norm_animator(hmi_test_map):
    seq = sunpy.map.Map([hmi_test_map, hmi_test_map], sequence=True)
    ani = seq.peek()
    ani._slider_changed(1, ani.sliders[ani.active_slider]._slider)
    return ani.fig


@figure_test
def test_map_sequence_plot(aia171_test_map, hmi_test_map):
    seq = sunpy.map.Map([aia171_test_map, hmi_test_map], sequence=True)
    seq.plot()


@figure_test
def test_map_sequence_plot_custom_cmap_norm(aia171_test_map, hmi_test_map):
    seq = sunpy.map.Map([aia171_test_map, hmi_test_map], sequence=True)
    animation = seq.plot(cmap='Greys', norm=ImageNormalize(vmin=0, vmax=100))
    animation._step()


def test_save(aia171_test_map, hmi_test_map, tmp_path):
    """
    Tests the MapSequence save function
    """
    seq = sunpy.map.Map([aia171_test_map, hmi_test_map], sequence=True)

    with pytest.raises(ValueError, match="'{index}' must be appear in the string"):
        seq.save("index", filetype='fits', overwrite=True)

    base_str = (tmp_path / "map").as_posix()
    seq.save(f"{base_str}_{{index:03}}.fits", filetype='auto', overwrite=True)

    test_seq = sunpy.map.Map(f"{base_str}_000.fits",
                             f"{base_str}_001.fits", sequence=True)

    assert isinstance(test_seq.maps[0], sunpy.map.sources.sdo.AIAMap)
    assert isinstance(test_seq.maps[1], sunpy.map.sources.sdo.HMIMap)

    assert test_seq.maps[0].meta.keys() == seq.maps[0].meta.keys()
    for k in seq.maps[0].meta:
        assert test_seq.maps[0].meta[k] == seq.maps[0].meta[k]
    assert_quantity_allclose(test_seq.maps[0].data, seq.maps[0].data)

    assert test_seq.maps[1].meta.keys() == seq.maps[1].meta.keys()
    for k in seq.maps[1].meta:
        assert test_seq.maps[1].meta[k] == seq.maps[1].meta[k]
    assert_quantity_allclose(test_seq.maps[1].data, seq.maps[1].data)


@figure_test
def test_map_sequence_plot_clip_interval(aia171_test_map):
    seq = sunpy.map.Map([aia171_test_map, aia171_test_map, aia171_test_map], sequence=True)
    # We want to blow the image out to make sure it is clipped on the test data
    animation = seq.plot(clip_interval=(5,75)*u.percent)
    animation._step()
