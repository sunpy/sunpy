import os

import numpy as np
import pytest

import asdf
import astropy.units as u
from astropy.coordinates import CartesianRepresentation

import sunpy.coordinates.frames as frames
from sunpy.tests.helpers import asdf_entry_points

from asdf.tests.helpers import assert_roundtrip_tree  # NOQA isort:skip

sunpy_frames = list(map(lambda name: getattr(frames, name), frames.__all__))
# Don't test the two base frames
sunpy_frames = [frame for frame in sunpy_frames if 'base' not in frame.name]


@pytest.fixture(params=sunpy_frames)
@asdf_entry_points
def coordframe_scalar(request):
    frame = request.param

    if frame._default_representation is CartesianRepresentation:
        data = np.random.random(3) * u.km
    else:
        data = np.random.random(2) * u.arcsec

    return frame(*data, obstime='2018-01-01T00:00:00')


@pytest.fixture(params=sunpy_frames)
@asdf_entry_points
def coordframe_array(request):
    frame = request.param

    if frame._default_representation is CartesianRepresentation:
        data = np.random.random((3, 10)) * u.km
    else:
        data = np.random.random((2, 10)) * u.arcsec

    return frame(*data, obstime='2018-01-01T00:00:00')


# Ignore warnings thrown when trying to load the ASDF in a different astropy
# version to that with which it was created.
@pytest.mark.filterwarnings('ignore:.*was created with extension.*')
def test_hgc_100():
    # Test that HeliographicCarrington is populated with Earth as the observer when loading a
    # older schema (1.0.0)
    test_file = os.path.join(os.path.dirname(__file__), "hgc_100.asdf")
    with asdf.open(test_file) as input_asdf:
        hgc = input_asdf['hgc']
        assert isinstance(hgc, frames.HeliographicCarrington)
        if hgc.obstime is None:
            assert hgc.observer == 'earth'
        else:
            assert hgc.observer.object_name == 'earth'


@asdf_entry_points
def test_saveframe(coordframe_scalar, tmpdir):
    tree = {'frame': coordframe_scalar}
    assert_roundtrip_tree(tree, tmpdir)


@asdf_entry_points
def test_saveframe_arr(coordframe_array, tmpdir):
    tree = {'frame': coordframe_array}
    assert_roundtrip_tree(tree, tmpdir)
