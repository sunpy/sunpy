import os

import numpy as np
import pytest

import asdf
import astropy.units as u
from asdf.testing.helpers import roundtrip_object
from astropy.coordinates import CartesianRepresentation

import sunpy.coordinates.frames as frames
from sunpy.tests.helpers import asdf_entry_points

sunpy_frames = list(map(lambda name: getattr(frames, name), frames.__all__))
# Don't test the two base frames
sunpy_frames = [frame for frame in sunpy_frames if 'base' not in frame.name]


def assert_round_trip_frame(old):
    new = roundtrip_object(old)
    assert isinstance(new, type(old))
    if new.has_data:
        assert new.data.components == old.data.components
        for comp in new.data.components:
            assert u.allclose(getattr(new.data, comp), getattr(old.data, comp))


@pytest.fixture(params=sunpy_frames)
def coordframe_scalar(request):
    frame = request.param

    if frame._default_representation is CartesianRepresentation:
        data = np.random.random(3) * u.km
    else:
        data = np.random.random(2) * u.arcsec

    return frame(*data, obstime='2018-01-01T00:00:00')


@pytest.fixture(params=sunpy_frames)
def coordframe_array(request):
    frame = request.param

    if frame._default_representation is CartesianRepresentation:
        data = np.random.random((3, 10)) * u.km
    else:
        data = np.random.random((2, 10)) * u.arcsec

    return frame(*data, obstime='2018-01-01T00:00:00')


# Ignore warnings thrown when trying to load the ASDF in a different
# astropy version to that with which it was created.
@pytest.mark.filterwarnings('ignore:.*was created with extension.*')
@asdf_entry_points
def test_hgc_100():
    # Test that HeliographicCarrington is populated with Earth as the
    # observer when loading a older schema (1.0.0)
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
    assert_round_trip_frame(coordframe_scalar)


@asdf_entry_points
def test_saveframe_arr(coordframe_array, tmpdir):
    assert_round_trip_frame(coordframe_array)


@asdf_entry_points
def test_hpc_observer_version(tmpdir):
    """
    This test verifies that the helioprojective frame has an upto date list of
    the HGS schema versions by ensuring that a HPC frame with a instantiated
    HGS frame as the observer round trips correctly.
    """
    time = "2021-10-13T11:08"
    obs = frames.HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU, obstime=time)
    coord = frames.Helioprojective(10*u.arcsec, 10*u.arcsec, obstime=time, observer=obs)
    assert_round_trip_frame(coord)


@asdf_entry_points
def test_hcc_observer_version(tmpdir):
    """
    This test verifies that the heliocentric frame has an upto date list of
    the HGS schema versions by ensuring that a HPC frame with a instantiated
    HGS frame as the observer round trips correctly.
    """
    time = "2021-10-13T11:08"
    obs = frames.HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU, obstime=time)
    coord = frames.Heliocentric(1*u.Mm, 1*u.Mm, 1*u.Mm, obstime=time, observer=obs)
    assert_round_trip_frame(coord)
