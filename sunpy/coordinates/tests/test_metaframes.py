import pytest
from hypothesis import given, settings

import astropy.units as u
from astropy.coordinates import HeliocentricMeanEcliptic, SkyCoord, frame_transform_graph
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

import sunpy.coordinates.frames as f
from sunpy.coordinates.metaframes import RotatedSunFrame, _rotatedsun_cache
from sunpy.physics.differential_rotation import diff_rot
from .helpers import assert_longitude_allclose
from .strategies import latitudes, longitudes, times

# NorthOffsetFrame is tested in test_offset_frame.py


###########################
# Tests for RotatedSunFrame
###########################

@pytest.fixture
def indirect_fixture(request):
    return request.getfixturevalue(request.param)


@pytest.fixture
def rot_hgs():
    return (f.HeliographicStonyhurst,
            RotatedSunFrame(lon=1*u.deg, lat=2*u.deg, radius=3*u.AU,
                            base=f.HeliographicStonyhurst(obstime='2001-01-01'),
                            duration=4*u.day))


@pytest.fixture
def rot_hgc():
    return (f.HeliographicCarrington,
            RotatedSunFrame(lon=1*u.deg, lat=2*u.deg, radius=3*u.AU,
                            base=f.HeliographicCarrington(observer='earth', obstime='2001-01-01'),
                            duration=4*u.day))


@pytest.fixture
def rot_hci():
    return (f.HeliocentricInertial,
            RotatedSunFrame(lon=1*u.deg, lat=2*u.deg, distance=3*u.AU,
                            base=f.HeliocentricInertial(obstime='2001-01-01'),
                            duration=4*u.day))


@pytest.fixture
def rot_hcc():
    return (f.Heliocentric,
            RotatedSunFrame(x=1*u.AU, y=2*u.AU, z=3*u.AU,
                            base=f.Heliocentric(observer='earth', obstime='2001-01-01'),
                            duration=4*u.day))


@pytest.fixture
def rot_hpc():
    return (f.Helioprojective,
            RotatedSunFrame(Tx=1*u.deg, Ty=2*u.deg, distance=3*u.AU,
                            base=f.Helioprojective(observer='earth', obstime='2001-01-01'),
                            duration=4*u.day))


# Test a non-SunPy frame
@pytest.fixture
def rot_hme():
    return (HeliocentricMeanEcliptic,
            RotatedSunFrame(lon=1*u.deg, lat=2*u.deg, distance=3*u.AU,
                            base=HeliocentricMeanEcliptic(obstime='2001-01-01', equinox='2001-01-01'),
                            duration=4*u.day))


@pytest.mark.parametrize("indirect_fixture",
                         ["rot_hgs", "rot_hgc", "rot_hci", "rot_hcc", "rot_hpc", "rot_hme"], indirect=True)
def test_class_creation(indirect_fixture):
    base_class, rot_frame = indirect_fixture

    rot_class = type(rot_frame)

    # Check that the RotatedSunFrame class name has both 'RotatedSun' and the name of the base
    assert 'RotatedSun' in rot_class.__name__
    assert base_class.__name__ in rot_class.__name__

    # Check that the base class is in fact the specified class
    assert type(rot_frame.base) == base_class

    # Check that the new class does *not* have the `obstime` frame attribute
    assert 'obstime' not in rot_frame.frame_attributes

    # Check that one-leg transformations have been created
    assert len(frame_transform_graph.get_transform(rot_class, rot_class).transforms) == 1
    assert len(frame_transform_graph.get_transform(base_class, rot_class).transforms) == 1
    assert len(frame_transform_graph.get_transform(rot_class, base_class).transforms) == 1

    # Check that the base frame is in the cache
    assert base_class in _rotatedsun_cache

    # Check that the component data has been migrated
    assert rot_frame.has_data
    assert not rot_frame.base.has_data


def test_no_obstime_frame_attribute():
    assert 'obstime' not in RotatedSunFrame.frame_attributes


def test_as_base(rot_hgs):
    # Check the as_base() method
    a = rot_hgs[1].as_base()

    assert type(a) == type(rot_hgs[1].base)

    assert_longitude_allclose(a.lon, rot_hgs[1].lon)
    assert_quantity_allclose(a.lat, rot_hgs[1].lat)
    assert_quantity_allclose(a.radius, rot_hgs[1].radius)


def test_no_base():
    with pytest.raises(TypeError):
        RotatedSunFrame()


def test_no_obstime():
    with pytest.raises(ValueError):
        RotatedSunFrame(base=f.HeliographicStonyhurst(obstime=None))


def test_default_duration():
    r = RotatedSunFrame(base=f.HeliographicStonyhurst(obstime='2001-01-01'))
    assert_quantity_allclose(r.duration, 0*u.day)


def test_rotated_time_to_duration():
    r1 = RotatedSunFrame(base=f.HeliographicStonyhurst(obstime='2001-01-02'),
                         rotated_time='2001-01-03')
    assert_quantity_allclose(r1.duration, 1*u.day)

    r2 = RotatedSunFrame(base=f.HeliographicStonyhurst(obstime='2001-01-02'),
                         rotated_time='2001-01-01')
    assert_quantity_allclose(r2.duration, -1*u.day)


def test_rotated_time_property():
    r1 = RotatedSunFrame(base=f.HeliographicStonyhurst(obstime='2001-01-02'), duration=1*u.day)
    assert r1.rotated_time == Time('2001-01-03')

    r2 = RotatedSunFrame(base=f.HeliographicStonyhurst(obstime='2001-01-02'), duration=-1*u.day)
    assert r2.rotated_time == Time('2001-01-01')


def test_scalar_base_and_array_duration():
    scalar_base = f.HeliographicStonyhurst(1*u.deg, 2*u.deg, obstime='2001-01-02')
    array_duration = [1, -1]*u.day
    r = RotatedSunFrame(base=scalar_base, duration=array_duration)

    assert not r.data.isscalar
    assert r.data.shape == array_duration.shape
    assert_quantity_allclose(r.cartesian[0].xyz, scalar_base.cartesian.xyz)
    assert_quantity_allclose(r.cartesian[1].xyz, scalar_base.cartesian.xyz)


def test_base_skycoord(rot_hgs):
    # Check that RotatedSunFrame can be instantiated from a SkyCoord
    s = SkyCoord(1*u.deg, 2*u.deg, 3*u.AU, frame=f.HeliographicStonyhurst, obstime='2001-01-01')
    r = RotatedSunFrame(base=s)

    assert type(r) == type(rot_hgs[1])
    assert r.has_data
    assert not r.base.has_data

    assert_longitude_allclose(r.lon, s.lon)
    assert_quantity_allclose(r.lat, s.lat)
    assert_quantity_allclose(r.radius, s.radius)


def test_default_rotation_model():
    r = RotatedSunFrame(base=f.HeliographicStonyhurst(obstime='2001-01-01'))
    assert r.rotation_model == "howard"


def test_alternate_rotation_model():
    r = RotatedSunFrame(base=f.HeliographicStonyhurst(obstime='2001-01-01'),
                        rotation_model="allen")
    assert r.rotation_model == "allen"


@pytest.mark.parametrize("frame", [f.HeliographicStonyhurst,
                                   f.HeliographicCarrington,
                                   f.HeliocentricInertial])
@given(lon=longitudes(), lat=latitudes(),
       obstime=times(), rotated_time1=times(), rotated_time2=times())
@settings(deadline=None)
def test_rotatedsun_transforms(frame, lon, lat, obstime, rotated_time1, rotated_time2):
    # Tests the transformations (to, from, and loopback) for consistency with `diff_rot` output

    if hasattr(frame, 'observer'):
        base = frame(lon=lon, lat=lat, observer='earth', obstime=obstime)
    else:
        base = frame(lon=lon, lat=lat, obstime=obstime)

    # Test the RotatedSunFrame->base transformation
    rsf1 = RotatedSunFrame(base=base, rotated_time=rotated_time1)
    result1 = rsf1.transform_to(base)

    desired_delta_lon1 = diff_rot((rotated_time1 - obstime).to(u.day), lat)

    assert_longitude_allclose(result1.lon, rsf1.lon + desired_delta_lon1, atol=1e-5*u.deg)
    assert_quantity_allclose(base.lat, result1.lat)
    # Use the `spherical` property since the name of the component varies with frame
    assert_quantity_allclose(base.spherical.distance, result1.spherical.distance)

    # Test the base->RotatedSunFrame transformation
    rsf2 = RotatedSunFrame(base=base, rotated_time=rotated_time2)
    result2 = base.transform_to(rsf2)

    desired_delta_lon2 = -diff_rot((rotated_time2 - obstime).to(u.day), lat)

    assert_longitude_allclose(result2.lon, rsf2.lon + desired_delta_lon2, atol=1e-5*u.deg)
    assert_quantity_allclose(base.lat, result2.lat)
    # Use the `spherical` property since the name of the component varies with frame
    assert_quantity_allclose(base.spherical.distance, result2.spherical.distance)

    # Test the RotatedSunFrame->RotatedSunFrame transformation
    result3 = rsf1.transform_to(rsf2)

    desired_delta_lon3 = desired_delta_lon1 + desired_delta_lon2

    assert_longitude_allclose(result3.lon, rsf1.lon + desired_delta_lon3, atol=1e-5*u.deg)
    assert_quantity_allclose(result3.lat, result1.lat)
    # Use the `spherical` property since the name of the component varies with frame
    assert_quantity_allclose(result3.spherical.distance, result1.spherical.distance)
