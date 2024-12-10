import pickle

import pytest
from hypothesis import given, settings

import astropy.units as u
from astropy.coordinates import HeliocentricMeanEcliptic, SkyCoord, frame_transform_graph
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time, TimeDelta

import sunpy.coordinates.frames as f
from sunpy.coordinates.metaframes import RotatedSunFrame, _rotatedsun_cache
from sunpy.coordinates.tests.helpers import assert_longitude_allclose
from sunpy.coordinates.tests.strategies import latitudes, longitudes, times
from sunpy.sun import constants
from sunpy.sun.models import differential_rotation

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
    assert type(rot_frame.base) == base_class  # NOQA: E721

    # Check that the new class does *not* have the `obstime` frame attribute
    assert 'obstime' not in rot_frame.frame_attributes

    # Check that one-leg transformations have been created
    assert len(frame_transform_graph.get_transform(rot_class, rot_class).transforms) == 1
    assert len(frame_transform_graph.get_transform(f.HeliographicStonyhurst, rot_class).transforms) == 1
    assert len(frame_transform_graph.get_transform(rot_class, f.HeliographicStonyhurst).transforms) == 1

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

    assert type(a) == type(rot_hgs[1].base)  # NOQA: E721

    assert_longitude_allclose(a.lon, rot_hgs[1].lon)
    assert_quantity_allclose(a.lat, rot_hgs[1].lat)
    assert_quantity_allclose(a.radius, rot_hgs[1].radius)


def test_no_base():
    with pytest.raises(TypeError):
        RotatedSunFrame()


def test_no_obstime():
    with pytest.raises(ValueError, match="The base coordinate frame must have a defined `obstime`"):
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


def test_duration_from_timedelta():
    base_frame = f.HeliographicStonyhurst(obstime='2001-01-01')

    duration_timedelta = TimeDelta(4 * u.day)
    r = RotatedSunFrame(base=base_frame, duration=duration_timedelta)

    # Verify that the duration is correctly converted to a quantity in days
    assert_quantity_allclose(r.duration, 4 * u.day)


def test_duration_with_quantity_hours():
    base_frame = f.HeliographicStonyhurst(obstime='2001-01-01')

    # Testing with Quantity in hours (conversion needed)
    duration_quantity = 96 * u.hour  # 4 days in hours
    r = RotatedSunFrame(base=base_frame, duration=duration_quantity)
    assert_quantity_allclose(r.duration, 4 * u.day)


def test_both_duration_and_rotated_time_provided():
    base_frame = f.HeliographicStonyhurst(obstime='2001-01-01')

    with pytest.raises(ValueError, match="Specify either `duration` or `rotated_time`, not both."):
        RotatedSunFrame(base=base_frame, duration=TimeDelta(1*u.day), rotated_time=Time('2001-01-02'))


def test_duration_calculation():

    base_time = Time('2001-01-01')
    base_frame = f.HeliographicStonyhurst(obstime=base_time)
    rotated_time = Time('2001-01-02')
    r = RotatedSunFrame(base=base_frame, rotated_time=rotated_time)
    expected_duration = (rotated_time.utc - base_time).to('day')
    assert r.duration == expected_duration


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

    assert type(r) == type(rot_hgs[1])  # NOQA: E721
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
@settings(deadline=None, max_examples=10)
def test_rotatedsun_transforms(frame, lon, lat, obstime, rotated_time1, rotated_time2):
    # Tests the transformations (to, from, and loopback) for consistency with `differential_rotation` output

    if hasattr(frame, 'observer'):
        base = frame(lon=lon, lat=lat, observer='earth', obstime=obstime)
    else:
        base = frame(lon=lon, lat=lat, obstime=obstime)

    # Map any 2D coordinate to the surface of the Sun
    if base.spherical.distance.unit is u.one and u.allclose(base.spherical.distance, 1*u.one):
        newrepr = base.spherical * constants.radius
        base = base.realize_frame(newrepr)

    # Test the RotatedSunFrame->base transformation
    rsf1 = RotatedSunFrame(base=base, rotated_time=rotated_time1)
    result1 = rsf1.transform_to(base)

    desired_delta_lon1 = differential_rotation((rotated_time1 - obstime).to(u.day), lat)

    assert_longitude_allclose(result1.lon, rsf1.lon + desired_delta_lon1, atol=1e-5*u.deg)
    assert_quantity_allclose(base.lat, result1.lat, atol=1e-10*u.deg)
    # Use the `spherical` property since the name of the component varies with frame
    assert_quantity_allclose(base.spherical.distance, result1.spherical.distance)

    # Test the base->RotatedSunFrame transformation
    rsf2 = RotatedSunFrame(base=base, rotated_time=rotated_time2)
    result2 = base.transform_to(rsf2)

    desired_delta_lon2 = -differential_rotation((rotated_time2 - obstime).to(u.day), lat)

    assert_longitude_allclose(result2.lon, rsf2.lon + desired_delta_lon2, atol=1e-5*u.deg)
    assert_quantity_allclose(base.lat, result2.lat, atol=1e-10*u.deg)
    # Use the `spherical` property since the name of the component varies with frame
    assert_quantity_allclose(base.spherical.distance, result2.spherical.distance)

    # Test the RotatedSunFrame->RotatedSunFrame transformation
    result3 = rsf1.transform_to(rsf2)

    desired_delta_lon3 = desired_delta_lon1 + desired_delta_lon2

    assert_longitude_allclose(result3.lon, rsf1.lon + desired_delta_lon3, atol=1e-5*u.deg)
    assert_quantity_allclose(result3.lat, result1.lat, atol=1e-10*u.deg)
    # Use the `spherical` property since the name of the component varies with frame
    assert_quantity_allclose(result3.spherical.distance, result1.spherical.distance)


@pytest.mark.parametrize("indirect_fixture",
                         ["rot_hgs", "rot_hci"], indirect=True)
def test_obstime_change(indirect_fixture):
    base_class, rot_frame = indirect_fixture

    # Transforming a RotatedSun coordinate to a different obstime should give the same answer
    # as first transforming to its base and then transforming to the different obstime
    new_frame = base_class(obstime='2001-01-31')
    new_implicit = rot_frame.transform_to(new_frame)
    new_explicit = rot_frame.transform_to(rot_frame.base).transform_to(new_frame)

    assert_longitude_allclose(new_implicit.lon, new_explicit.lon, atol=1e-10*u.deg)


@pytest.mark.parametrize("indirect_fixture",
                         ["rot_hgs", "rot_hci"], indirect=True)
def test_obstime_change_loopback(indirect_fixture):
    base_class, rot_frame = indirect_fixture

    # Doing a explicit loopback transformation through an intermediate frame with a different
    # obstime should get back to the original coordinate
    int_frame = base_class(obstime='2001-01-31')
    loopback = rot_frame.transform_to(int_frame).transform_to(rot_frame)

    assert_longitude_allclose(loopback.lon, rot_frame.lon, atol=1e-10*u.deg)


@pytest.mark.parametrize("indirect_fixture",
                         ["rot_hgs", "rot_hgc", "rot_hci", "rot_hcc", "rot_hpc", "rot_hme"], indirect=True)
def test_transformation_to_nonobserver_frame(indirect_fixture):
    base_class, rot_frame = indirect_fixture

    hgs_frame = f.HeliographicStonyhurst(obstime='2020-01-01')
    hgs_coord = rot_frame.transform_to(hgs_frame)

    assert hgs_coord.obstime == hgs_frame.obstime


def test_pickle_rotatedsunframe():
    base_coord = SkyCoord(1*u.deg, 2*u.deg, obstime="2003-04-05", rsun=600*u.Mm,
                          frame='heliographic_stonyhurst')
    rotated_coord = RotatedSunFrame(base=base_coord, duration=7*u.day)
    pickled_coord = pickle.loads(pickle.dumps(rotated_coord))
    assert pickled_coord == rotated_coord
