import re

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import ConvertError, SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time, TimeDelta

from sunpy.coordinates import frames, get_earth, sun, transform_with_sun_center
from sunpy.coordinates.metaframes import RotatedSunFrame
from sunpy.coordinates.utils import (
    GreatArc,
    get_limb_coordinates,
    get_new_observer,
    get_rectangle_coordinates,
    solar_angle_equivalency,
    solar_coordinate_rotation,
)
from sunpy.sun import constants

# Test the great arc code against calculable quantities
# The inner angle is the same between each pair of coordinates. You can
# calculate these coordinates using the inner angle formulae as listed here:
# https://en.wikipedia.org/wiki/Great-circle_distance


@pytest.mark.parametrize(("start", "end"), [((0, 0), (0, 45)),
                                            ((0, 0), (45, 0)),
                                            ((0, 45), (0, 0)),
                                            ((45, 0), (0, 0)),
                                            ((12, 13), (12, 58)),
                                            ((-10, 6), (-10, 51)),
                                            ((-20, -50), (-20, -5)),
                                            ((10, -50), (87.53163324626676, -55))])
def test_great_arc_calculable(start, end):
    c = SkyCoord(start[0]*u.degree, start[1]*u.degree, frame=frames.HeliographicStonyhurst,
                 observer=frames.HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU))
    d = SkyCoord(end[0]*u.degree, end[1]*u.degree, frame=frames.HeliographicStonyhurst,
                 observer=frames.HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU))
    gc = GreatArc(c, d)

    c_trans = c.transform_to(frames.Heliocentric)
    assert gc.start.x == c_trans.x
    assert gc.start.y == c_trans.y
    assert gc.start.z == c_trans.z
    assert gc.start.observer.lat == 0*u.deg
    assert gc.start.observer.lon == 0*u.deg
    assert gc.start.observer.radius == 1 * u.AU

    d_trans = d.transform_to(frames.Heliocentric(observer=c.observer))
    assert gc.end.x == d_trans.x
    assert gc.end.y == d_trans.y
    assert gc.end.z == d_trans.z
    assert gc.end.observer.lat == 0*u.deg
    assert gc.end.observer.lon == 0*u.deg
    assert gc.end.observer.radius == 1 * u.AU

    np.testing.assert_almost_equal(gc.inner_angle.to('deg').value, 45.0)
    np.testing.assert_almost_equal(gc.radius.to('km').value, sun.constants.radius.to('km').value)
    np.testing.assert_almost_equal(gc.distance.to(
        'km').value, sun.constants.radius.to('km').value * 2 * np.pi/8, decimal=1)


# Test the calculation of coordinates using varying numbers of points on
# initialization of the GreatArc object.
@pytest.mark.parametrize(("points_requested", "points_expected", "first_point", "last_point", "last_inner_angle", "last_distance"),
                         # Test default
                         [(None, 100, (600, -600), (-100, 800), 1.8683580432741789, 1300377.1981299),
                          # Test int as an option
                          (3, 3, (600, -600), (-100, 800), 1.8683580432741789, 1300377.1981299),
                          # Test equally spaced monotonically increasing numpy
                          # array
                          (np.linspace(0, 1, 43), 43, (600, -600),
                           (-100, 800), 1.8683580432741789, 1300377.1981299),
                          # Test unequally spaced monotonically increasing numpy
                          # array
                          (np.asarray([0.1, 0.2, 0.6, 0.67, 0.99]), 5, (604.68091703, -468.64217597),
                           (-88.83212616, 792.76284375), 1.84967446, 1287373.4261486),
                          # Test unequally spaced monotonically decreasing numpy
                          # array
                          (np.asarray([0.93, 0.78, 0.3, 0.001]), 4, (-21.28208654, 743.58866798),
                           (600.1512768, -598.78376614), 0.00186836, 1300.37719813),
                          # Test numpy array that increases and decreases
                          (np.asarray([0.94, 0.73, 0.8, 0.21]), 4, (-32.5852606, 752.45507707),
                           (585.45829119, -305.26965043), 0.39235519, 273079.2116073)])
def test_great_arc_coordinates(points_requested, points_expected, first_point,
                               last_point, last_inner_angle, last_distance, aia171_test_map):
    coordinate_frame = aia171_test_map.coordinate_frame
    a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=coordinate_frame)
    b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=coordinate_frame)
    gc = GreatArc(a, b, points=points_requested)
    coordinates = gc.coordinates()
    inner_angles = gc.inner_angles()
    distances = gc.distances()

    # Ensure a GreatArc object is returned
    assert isinstance(gc, GreatArc)

    # Test the properties of the GreatArc object
    a_trans = a.transform_to(frames.Heliocentric)
    assert gc.start.x == a_trans.x
    assert gc.start.y == a_trans.y
    assert gc.start.z == a_trans.z
    b_trans = b.transform_to(frames.Heliocentric(observer=a.observer))
    assert gc.end.x == b_trans.x
    assert gc.end.y == b_trans.y
    assert gc.end.z == b_trans.z
    assert gc.distance_unit == u.m
    assert gc.observer == a.observer
    assert gc.center.x == 0 * u.m
    assert gc.center.y == 0 * u.m
    assert gc.center.z == 0 * u.m

    assert u.allclose(gc.start_cartesian * u.m, np.asarray(
        [428721.0913539, -428722.9051924, 341776.0910214]) * u.km)
    assert u.allclose(gc.end_cartesian * u.m, np.asarray(
        [-71429.5229381, 571439.071248, 390859.5797815]) * u.km)
    assert u.allclose(gc.center_cartesian * u.m, np.asarray([0, 0, 0]) * u.km)

    assert u.allclose(gc.v1 * u.m, np.asarray(
        [428721.0913539, -428722.9051924, 341776.0910214]) * u.km)
    assert u.allclose(gc._r, 696000000.0015007)
    assert u.allclose(gc.v2 * u.m, np.asarray(
        [-71429.5229381, 571439.071248, 390859.5797815]) * u.km)
    assert u.allclose(gc.v3 * u.m, np.asarray(
        [56761.6265851, 466230.7005856, 513637.0815867]) * u.km)

    # Inner angle
    assert gc.inner_angle.unit == u.rad
    np.testing.assert_almost_equal(gc.inner_angle.value, 1.8683580432741789)

    # Distance
    assert gc.distance.unit == u.m
    np.testing.assert_approx_equal(gc.distance.value, 1300377198.1299164)

    # Radius of the sphere
    assert gc.radius.unit == u.m
    assert u.isclose(gc.radius.value * u.m, 696000.000001501 * u.km)

    # Test the calculation of the SkyCoords
    # Coordinates method
    # Number of points
    assert len(coordinates) == points_expected

    # Start and end coordinates
    np.testing.assert_almost_equal(coordinates[0].Tx.value, first_point[0])
    np.testing.assert_almost_equal(coordinates[0].Ty.value, first_point[1])
    np.testing.assert_almost_equal(coordinates[-1].Tx.value, last_point[0])
    np.testing.assert_almost_equal(coordinates[-1].Ty.value, last_point[1])

    # Inner angles method
    # Inner angles
    assert len(inner_angles) == points_expected
    np.testing.assert_almost_equal(inner_angles[-1].value, last_inner_angle)

    # Distances method
    assert len(distances) == points_expected
    assert u.isclose(distances[-1].value * u.m, last_distance * u.km)

# Test that the great arc code rejects wrongly formatted points
@pytest.mark.parametrize(
    ("points", "expected_error"),
    [
        (np.asarray([[0, 0.1], [0.2, 0.3]]), "One dimensional numpy ndarrays only"),
        (np.asarray([0.1, 0.2, -0.1, 0.4]), "All value in points array must be strictly >=0 and <=1."),
        (np.asarray([0.3, 1.1, 0.6, 0.7]), "All value in points array must be strictly >=0 and <=1."),
        ('strings_not_permitted', "Incorrectly specified \"points\" keyword value."),
    ]
    )

def test_great_arc_wrongly_formatted_points(points, expected_error, aia171_test_map):
    coordinate_frame = aia171_test_map.coordinate_frame
    a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=coordinate_frame)
    b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=coordinate_frame)
    with pytest.raises(ValueError,match=re.escape(expected_error)):
        GreatArc(a, b, points=points)

    with pytest.raises(ValueError, match=re.escape(expected_error)):
        GreatArc(a, b).coordinates(points=points)

    with pytest.raises(ValueError, match=re.escape(expected_error)):
        GreatArc(a, b).inner_angles(points=points)

    with pytest.raises(ValueError, match=re.escape(expected_error)):
        GreatArc(a, b).distances(points=points)

    with pytest.raises(ValueError, match=re.escape(expected_error)):
        GreatArc(a, b).distances(points=points)


# Test that the great arc code properly differentiates between the default
# points and the requested points
def test_great_arc_points_differentiates(aia171_test_map):
    coordinate_frame = aia171_test_map.coordinate_frame
    a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=coordinate_frame)
    b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=coordinate_frame)
    gc = GreatArc(a, b)
    coordinates = gc.coordinates(10)
    inner_angles = gc.inner_angles(11)
    distances = gc.distances(12)
    assert len(coordinates) == 10
    assert len(gc.coordinates()) == 100
    assert len(inner_angles) == 11
    assert len(gc.inner_angles()) == 100
    assert len(distances) == 12
    assert len(gc.distances()) == 100


# Test that the great arc code properly understands different observers
# for the start and end points
def test_great_arc_different_observer(aia171_test_map):
    a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=aia171_test_map.coordinate_frame)

    observer = SkyCoord(-10.0*u.deg, 83*u.deg, radius=0.9*u.au,
                        frame=frames.HeliographicStonyhurst, obstime=aia171_test_map.date)
    b = SkyCoord(400*u.arcsec, 600*u.arcsec, observer=observer, frame=frames.Helioprojective)

    # Test that the input observers are indeed different
    assert a.observer.lon != b.observer.lon
    assert a.observer.lat != b.observer.lat
    assert a.observer.radius != b.observer.radius

    # Create the great arc
    gc = GreatArc(a, b)

    # The start and end points stored internally are Heliocentric
    start = gc.start
    assert isinstance(start.frame, frames.Heliocentric)
    end = gc.end
    assert isinstance(end.frame, frames.Heliocentric)

    # The start and end points stored internally have the same observer
    assert start.observer.lon == end.observer.lon
    assert start.observer.lat == end.observer.lat
    assert start.observer.radius == end.observer.radius

    # The start point stored internally has the Heliocentric coordinates of the initial coordinate passed in.
    a2h = a.transform_to(frames.Heliocentric)
    assert start.x == a2h.x
    assert start.y == a2h.y
    assert start.z == a2h.z

    # The end point stored internally has the Heliocentric coordinates of the initial coordinate passed in.
    b2h = b.transform_to(frames.Heliocentric(observer=aia171_test_map.observer_coordinate))

    # Missing an dp on b2h compared to end (TODO BUG?)
    np.testing.assert_almost_equal(end.x.value, b2h.x.value)
    np.testing.assert_almost_equal(end.y.value, b2h.y.value)
    np.testing.assert_almost_equal(end.z.value, b2h.z.value)


@pytest.fixture
def rectangle_args():
    bottom_left = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame='heliographic_stonyhurst')
    top_right = SkyCoord(10 * u.arcsec, 10 * u.arcsec, frame='heliographic_stonyhurst')
    width = 10 * u.arcsec
    height = 10 * u.arcsec

    return bottom_left, top_right, width, height


def test_rectangle_incomplete_input(rectangle_args):
    bottom_left, _, _, height = rectangle_args

    with pytest.raises(ValueError, match="Invalid input, either bottom_left and top_right or bottom_left and height and width should be provided."):
        get_rectangle_coordinates(bottom_left, height=height)


def test_rectangle_invalid_input(rectangle_args):
    _, _, width, height = rectangle_args

    with pytest.raises(TypeError):
        get_rectangle_coordinates(width, height=height)


def test_rectangle_all_parameters_passed(rectangle_args):
    bottom_left, top_right, width, height = rectangle_args

    with pytest.raises(ValueError, match="Invalid input, width, height and top_right parameters should not be passed simultaneously."):
        get_rectangle_coordinates(bottom_left, width=width, top_right=top_right, height=height)


def test_rectangle_width_height(rectangle_args):
    bottom_left, _, width, height = rectangle_args

    _, top_right_1 = get_rectangle_coordinates(bottom_left, width=width, height=height)

    assert bottom_left.spherical.lon + width == top_right_1.spherical.lon
    assert bottom_left.spherical.lat + height == top_right_1.spherical.lat


def test_rectangle_mismatching_frames_missing_parameters(rectangle_args):
    bottom_left, top_right, _, _ = rectangle_args
    top_right = SkyCoord(10 * u.arcsec, 10 * u.arcsec, frame='heliographic_carrington')

    with pytest.raises(ConvertError):
        bottom_left, top_right = get_rectangle_coordinates(bottom_left, top_right=top_right)


def test_rectangle_top_right(rectangle_args):
    bottom_left, top_right, _, _ = rectangle_args

    bottom_left_1, top_right_1 = get_rectangle_coordinates(bottom_left, top_right=top_right)

    assert bottom_left.spherical.lon == bottom_left_1.spherical.lon
    assert bottom_left.spherical.lat == bottom_left_1.spherical.lat
    assert top_right.spherical.lon == top_right_1.spherical.lon
    assert top_right.spherical.lat == top_right_1.spherical.lat


def test_rectangle_bottom_left_different_types(rectangle_args):
    bottom_left, _, width, height = rectangle_args

    bottom_left_1, top_right_1 = get_rectangle_coordinates(
        bottom_left.frame, width=width, height=height)

    assert bottom_left.spherical.lon + width == top_right_1.spherical.lon
    assert bottom_left.spherical.lat + height == top_right_1.spherical.lat
    assert type(bottom_left_1) == type(top_right_1) == type(bottom_left.frame)  # NOQA: E721

    bottom_left_2, top_right_2 = get_rectangle_coordinates(bottom_left, width=width, height=height)

    assert bottom_left.spherical.lon + width == top_right_2.spherical.lon
    assert bottom_left.spherical.lat + height == top_right_2.spherical.lat
    assert type(bottom_left_2) == type(top_right_2) == type(bottom_left)   # NOQA: E721


def test_rectangle_bottom_left_vector():
    bottom_left_vector = SkyCoord([0 * u.arcsec, 10 * u.arcsec], [0 * u.arcsec, 10 * u.arcsec],
                                  frame='heliographic_stonyhurst')

    bottom_left, top_right = get_rectangle_coordinates(bottom_left_vector)

    assert bottom_left.spherical.lon == bottom_left_vector[0].spherical.lon
    assert bottom_left.spherical.lat == bottom_left_vector[0].spherical.lat
    assert top_right.spherical.lon == bottom_left_vector[1].spherical.lon
    assert top_right.spherical.lat == bottom_left_vector[1].spherical.lat


def test_solar_angle_equivalency_inputs():

    with pytest.raises(TypeError):
        solar_angle_equivalency("earth")

    test_coord = SkyCoord(0*u.arcsec, 0*u.arcsec)
    with pytest.raises(ValueError, match="Observer must have an observation time, `obstime`."):
        solar_angle_equivalency(test_coord)


def test_solar_angle_equivalency_outputs():
    observer = get_earth("2020-11-16")

    distance_in_arcsec = 1*u.arcsec
    distance_in_km = distance_in_arcsec.to(u.km, equivalencies=solar_angle_equivalency(observer))
    distance_back_to_arcsec = distance_in_km.to(u.arcsec, equivalencies=solar_angle_equivalency(observer))

    assert_quantity_allclose(distance_in_km, 717.25668*u.km)
    assert_quantity_allclose(distance_back_to_arcsec, distance_in_arcsec)


def test_limb_coords():
    observer = SkyCoord(0*u.deg, 0*u.deg, 1*u.au,
                        obstime='2021-01-01',
                        frame='heliographic_stonyhurst')

    limb_coords = get_limb_coordinates(observer)
    assert isinstance(limb_coords, SkyCoord)
    assert limb_coords.obstime == observer.obstime
    assert u.allclose(limb_coords.heliographic_stonyhurst.radius,
                      constants.radius)

    resolution = 2000
    limb_coords = get_limb_coordinates(observer, resolution=resolution)
    assert len(limb_coords) == resolution

    rsun = 1.1 * constants.radius
    limb_coords = get_limb_coordinates(observer, rsun=rsun)
    assert u.allclose(limb_coords.heliographic_stonyhurst.radius, rsun)

    with pytest.raises(ValueError, match='Observer distance must be greater than rsun'):
        get_limb_coordinates(observer, 1.1 * u.au)


# Tests of the helper functions
def test_get_new_observer(aia171_test_map):
    initial_obstime = aia171_test_map.date
    rotation_interval = 2 * u.day
    new_time = initial_obstime + rotation_interval
    time_delta = new_time - initial_obstime
    observer = get_earth(initial_obstime + rotation_interval)

    # The observer time is set along with other definitions of time
    for time in (rotation_interval, new_time, time_delta):
        with pytest.raises(ValueError, match="Either the 'observer' or the 'time' keyword must be specified, but not both simultaneously."):
            new_observer = get_new_observer(initial_obstime, observer, time)

    # Obstime property is present but the value is None
    observer_obstime_is_none = SkyCoord(12*u.deg, 46*u.deg, frame=frames.HeliographicStonyhurst)
    with pytest.raises(ValueError, match="The observer 'obstime' property must not be None."):
        new_observer = get_new_observer(None, observer_obstime_is_none, None)

    # When the observer is set, it gets passed back out
    new_observer = get_new_observer(initial_obstime, observer, None)
    assert isinstance(new_observer, SkyCoord)
    np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).lon.to(u.deg).value,
                                   observer.transform_to(frames.HeliographicStonyhurst).lon.to(u.deg).value, decimal=3)
    np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).lat.to(u.deg).value,
                                   observer.transform_to(frames.HeliographicStonyhurst).lat.to(u.deg).value, decimal=3)
    np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).radius.to(u.au).value,
                                   observer.transform_to(frames.HeliographicStonyhurst).radius.to(u.au).value, decimal=3)

    # When the time is set, a coordinate for Earth comes back out
    for time in (rotation_interval, new_time, time_delta):
        with pytest.warns(UserWarning, match="Using 'time' assumes an Earth-based observer"):
            new_observer = get_new_observer(initial_obstime, None, time)
        assert isinstance(new_observer, SkyCoord)

        np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).lon.to(u.deg).value,
                                       observer.transform_to(frames.HeliographicStonyhurst).lon.to(u.deg).value, decimal=3)
        np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).lat.to(u.deg).value,
                                       observer.transform_to(frames.HeliographicStonyhurst).lat.to(u.deg).value, decimal=3)
        np.testing.assert_almost_equal(new_observer.transform_to(frames.HeliographicStonyhurst).radius.to(u.au).value,
                                       observer.transform_to(frames.HeliographicStonyhurst).radius.to(u.au).value, decimal=3)

    # The observer and the time cannot both be None
    with pytest.raises(ValueError, match="Either the 'observer' or the 'time' keyword must not be None."):
        new_observer = get_new_observer(initial_obstime, None, None)


# Please note the numbers in these tests are not checked for physical
# accuracy, only that they are the values the function was outputting upon
# implementation. This is not a significant issue for the diff_rot function
# since it is relatively simple and the values it produces can be easily
# compared to other implementations of the same simple function. The same
# cannot be said for the solar_coordinate_rotation function. This functionality
# relies accurate knowledge of the solar ephemeris in particular.
# There is no reference implementation of the solar_coordinate_rotation function
# of demonstrated trustworthiness at time of writing in any language. There
# are no known independent values or tests that can be used to test the
# veracity of the solar_coordinate_rotation function. This being the case, the
# solar_coordinate_rotation function is tested against values that it generated.
# Therefore these tests test for consistency, not accuracy. Note that when the
# 0.8.0 branch was released, the solar ephemeris calculation was handed off to
# the relevant Astropy code. The solar_coordinate_rotation tests were changed
# for self-consistency. Note that the change in position comparing the results
# of pre- and 0.8.0 sunpy solar coordinate rotation functionality (rot_hpc
# and solar_coordinate_rotation respectively) was on the order of 0.5 arcseconds.
# At time of writing, the difference between the rotation
# calculated using the pre-0.8.0 rot_hpc function and the SSWIDL equivalent
# rot_xy.pro for the tests given in pre-0.8.0 were on the order of hundredths
# of an arcsecond. I suspect that the reason for the small differences is
# because the sunpy's ephemeris and coordinate transformation infrastructure
# was largely based on that in SSWIDL.

def test_solar_coordinate_rotation():
    # Testing along the Sun-Earth line, observer is on the Earth
    obs_time = '2010-09-10 12:34:56'
    observer = get_earth(obs_time)
    c = SkyCoord(-570*u.arcsec, 120*u.arcsec, obstime=obs_time,
                 observer=observer, frame=frames.Helioprojective)
    new_time = '2010-09-11 12:34:56'
    new_observer = get_earth(new_time)

    # Test that when both the observer and the time are specified, an error is raised.
    with pytest.raises(ValueError, match="Either the 'observer' or the 'time' keyword must be specified, but not both simultaneously."):
        d = solar_coordinate_rotation(c, observer=observer, time=new_time)

    # Test that the code properly filters the observer keyword
    with pytest.raises(ValueError, match="The 'observer' must be an astropy.coordinates.BaseCoordinateFrame or an astropy.coordinates.SkyCoord."):
        d = solar_coordinate_rotation(c, observer='earth')

    # Test that the code properly filters the time keyword
    with pytest.raises(ValueError, match="Input values did not match any of the formats where the format keyword is optional"):
        with pytest.warns(UserWarning, match="Using 'time' assumes an Earth-based observer"):
            d = solar_coordinate_rotation(c, time='noon')

    # Test that the code gives the same output for multiple different inputs
    # that define the same observer location and time.
    for i, definition in enumerate((1 * u.day, TimeDelta(1*u.day), new_time, new_observer)):
        if i in (0, 1, 2):
            with pytest.warns(UserWarning, match="Using 'time' assumes an Earth-based observer"):
                d = solar_coordinate_rotation(c, time=definition)
        else:
            d = solar_coordinate_rotation(c, observer=definition)

        # Test that a SkyCoordinate is created
        assert isinstance(d, SkyCoord)

        # Test the coordinate
        np.testing.assert_almost_equal(d.Tx.to(u.arcsec).value, -386.4519332773052, decimal=1)
        np.testing.assert_almost_equal(d.Ty.to(u.arcsec).value, 106.1647811048218, decimal=1)
        np.testing.assert_allclose(d.distance.to(u.km).value, 1.499689e+08, rtol=1e-5)

        # Test that the SkyCoordinate is Helioprojective
        assert isinstance(d.frame, frames.Helioprojective)

    # Test that the function works correctly with a HGS coordinate.
    earth_coord = get_earth(Time("2022-03-30"))
    coord_hpc = SkyCoord(100*u.arcsec, 100*u.arcsec, frame=frames.Helioprojective(observer=earth_coord))

    coord_hgs = coord_hpc.transform_to(frames.HeliographicStonyhurst)
    with pytest.warns(UserWarning, match="Using 'time' assumes an Earth-based observer"):
        rotated_coord_hgs = solar_coordinate_rotation(coord_hgs, time=Time("2022-03-31"))

    assert isinstance(rotated_coord_hgs.frame, frames.HeliographicStonyhurst)


def test_consistency_with_rotatedsunframe():
    old_observer = frames.HeliographicStonyhurst(10*u.deg, 20*u.deg, 1*u.AU, obstime='2001-01-01')
    new_observer = frames.HeliographicStonyhurst(30*u.deg, 40*u.deg, 2*u.AU, obstime='2001-01-08')

    hpc_coord = SkyCoord(100*u.arcsec, 200*u.arcsec, frame='helioprojective',
                         observer=old_observer, obstime=old_observer.obstime)

    # Perform the differential rotation using solar_coordinate_rotation()
    result1 = solar_coordinate_rotation(hpc_coord, observer=new_observer)

    # Perform the differential rotation using RotatedSunFrame, with translational motion of the Sun
    # ignored using transform_with_sun_center()
    rsf_coord = RotatedSunFrame(base=hpc_coord, rotated_time=new_observer.obstime)
    with transform_with_sun_center():
        result2 = rsf_coord.transform_to(result1.replicate_without_data())

    assert_quantity_allclose(result1.Tx, result2.Tx)
    assert_quantity_allclose(result1.Ty, result2.Ty)
    assert_quantity_allclose(result1.distance, result2.distance)
