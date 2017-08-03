from __future__ import absolute_import, division, print_function

import pytest

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import sunpy.coordinates
import sunpy.map
import sunpy.data.test
from sunpy.coordinates import frames
from sunpy.coordinates.utils import GreatArc
from sunpy.sun import sun


# Test the great arc code against calculable quantities
# The inner angle is the same between each pair of co-ordinates.  You can
# calculate these co-ordinates using the inner angle formulae as listed here:
# https://en.wikipedia.org/wiki/Great-circle_distance
#
@pytest.mark.parametrize("start, end", [((0, 0), (0, 45)),
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

    assert gc.start == c
    assert gc.end == d
    np.testing.assert_almost_equal(gc.inner_angle.to('deg').value, 45.0)
    np.testing.assert_almost_equal(gc.radius.to('km').value, sun.constants.radius.to('km').value)
    np.testing.assert_almost_equal(gc.distance.to('km').value, sun.constants.radius.to('km').value * 2 * np.pi/8, decimal=1)


# Test the calculation of coordinates using varying numbers of points on
# initialization of the GreatArc object.
@pytest.mark.parametrize("points_requested, points_expected, first_point, last_point, last_inner_angle, last_distance",
                         # Test default
                         [(None, 100, (600, -600), (-100, 800), 1.8683580432741789, 1300377.1981272686),
                          # Test int as an option
                          (3, 3, (600, -600), (-100, 800), 1.8683580432741789, 1300377.1981272686),
                          # Test equally spaced monotonically increasing numpy
                          # array
                          (np.linspace(0, 1, 43), 43, (600, -600), (-100, 800), 1.8683580432741789, 1300377.1981272686),
                          # Test unequally spaced monotonically increasing numpy
                          # array
                          (np.asarray([0.1, 0.2, 0.6, 0.67, 0.99]), 5, (604.68091703, -468.64217597), (-88.83212616, 792.76284375), 1.84967446, 1287373.426146),
                          # Test unequally spaced monotonically decreasing numpy
                          # array
                          (np.asarray([0.93, 0.78, 0.3, 0.001]), 4, (-21.28208654, 743.58866798), (600.1512768 , -598.78376614), 0.00186836, 1300.37719813),
                          # Test numpy array that increases and decreases
                          (np.asarray([0.94, 0.73, 0.8, 0.21]), 4, (-32.5852606 , 752.45507707), (585.45829119, -305.26965043), 0.39235519, 273079.21160673)])
def test_great_arc_coordinates(points_requested, points_expected, first_point,
                               last_point, last_inner_angle, last_distance):

    m = sunpy.map.Map(sunpy.map.Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits')))
    coordinate_frame = m.coordinate_frame
    a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=coordinate_frame)
    b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=coordinate_frame)
    gc = GreatArc(a, b, points=points_requested)
    coordinates = gc.coordinates()
    inner_angles = gc.inner_angles()
    distances = gc.distances()

    # Ensure a GreatArc object is returned
    assert isinstance(gc, GreatArc)

    # Test the properties of the GreatArc object
    assert gc.start == a
    assert gc.end == b
    assert gc.distance_unit == u.km
    assert gc.observer == a.observer
    assert gc.center.x == 0 * u.km
    assert gc.center.y == 0 * u.km
    assert gc.center.z == 0 * u.km

    np.testing.assert_almost_equal(gc.start_cartesian, np.asarray([428721.09135385, -428722.90519236, 341776.09102756]))
    np.testing.assert_almost_equal(gc.end_cartesian, np.asarray([-71429.52293813, 571439.07124801, 390859.57978693]))
    np.testing.assert_almost_equal(gc.center_cartesian, np.asarray([0, 0, 0]))

    np.testing.assert_almost_equal(gc.v1, np.asarray([428721.09135385, -428722.90519236, 341776.09102756]))
    np.testing.assert_almost_equal(gc._r, 696000.00000451738)
    np.testing.assert_almost_equal(gc.v2, np.asarray([-71429.52293813, 571439.07124801, 390859.57978693]))
    np.testing.assert_almost_equal(gc.v3, np.asarray([56761.62657985, 466230.70058902, 513637.08158833]))

    # Inner angle
    assert gc.inner_angle.unit == u.rad
    np.testing.assert_almost_equal(gc.inner_angle.value, 1.8683580432741789)

    # Distance
    assert gc.distance.unit == u.km
    np.testing.assert_almost_equal(gc.distance.value, 1300377.1981272686)

    # Radius of the sphere
    assert gc.radius.unit == u.km
    np.testing.assert_almost_equal(gc.radius.value, 696000.0000045174)

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
    np.testing.assert_almost_equal(distances[-1].value, last_distance)


# Test that the great arc code rejects wrongly formatted points
@pytest.mark.parametrize("points", [np.asarray([[0, 0.1], [0.2, 0.3]]),
                                    np.asarray([0.1, 0.2, -0.1, 0.4]),
                                    np.asarray([0.3, 1.1, 0.6, 0.7]),
                                    'strings_not_permitted'])
def test_great_arc_wrongly_formatted_points(points):
    m = sunpy.map.Map(sunpy.map.Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits')))
    coordinate_frame = m.coordinate_frame
    a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=coordinate_frame)
    b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=coordinate_frame)
    with pytest.raises(ValueError):
        dummy = GreatArc(a, b, points=points)

    with pytest.raises(ValueError):
        dummy = GreatArc(a, b).coordinates(points=points)

    with pytest.raises(ValueError):
        dummy = GreatArc(a, b).inner_angles(points=points)

    with pytest.raises(ValueError):
        dummy = GreatArc(a, b).distances(points=points)


# Test that the great arc code properly differentiates between the default
# points and the requested points.
def test_great_arc_points_differentiates():
    m = sunpy.map.Map(sunpy.map.Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits')))
    coordinate_frame = m.coordinate_frame
    a = SkyCoord(600*u.arcsec, -600*u.arcsec, frame=coordinate_frame)
    b = SkyCoord(-100*u.arcsec, 800*u.arcsec, frame=coordinate_frame)
    gc = GreatArc(a, b)
    coordinates = gc.coordinates(10)
    inner_angles = gc.inner_angles(11)
    distances = gc.distances(12)
    assert len(coordinates) == 10 and len(gc.coordinates()) == 100
    assert len(inner_angles) == 11 and len(gc.inner_angles()) == 100
    assert len(distances) == 12 and len(gc.distances()) == 100
