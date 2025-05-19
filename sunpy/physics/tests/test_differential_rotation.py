
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.coordinates import frames
from sunpy.coordinates.ephemeris import get_earth
from sunpy.physics.differential_rotation import diff_rot, differential_rotate, solar_rotate_coordinate
from sunpy.util.exceptions import SunpyDeprecationWarning

# Please note the numbers in these tests are not checked for physical
# accuracy, only that they are the values the function was outputting upon
# implementation. This is not a significant issue for the diff_rot function
# since it is relatively simple and the values it produces can be easily
# compared to other implementations of the same simple function. The same
# cannot be said for the solar_rotate_coordinate function. This functionality
# relies accurate knowledge of the solar ephemeris in particular.
# There is no reference implementation of the solar_rotate_coordinate function
# of demonstrated trustworthiness at time of writing in any language. There
# are no known independent values or tests that can be used to test the
# veracity of the solar_rotate_coordinate function. This being the case, the
# solar_rotate_coordinate function is tested against values that it generated.
# Therefore these tests test for consistency, not accuracy. Note that when the
# 0.8.0 branch was released, the solar ephemeris calculation was handed off to
# the relevant Astropy code. The solar_rotate_coordinate tests were changed
# for self-consistency. Note that the change in position comparing the results
# of pre- and 0.8.0 sunpy solar coordinate rotation functionality (rot_hpc
# and solar_rotate_coordinate respectively) was on the order of 0.5 arcseconds.
# At time of writing, the difference between the rotation
# calculated using the pre-0.8.0 rot_hpc function and the SSWIDL equivalent
# rot_xy.pro for the tests given in pre-0.8.0 were on the order of hundredths
# of an arcsecond. I suspect that the reason for the small differences is
# because the sunpy's ephemeris and coordinate transformation infrastructure
# was largely based on that in SSWIDL.


@pytest.fixture
def seconds_per_day():
    return 24 * 60 * 60.0 * u.s


def test_diff_rot_deprecated_warning(seconds_per_day):
    with pytest.warns(SunpyDeprecationWarning, match='The diff_rot function is deprecated'):
        diff_rot(10 * seconds_per_day, 30 * u.deg)


def test_solar_rotate_coordinate_deprecated_warning():
    old_observer = frames.HeliographicStonyhurst(10*u.deg, 20*u.deg, 1*u.AU, obstime='2001-01-01')
    new_observer = frames.HeliographicStonyhurst(30*u.deg, 40*u.deg, 2*u.AU, obstime='2001-01-08')

    hpc_coord = SkyCoord(100*u.arcsec, 200*u.arcsec, frame='helioprojective',
                         observer=old_observer, obstime=old_observer.obstime)

    with pytest.warns(SunpyDeprecationWarning, match='The solar_rotate_coordinate function is deprecated'):
        solar_rotate_coordinate(hpc_coord, observer=new_observer)

def test_differential_rotate_deprecated_warning(aia171_test_map, seconds_per_day):
    new_observer = get_earth(aia171_test_map.date + 6*u.hr)
    with pytest.warns(SunpyDeprecationWarning, match='The differential_rotate function is deprecated'):
        differential_rotate(aia171_test_map, observer=new_observer)
