from numpy.testing import assert_array_almost_equal as a
from sunpy.physics.transforms.differential_rotation import diff_rot
from astropy import units as u
from astropy.coordinates import Latitude

def test_diff_rot():
    """ Please note the numbers in these tests are not checked for physical
    accuracy, only that they are the values the function was outputting upon
    implementation. """

    a(diff_rot(10, Latitude(30, u.deg)), 136.8216 * u.deg, decimal = 4)
    a(diff_rot(10, Latitude(30, u.deg), rot_type='howard', frame_time='synodic'), 126.9656 * u.deg,decimal = 4)
    a(diff_rot(10, Latitude(30, u.deg), rot_type='howard', frame_time='sidereal'), 136.8216 * u.deg, decimal=4)				
    a(diff_rot(10, Latitude(30, u.deg), rot_type='howard'), 136.8216 * u.deg, decimal =4)
    a(diff_rot(10, Latitude(30, u.deg), rot_type='allen'), 136.9 * u.deg, decimal = 1)
    a(diff_rot(10, Latitude(30, u.deg), rot_type='snodgrass'), 135.4232 * u.deg, decimal = 4)

