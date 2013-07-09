from __future__ import absolute_import
import pytest
import numpy as np
from sunpy.coords import diff_rot, sun_pos, calc_P_B0_SD, rot_hpc
#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103


class TestDiffRot():
    """ Please note the numbers in these tests are not checked for physical
    accuracy, only that they are the values the function was outputting upon
    implementation. """

    def test_single(self):
        rot = diff_rot(10, 30)
        np.testing.assert_almost_equal(rot,136.8216, decimal=4)

    def test_array(self):
        rot = diff_rot(10, np.linspace(-70, 70, 2))
        np.testing.assert_almost_equal(rot, np.array([110.2725, 110.2725]),
                                       decimal=4)

    def test_synodic(self):
        rot = diff_rot(10, 30, rot_type='howard', frame_time='synodic')
        np.testing.assert_almost_equal(rot, 126.9656, decimal=4)

    def test_sidereal(self):
        rot = diff_rot(10, 30, rot_type='howard', frame_time='sidereal')
        np.testing.assert_almost_equal(rot, 136.8216, decimal=4)

    def test_howard(self):
        rot = diff_rot(10, 30, rot_type='howard')
        np.testing.assert_almost_equal(rot, 136.8216, decimal=4)

    def test_allen(self):
        rot = diff_rot(10, 30, rot_type='allen')
        np.testing.assert_almost_equal(rot, 136.9, decimal=4)

    def test_snodgrass(self):
        rot = diff_rot(10, 30, rot_type='snodgrass')
        np.testing.assert_almost_equal(rot,135.4232, decimal=4)

    def test_fail(self):
        pytest.raises(ValueError, diff_rot, 10, 30, rot_type='garbage')

    def test_sunpos(self):
        result = sun_pos('2013-05-14')
        assertion = {'obliq': 23.4358, 'app_long': 53.3683, 'dec': 18.6125,
                     'ra': 50.9796, 'longitude': 53.3705}
        for k in assertion:
            np.testing.assert_almost_equal(result[k], assertion[k], decimal=4)

    def test_pb0r(self):
        result = calc_P_B0_SD('2012-12-14')
        assertion = {'p': 10.4868, 'b0': -0.8127, 'l0': 0.0000, 'sd': 16.2364}
        for k in assertion:
            np.testing.assert_almost_equal(result[k], assertion[k], decimal=4)

    def test_rot_hpc(self):
        # testing along the Sun-Earth line, observer is on the Earth
        x, y = rot_hpc(451.4, -108.9, '2012-06-15', '2012-06-15 16:05:23')
        np.testing.assert_almost_equal(x, 574.2, decimal=1)
        np.testing.assert_almost_equal(y, -108.4, decimal=1)