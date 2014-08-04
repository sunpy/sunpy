from __future__ import absolute_import
import unittest
import numpy as np
from astropy import units as u
from astropy.coordinates import Longitude
from sunpy.physics.transforms.differential_rotation import diff_rot
#pylint: disable=C0103,R0904,W0201,W0212,W0232,E1103

class DiffRotTest(unittest.TestCase):
    """ Please note the numbers in these tests are not checked for physical
    accuracy, only that they are the values the function was outputting upon
    implementation. """
    
    def test_single(self):
        rot = diff_rot(10, 30 * u.deg)
        self.failUnless(rot == 136.8216 * u.deg)
        
    def test_array(self):
        rot = diff_rot(10, np.linspace(-70, 70, 2) * u.deg)
        self.failUnless(np.array_equal(rot, Longitude(np.array([110.2725,  110.2725]) * u.deg)))
        
    def test_synodic(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='howard', frame_time='synodic')
        self.failUnless(rot == 126.9656 * u.deg)
        
    def test_sidereal(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='howard', frame_time='sidereal')
        self.failUnless(rot == 136.8216 * u.deg)
        
    def test_howard(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='howard')
        self.failUnless(rot == 136.8216 * u.deg)
 
    def test_allen(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='allen')
        self.failUnless(rot == 136.9 * u.deg)
    
    def test_snodgrass(self):
        rot = diff_rot(10, 30 * u.deg, rot_type='snodgrass')
        self.failUnless(rot == 135.4232 * u.deg)
    
    def test_fail(self):
        self.assertRaises(ValueError, diff_rot, 10, 30 * u.deg, rot_type='garbage')

