"""
General ANA Tests
"""
import unittest
from sunpy.io import ana
import numpy as np
            
class anaTests(unittest.TestCase):
    """Basic ANA tests"""
    def setUp(self):
            # Create a test image, store it, reread it and compare
            self.img_size = (456, 345)
            self.img_src = np.arange(np.product(self.img_size))
            self.img_src.shape = self.img_size
            self.img_i8 = self.img_src*2**8/self.img_src.max()
            self.img_i8 = self.img_i8.astype(np.int8)
            self.img_i16 = self.img_src*2**16/self.img_src.max()
            self.img_i16 = self.img_i16.astype(np.int16)
            self.img_f32 = self.img_src*1.0/self.img_src.max()
            self.img_f32 = self.img_f32.astype(np.float32)
    
    def runTests(self):
        unittest.TextTestRunner(verbosity=2).run(self.suite())
    
    def suite(self):
        suite = unittest.TestLoader().loadTestsFromTestCase(anaTests)
        return suite
    
    def testi8c(self):
        # Test int 8 compressed functions
        ana.write('/tmp/pyana-testi8c', self.img_i8, 'testcase', 0)
        self.img_i8c_rec = ana.read('/tmp/pyana-testi8c')
        self.assert_(np.sum(self.img_i8c_rec[0][0] - self.img_i8) == 0,
             msg="Storing 8 bits integer data with compression failed (diff: %d)" % (np.sum(self.img_i8c_rec[0][0] - self.img_i8)))
    
    def testi8u(self):
        # Test int 8 uncompressed functions
        ana.write('/tmp/pyana-testi8u', self.img_i8, 'testcase', 0)
        self.img_i8u_rec = ana.read('/tmp/pyana-testi8u')
        self.assert_(np.sum(self.img_i8u_rec[0][0] - self.img_i8) == 0,
            msg="Storing 8 bits integer data without compression failed (diff: %d)" % (np.sum(self.img_i8u_rec[0][0] - self.img_i8)))
    
    def testi16c(self):
        # Test int 16 compressed functions
        ana.write('/tmp/pyana-testi16c', self.img_i16, 'testcase', 0)
        self.img_i16c_rec = ana.read('/tmp/pyana-testi16c')
        self.assert_(np.sum(self.img_i16c_rec[0][0] - self.img_i16) == 0,
            msg="Storing 16 bits integer data with compression failed (diff: %d)" % (np.sum(self.img_i16c_rec[0][0] - self.img_i16)))
    
    def testi16u(self):
        # Test int 16 uncompressed functions
        ana.write('/tmp/pyana-testi16u', self.img_i16, 'testcase', 0)
        self.img_i16u_rec = ana.read('/tmp/pyana-testi16u')
        self.assert_(np.sum(self.img_i16u_rec[0][0] - self.img_i16) == 0,
            msg="Storing 16 bits integer data without compression failed (diff: %d)" % (np.sum(self.img_i16u_rec[0][0] - self.img_i16)))
    
    def testf32u(self):
        # Test float 32 uncompressed functions
        ana.write('/tmp/pyana-testf32u', self.img_f32, 'testcase', 0)
        self.img_f32u_rec = ana.read('/tmp/pyana-testf32u')
        self.assert_(np.sum(self.img_f32u_rec[0][0]- self.img_f32) == 0,
            msg="Storing 32 bits float data without compression failed (diff: %g)" % (1.0*np.sum(self.img_f32u_rec[0][0] - self.img_f32)))
    
    def testf32c(self):
        # Test if float 32 compressed functions
        #TODO: Bug with same code. Needs to be tracked down.
    
#        ana.write('/tmp/pyana-testf32c', self.img_f32, 1, 'testcase', 0)
#        self.img_f32c_rec = ana.read('/tmp/pyana-testf32c', 1)
#        self.assert_(np.sum(self.img_f32c_rec[0][1]- self.img_f32) == 0,
#            msg="Storing 32 bits float data without compression failed (diff: %g)" % (1.0*np.sum(self.img_f32c_rec[0][1] - self.img_f32)))    
           self.assertRaises(RuntimeError, ana.write, '/tmp/pyana-testf32c', self.img_f32, 'testcase', 1)