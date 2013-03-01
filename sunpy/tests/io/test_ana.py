"""
General ANA Tests
"""
import unittest
from sunpy.io import ana
            
class anaTests(unittest.TestCase):
    """Basic ANA tests"""
    def setUp(self):
            # Create a test image, store it, reread it and compare
            import numpy as N
            self.numpy = N
            self.img_size = (456, 345)
            self.img_src = N.arange(N.product(self.img_size))
            self.img_src.shape = self.img_size
            self.img_i8 = self.img_src*2**8/self.img_src.max()
            self.img_i8 = self.img_i8.astype(N.int8)
            self.img_i16 = self.img_src*2**16/self.img_src.max()
            self.img_i16 = self.img_i16.astype(N.int16)
            self.img_f32 = self.img_src*1.0/self.img_src.max()
            self.img_f32 = self.img_f32.astype(N.float32)
	
    def runTests(self):
		unittest.TextTestRunner(verbosity=2).run(self.suite())
	
    def suite(self):
		suite = unittest.TestLoader().loadTestsFromTestCase(anaTests)
		return suite
	
    def testi8c(self):
		# Test int 8 compressed functions
		ana.fz_write('/tmp/pyana-testi8c', self.img_i8, 1, 'testcase', 1)
		self.img_i8c_rec = ana.fz_read('/tmp/pyana-testi8c', 1)
		self.assert_(self.numpy.sum(self.img_i8c_rec['data'] - self.img_i8) == 0,
		 	msg="Storing 8 bits integer data with compression failed (diff: %d)" % (self.numpy.sum(self.img_i8c_rec['data'] - self.img_i8)))
	
    def testi8u(self):
		# Test int 8 uncompressed functions
		ana.fz_write('/tmp/pyana-testi8u', self.img_i8, 0, 'testcase', 1)
		self.img_i8u_rec = ana.fz_read('/tmp/pyana-testi8u', 1)
		self.assert_(self.numpy.sum(self.img_i8u_rec['data'] - self.img_i8) == 0,
			msg="Storing 8 bits integer data without compression failed (diff: %d)" % (self.numpy.sum(self.img_i8u_rec['data'] - self.img_i8)))
	
    def testi16c(self):
		# Test int 16 compressed functions
		ana.fz_write('/tmp/pyana-testi16c', self.img_i16, 1, 'testcase', 1)
		self.img_i16c_rec = ana.fz_read('/tmp/pyana-testi16c', 1)
		self.assert_(self.numpy.allclose(self.img_i16c_rec['data'], self.img_i16),
			msg="Storing 16 bits integer data with compression failed (diff: %d)" % (self.numpy.sum(self.img_i16c_rec['data'] - self.img_i16)))
	
    def testi16u(self):
		# Test int 16 uncompressed functions
		ana.fz_write('/tmp/pyana-testi16u', self.img_i16, 0, 'testcase', 1)
		self.img_i16u_rec = ana.fz_read('/tmp/pyana-testi16u', 1)
		self.assert_(self.numpy.allclose(self.img_i16u_rec['data'], self.img_i16),
			msg="Storing 16 bits integer data without compression failed (diff: %d)" % (self.numpy.sum(self.img_i16u_rec['data'] - self.img_i16)))
	
    def testf32u(self):
		# Test float 32 uncompressed functions
		ana.fz_write('/tmp/pyana-testf32', self.img_f32, 0, 'testcase', 1)
		self.img_f32_rec = ana.fz_read('/tmp/pyana-testf32', 1)
		self.assert_(self.numpy.allclose(self.img_f32_rec['data'], self.img_f32),
			msg="Storing 32 bits float data without compression failed (diff: %g)" % (1.0*self.numpy.sum(self.img_f32_rec['data'] - self.img_f32)))
	
    def testf32c(self):
		# Test if float 32 compressed fails
		self.assertRaises(RuntimeError, ana.fz_write, '/tmp/pyana-testf32', self.img_f32, 1, 'testcase', 1)
		

#if __name__ == "__main__":	
#	suite = unittest.TestLoader().loadTestsFromTestCase(anaTests)
#	unittest.TextTestRunner(verbosity=2).run(suite)