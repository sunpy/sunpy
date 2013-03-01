"""
ana.py

A C extension for Python to read ana f0 files. Based on Michiel van Noort's
IDL DLM library 'f0' which contains a cleaned up version of the original
anarw routines.

To read a file:
> anadata = pyana.fzread(<filename>, [debug=0])
which will return a dict with the data in anadata['data'] and some meta info
in anadata['header']. To return only the data or header, use pyana.getdata()
and pyana.getheader() respectively.

To write a file:
> pyana.fzwrite(<filename>, <data>, [compress=1, [comments=False, [debug=0]]]):
or use pyana.writeto(), which is an alias to fzwrite().

Created by Tim van Werkhoven (t.i.m.vanwerkhoven@gmail.com) on 2009-02-11.
Copyright (c) 2009--2011 Tim van Werkhoven. All rights reserved.
"""

from __future__ import absolute_import

import os
import unittest
import sunpy.io.ana._pyana as _pyana

## Functions for loading in ANA files

def fzread(filename, debug=0):
	"""
	Load an ANA file and return the data, size, dimensions and comments in a
	dict.
	
	data = pyana.load(filename)
	"""
	if not os.path.isfile(filename):
		raise IOError("File does not exist!")
	
	data = _pyana.fzread(filename, debug)
	return data


def getdata(filename, debug=0):
	"""
	Load an ANA file and only return the data as a numpy array.
	
	data = pyana.getdata(filename)
	"""
	return (fzread(filename, debug))['data']
	# data = _pyana.fzread(filename, debug)
	# return data['data']

	
def getheader(filename, debug=0):
	"""
	Load an ANA file and only return the header consisting of the dimensions,
	size (defined as the product of all dimensions times the size of the
	datatype, this not relying on actual filesize) and comments.
	
	header = pyana.getheader(filename)
	"""
	return (fzread(filename, debug))['header']
	# data = _pyana.fzread(filename, debug)
	# return data['header']

## Functions for storing ANA files
def fzwrite(filename, data, compress=1, comments=False, debug=0):
	"""
	Save a 2d numpy array as an ANA file and return the bytes written, or NULL
	
	written = pyana.fzwrite(filename, data, compress=1, comments=False)
	"""
	if (comments):
		return _pyana.fzwrite(filename, data, compress, comments, debug)
	else:
		return _pyana.fzwrite(filename, data, compress, '', debug)


def writeto(filename, data, compress=1, comments=False, debug=0):
	"""
	Similar as pyana.fzwrite().
	"""
	return fzwrite(filename, data, compress, comments, debug)


## Selftesting using unittest starts below this line
class pyanaTests(unittest.TestCase):
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
		suite = unittest.TestLoader().loadTestsFromTestCase(pyanaTests)
		return suite
	
	def testi8c(self):
		# Test int 8 compressed functions
		fzwrite('/tmp/pyana-testi8c', self.img_i8, 1, 'testcase', 1)
		self.img_i8c_rec = fzread('/tmp/pyana-testi8c', 1)
		self.assert_(self.numpy.sum(self.img_i8c_rec['data'] - self.img_i8) == 0,
		 	msg="Storing 8 bits integer data with compression failed (diff: %d)" % (self.numpy.sum(self.img_i8c_rec['data'] - self.img_i8)))
	
	def testi8u(self):
		# Test int 8 uncompressed functions
		fzwrite('/tmp/pyana-testi8u', self.img_i8, 0, 'testcase', 1)
		self.img_i8u_rec = fzread('/tmp/pyana-testi8u', 1)
		self.assert_(self.numpy.sum(self.img_i8u_rec['data'] - self.img_i8) == 0,
			msg="Storing 8 bits integer data without compression failed (diff: %d)" % (self.numpy.sum(self.img_i8u_rec['data'] - self.img_i8)))
	
	def testi16c(self):
		# Test int 16 compressed functions
		fzwrite('/tmp/pyana-testi16c', self.img_i16, 1, 'testcase', 1)
		self.img_i16c_rec = fzread('/tmp/pyana-testi16c', 1)
		self.assert_(self.numpy.allclose(self.img_i16c_rec['data'], self.img_i16),
			msg="Storing 16 bits integer data with compression failed (diff: %d)" % (self.numpy.sum(self.img_i16c_rec['data'] - self.img_i16)))
	
	def testi16u(self):
		# Test int 16 uncompressed functions
		fzwrite('/tmp/pyana-testi16u', self.img_i16, 0, 'testcase', 1)
		self.img_i16u_rec = fzread('/tmp/pyana-testi16u', 1)
		self.assert_(self.numpy.allclose(self.img_i16u_rec['data'], self.img_i16),
			msg="Storing 16 bits integer data without compression failed (diff: %d)" % (self.numpy.sum(self.img_i16u_rec['data'] - self.img_i16)))
	
	def testf32u(self):
		# Test float 32 uncompressed functions
		fzwrite('/tmp/pyana-testf32', self.img_f32, 0, 'testcase', 1)
		self.img_f32_rec = fzread('/tmp/pyana-testf32', 1)
		self.assert_(self.numpy.allclose(self.img_f32_rec['data'], self.img_f32),
			msg="Storing 32 bits float data without compression failed (diff: %g)" % (1.0*self.numpy.sum(self.img_f32_rec['data'] - self.img_f32)))
	
	def testf32c(self):
		# Test if float 32 compressed fails
		self.assertRaises(RuntimeError, fzwrite, '/tmp/pyana-testf32', self.img_f32, 1, 'testcase', 1)
		

if __name__ == "__main__":	
	suite = unittest.TestLoader().loadTestsFromTestCase(pyanaTests)
	unittest.TextTestRunner(verbosity=2).run(suite)
	

