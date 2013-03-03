"""
ANA is a simple script that allows people to access compressed ana files.
It accesses a C library, based on Michiel van Noort's
IDL DLM library 'f0' which contains a cleaned up version of the original
anarw routines.

Created by Tim van Werkhoven (t.i.m.vanwerkhoven@gmail.com) on 2009-02-11.
Copyright (c) 2009--2011 Tim van Werkhoven. All rights reserved.   
"""

#To read a file:
#> anadata = pyana.fzread(<filename>, [debug=0])
#which will return a dict with the data in anadata['data'] and some meta info
#in anadata['header']. To return only the data or header, use pyana.getdata()
#and pyana.getheader() respectively.
#
#To write a file:
#> pyana.fzwrite(<filename>, <data>, [compress=1, [comments=False, [debug=0]]]):
#or use pyana.writeto(), which is an alias to fzwrite().
 
from __future__ import absolute_import

import os

from . import _pyana

__all__ = ['fz_read', 'get_data', 'get_header', 'fz_write', 'write_to']

def fz_read(filename, debug=0):
	"""
	Load an ANA file and return the data, size, dimensions and comments in a
	dict.
	
	data = pyana.load(filename)
	"""
	if not os.path.isfile(filename):
		raise IOError("File does not exist!")
	
	data = _pyana.fzread(filename, debug)
	return data


def get_data(filename, debug=0):
	"""
	Load an ANA file and only return the data as a numpy array.
	
	data = pyana.getdata(filename)
	"""
	return (fz_read(filename, debug))['data']
	# data = _pyana.fzread(filename, debug)
	# return data['data']

	
def get_header(filename, debug=0):
	"""
	Load an ANA file and only return the header consisting of the dimensions,
	size (defined as the product of all dimensions times the size of the
	datatype, this not relying on actual filesize) and comments.
	
	header = pyana.getheader(filename)
	"""
	return (fz_read(filename, debug))['header']
	# data = _pyana.fzread(filename, debug)
	# return data['header']

## Functions for storing ANA files
def fz_write(filename, data, compress=1, comments=False, debug=0):
	"""
	Save a 2d numpy array as an ANA file and return the bytes written, or NULL
	
	written = pyana.fzwrite(filename, data, compress=1, comments=False)
	"""
	if (comments):
		return _pyana.fzwrite(filename, data, compress, comments, debug)
	else:
		return _pyana.fzwrite(filename, data, compress, '', debug)


def write_to(filename, data, compress=1, comments=False, debug=0):
	"""
	Similar as pyana.fzwrite().
	"""
	return fz_write(filename, data, compress, comments, debug)