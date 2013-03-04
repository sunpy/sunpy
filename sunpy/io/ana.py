"""
ANA is a script that allows people to access compressed ana files.
It accesses a C library, based on Michiel van Noort's
IDL DLM library 'f0' which contains a cleaned up version of the original
anarw routines.

Created by Tim van Werkhoven (t.i.m.vanwerkhoven@gmail.com) on 2009-02-11.
Copyright (c) 2009--2011 Tim van Werkhoven. All rights reserved.   

Examples
--------

    To read a file:
        anadata = sunpy.io.ana.read(<filename>, [debug=0])
    which will return a dict with the data in anadata['data'] and 
    some meta info in anadata['header']. To return only the data or header, 
    use sunpy.io.ana.getdata() and sunpy.io.ana.getheader() respectively.

    To write a file:
        sunpy.io.ana.write(<filename>, <data>, [compress=1, [comments=False, [debug=0]]]):
"""
 
from __future__ import absolute_import

import os

from . import _pyana

__all__ = ['read', 'get_data', 'get_header', 'write']

def read(filename, debug=0):
    """
    Loads an ANA file and returns the data, size, dimensions and comments in a
    dictionary.
    
    Parameters
    ----------
    filename: string
        Name of file to be read.
    debug: int
        Does something 0 or 1.
    
    Returns
    -------
    out: dictionary
        Contains data, size, diemnsions and comments of an ANA file.    
    
    Examples
    --------
    data = sunpy.io.ana.read(filename)
    
    """
    if not os.path.isfile(filename):
        raise IOError("File does not exist!")
	
    data = _pyana.fzread(filename, debug)
    return data

def get_data(filename, debug=0):
    """
    Loads an ANA file and only returns the data as a numpy array.

    Parameters
    ----------
    filename: string
        Name of file to be read.
    debug: int, optional
        Does something 0 or 1.
    
    Returns
    -------
    out: ndarray
        Contains the data only of an ANA file.    
    
    Examples
    --------
    data = sunpy.io.ana.getdata(filename)
    """
    data = read(filename, debug)
    return data['data']

	
def get_header(filename, debug=0):
    """
    Load an ANA file and only return the header consisting of the dimensions,
    size (defined as the product of all dimensions times the size of the
    datatype, this not relying on actual filesize) and comments.

    Parameters
    ----------
    filename: string
        Name of file to be read.
    debug: int, optional
        Does something 0 or 1.
    
    Returns
    -------
    out: string
        Contains the header only of an ANA file.    

    Examples
    --------    
    header = sunpy.io.ana.getheader(filename)
    """
    data = read(filename, debug)
    return data['header']

def write(filename, data, compress=1, comments=False, debug=0):
    """
    Saves a 2D numpy array as an ANA file and returns the bytes written or NULL

    Parameters
    ----------
    filename: string
        Name of file to be created.
    data: numpy array
        Name of data to be stored.
    compress: int, optional
        To compress the data or not.
        1 is to compress, 0 is uncompressed
    commnets: string, optional
        The comments to be stored as a header.
    debug: int, optional
        Does something 0 or 1.
    
    Returns
    -------
    out: ANA compressed archive
        A new ANA compressed archive containing the data and commments.    

    Examples
    --------    
    written = sunpy.io.ana.write(filename, data, compress=1, comments=False)
    """
    if comments:
        return _pyana.fzwrite(filename, data, compress, comments, debug)
    else:
        return _pyana.fzwrite(filename, data, compress, '', debug)