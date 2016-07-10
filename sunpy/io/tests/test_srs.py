"""
This module implements tests for SRS Reader.
"""
import sunpy.io.srs
from sunpy.io.srs import read

import sunpy.data.test
import os


testpath = sunpy.data.test.rootdir

srs_20150906 = os.path.join(testpath, '20150906SRS.txt')
srs_20150306 = os.path.join(testpath, '20150306SRS.txt')
srs_20150101 = os.path.join(testpath, '20150101SRS.txt')

def test_20150906():
    table = sunpy.io.srs.read(srs_20150906)
    assert len(table) == 5

def test_20150306():
    table = sunpy.io.srs.read(srs_20150306)
    assert len(table) == 4

def test_20150101():
    table = sunpy.io.srs.read(srs_20150101)
    assert len(table) == 9
    
