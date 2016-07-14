"""
This module implements tests for SRS Reader.
"""
import sunpy.io.srs, pytest
from sunpy.io.srs import read

import sunpy.data.test
import os


testpath = sunpy.data.test.rootdir

srs_20150906 = os.path.join(testpath, '20150906SRS.txt')
srs_20150306 = os.path.join(testpath, '20150306SRS.txt')
srs_20150101 = os.path.join(testpath, '20150101SRS.txt')
@pytest.mark.paramterize("path_, number_of_rows",
                    [(srs_20150906, 5),
                     (srs_20150306, 4),
                     (srs_20150101, 9)])
def test_number_of_rows(path_, number_of_rows):
    table = sunpy.io.srs.read(path_)
    assert len(table) == number_of_rows
    
