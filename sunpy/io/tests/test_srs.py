"""
This module implements tests for SRS Reader.
"""
import os

import pytest

from sunpy.io.special import srs
import sunpy.data.test

testpath = sunpy.data.test.rootdir

filenames = [{'file': '20150906SRS.txt', 'rows': 5},
             {'file': '20150306SRS.txt', 'rows': 4},
             {'file': '20150101SRS.txt', 'rows': 9}]

@pytest.mark.parametrize("path, number_of_rows",
                         [(os.path.join(testpath, elem['file']), elem['rows'])
                          for elem in filenames])
def test_number_of_rows(path, number_of_rows):
    table = srs.read(path)
    assert len(table) == number_of_rows

