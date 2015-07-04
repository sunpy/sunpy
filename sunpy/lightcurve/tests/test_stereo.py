"""
STEREO LightCurve Tests
"""

# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

from __future__ import absolute_import

import pytest

import sunpy.lightcurve
from sunpy.time import TimeRange


class TestLETLightCurve(object):

    @pytest.mark.online
    def test_header(self):
        """Test header parsing from file"""
        lc = sunpy.lightcurve.LETLightCurve
        assert (lc._parse_txt('<filepath_to_downloaded_file>')[0])[1][:14] in ['Flux for Bin 0', 'Column 6: LET '] 

    @pytest.mark.online
    def test_data(self):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.LETLightCurve
        assert (lc._parse_txt('<filepath_to_downloaded_file>')[1]).empty == False


class TestSITLightCurve(object):

    @pytest.mark.online
    def test_header(self):
        """Test parsed header and data columns list for equal lengths """
        lc = sunpy.lightcurve.SITLightCurve._parse_txt('<filepath_to_downloaded_file>')
        assert len(lc[0]) == len(lc[1].columns) #length of header list equals number of columns

    @pytest.mark.online
    def test_data(self):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.SITLightCurve
        assert (lc._parse_txt('<filepath_to_downloaded_file>')[1]).empty == False


class TestHETLightCurve(object):

    @pytest.mark.online
    def test_header(self):
        """Test parsed header and data columns list for equal lengths """
        lc = sunpy.lightcurve.HETLightCurve._parse_txt('<filepath_to_downloaded_file>')
        assert len(lc[0]) == len(lc[1].columns) #length of header list equals number of columns

    @pytest.mark.online
    def test_data(self):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.HETLightCurve
        assert (lc._parse_txt('<filepath_to_downloaded_file>')[1]).empty == False


class TestPLASTICLightCurve(object):

    @pytest.mark.online
    def test_header(self):
        """Test header parsing from file"""
        lc = sunpy.lightcurve.PLASTICLightCurve
        assert lc._parse_txt('<filepath_to_downloaded_file>')[0][-2:] == [ 'Spcrft. Long. [HCI, degrees]', 'Spcrft. Lat. [HCI, degrees]' ]

    @pytest.mark.online
    def test_data(self):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.PLASTICLightCurve
        assert (lc._parse_txt('<filepath_to_downloaded_file>')[1]).empty == False

    @pytest.mark.online
    def test_header(self):
        """Test parsed header and data columns list for equal lengths """
        lc = sunpy.lightcurve.PLASTICLightCurve._parse_txt('<filepath_to_downloaded_file>')
        assert len(lc[0]) == len(lc[1].columns) #length of header list equals number of columns

class TestSEPTLightCurve(object):

    @pytest.mark.online
    def test_data(self):
        """Test for non empty data parsing from file"""l
        lc = sunpy.lightcurve.SEPTLightCurve
        assert (lc._parse_txt('<filepath_to_downloaded_file>')[1]).empty == False

    @pytest.mark.online
    def test_header(self):
        """Test parsed header and data columns list for equal lengths """
        lc = sunpy.lightcurve.SEPTLightCurve._parse_txt('<filepath_to_downloaded_file>')
        assert len(lc[0]) == len(lc[1].columns) #length of header list equals number of columns


