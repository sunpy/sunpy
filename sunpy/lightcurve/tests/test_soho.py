"""
SOHO/ERNE LightCurve Tests
"""

# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

from __future__ import absolute_import

import pytest
import os

import sunpy.lightcurve
from sunpy.time import TimeRange
import sunpy.data.test

filepath = sunpy.data.test.rootdir

class TestERNELightCurve(object):

    #Test header parsing from file
    @pytest.mark.parametrize("data_file",[('erne/cr1907a.txt'),('erne/cr1907p.txt')])
    @pytest.mark.online
    def test_header(self, data_file):
        lc = sunpy.lightcurve.ERNELightCurve
        assert lc._parse_txt(os.path.join(filepath , data_file))[0] == ['TimeRange', 
                            'energy channel 1.8-3.3 MeV', 
                            'energy channel 3.3-6.4 MeV', 
                            'energy channel 6.4-12.7 MeV', 
                            'energy channel 13.5-25.8 MeV', 
                            'energy channel 25.8-50.7 MeV']


    #Test for non empty data parsing from file
    @pytest.mark.parametrize("data_file",[('erne/cr1907a.txt'),('erne/cr1907p.txt')])
    @pytest.mark.online
    def test_data(self, data_file):
        lc = sunpy.lightcurve.ERNELightCurve
        assert (lc._parse_txt(os.path.join(filepath , data_file))[1]).empty == False


    #Test parsed header and data columns list for equal lengths
    @pytest.mark.parametrize("data_file",[('erne/cr1907a.txt'),('erne/cr1907p.txt')])
    @pytest.mark.online
    def test_header(self, data_file):
        lc = sunpy.lightcurve.ERNELightCurve._parse_txt(os.path.join(filepath , data_file))
        assert len(lc[0]) == len(lc[1].columns) #length of header list equals number of columns


