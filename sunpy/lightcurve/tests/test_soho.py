"""
SOHO/ERNE LightCurve Tests
"""

# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

from __future__ import absolute_import

import pytest

import sunpy.lightcurve
from sunpy.time import TimeRange


class TestERNELightCurve(object):

    @pytest.mark.online
    def test_header(self):
        """Test header parsing from file"""
        lc = sunpy.lightcurve.ERNELightCurve
        assert lc._parse_txt('<filepath_to_downloaded_file>')[0] == [ 'TimeRange', 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 1.8-3.3  [MeV] '  , 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 3.3-6.4  [MeV] '  , 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 6.4-12.7  [MeV] ' , 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 13.5-25.8  [MeV]' , 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 25.8-50.7  [MeV]' ]

    @pytest.mark.online
    def test_data(self):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.ERNELightCurve
        assert (lc._parse_txt('<filepath_to_downloaded_file>')[1]).empty == False


