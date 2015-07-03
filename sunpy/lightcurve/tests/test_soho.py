"""
SOHO/ERNE LightCurve Tests
"""
from __future__ import absolute_import

import pytest

import sunpy.lightcurve
from sunpy.time import TimeRange

timerange_a = TimeRange('2008/01/01', '2010/01/01')
specie_a = 'proton'
specie_b = 'alpha'
test_file_a = 'cr1907p.txt'
test_file_b = 'cr1907a.txt'


class TestERNELightCurve(object):

    @pytest.mark.online
    def test_input_file(self):
        """Test creation with a TimeRange"""
        lc1 = sunpy.lightcurve.ERNELightCurve.create(timerange_a, specie_a)
        assert isinstance(lc1, sunpy.lightcurve.ERNELightCurve)

    @pytest.mark.online
    def test_isempty(self):
        lc = sunpy.lightcurve.ERNELightCurve.create(timerange_a, specie_b)
        assert lc.data.empty == False

    @pytest.mark.online
    def test_header(self):
        """Test header parsing from sample file"""
        info = sunpy.lightcurve.ERNELightCurve._parse_txt('<filepath_to_test_file>')
        assert info[1] == [ 'TimeRange', 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 1.8-3.3  [MeV] '  , 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 3.3-6.4  [MeV] '  , 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 6.4-12.7  [MeV] ' , 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 13.5-25.8  [MeV]' , 
                            'Intensities [1/(cm^2*sr*s*MeV)] in energy channel 25.8-50.7  [MeV]' ]


