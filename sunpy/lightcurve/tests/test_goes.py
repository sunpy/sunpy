
"""
GOES LightCurve Tests
"""
from __future__ import absolute_import

import pytest
import sunpy.lightcurve
from sunpy.time import TimeRange

timerange_a = TimeRange('2008/06/01', '2008/06/02')
timerange_b = TimeRange('1995/06/03', '1995/06/04')

class TestGOESLightCurve():
    
    @pytest.mark.online
    def test_goes_range(self):
        """Test creation with two times"""
        lc1 = sunpy.lightcurve.GOESLightCurve.create(timerange_a.start(), timerange_a.end())
        assert isinstance(lc1, sunpy.lightcurve.GOESLightCurve)

    @pytest.mark.online
    def test_goes_timerange(self):
        """Test creation with a TimeRange"""
        lc1 = sunpy.lightcurve.GOESLightCurve.create(timerange_a)
        assert isinstance(lc1, sunpy.lightcurve.GOESLightCurve)
    
    @pytest.mark.online
    def test_goes_default(self):
        """Test creation with no input"""
        lc1 = sunpy.lightcurve.GOESLightCurve.create()
        assert isinstance(lc1, sunpy.lightcurve.GOESLightCurve)
    
    @pytest.mark.online
    def test_data(self):
        """Test presence of data"""
        lc1 = sunpy.lightcurve.GOESLightCurve.create(timerange_b)
        lc2 = sunpy.lightcurve.GOESLightCurve.create(timerange_a)
        assert lc1.data.empty == False
        assert lc2.data.empty == False

    @pytest.mark.online
    def test_header(self):
        """Test presence of GOES satellite number in header"""
        lc1 = sunpy.lightcurve.GOESLightCurve.create(timerange_b)
        lc2 = sunpy.lightcurve.GOESLightCurve.create(timerange_a)
        assert lc1.header['TELESCOP'] == 'GOES 7'
        assert lc2.header['TELESCOP'] == 'GOES 10'
    
    def test_goes_url(self):
        """Test creation with url"""
        url = 'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits'
        lc1 = sunpy.lightcurve.GOESLightCurve.create(url)
        assert isinstance(lc1, sunpy.lightcurve.GOESLightCurve)
    
    @pytest.mark.online
    def compare(self,lc1,lc2):
        try:
            (lc1.data == lc2.data)
        except:
            raise Exception

    @pytest.mark.online
    def test_filename(self):
        """Compare data from two different time ranges to make 
        sure they are not the same"""
        lc1 = sunpy.lightcurve.GOESLightCurve.create(timerange_a)
        lc2 = sunpy.lightcurve.GOESLightCurve.create(timerange_b)
        #If the dataframes are non-idential it raises an error, if they are
        #identical it returns True
        with pytest.raises((Exception)):
            self.compare(lc1, lc2)
    
    def test_goes_sat_numbers(self):
        """Test the ability to return GOES satellite availability"""
        g = sunpy.lightcurve.GOESLightCurve
        assert g._get_goes_sat_num(timerange_a.start(), timerange_a.end()) == [10]
        assert g._get_goes_sat_num(timerange_b.start(), timerange_b.end()) == [7]
    
    def test_get_url(self):
        """Test the getting of urls"""
        g = sunpy.lightcurve.GOESLightCurve
        # time ranges create urls with either 4 digit or 2 digit years
        assert g._get_url_for_date_range(timerange_b) == 'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits'
        assert g._get_url_for_date_range(timerange_a) == 'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080601.fits'