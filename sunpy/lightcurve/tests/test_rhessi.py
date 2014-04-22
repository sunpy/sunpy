
"""
GOES LightCurve Tests
"""
from __future__ import absolute_import

import pytest
import sunpy.lightcurve
from sunpy.time import TimeRange

timerange_a = TimeRange('2008/06/01', '2008/06/02')
timerange_b = TimeRange('2004/06/03', '2004/06/04')

class TestRHESSISummaryLightCurve():
    
    @pytest.mark.online
    def test_hsi_range(self):
        """Test creation with two times"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a.start(), timerange_a.end())
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)

    @pytest.mark.online
    def test_hsi_timerange(self):
        """Test creation with a TimeRange"""
        lc1 = sunpy.lightcurve.GOESLightCurve.create(timerange_a)
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)
    
    @pytest.mark.online
    def test_hsi_default(self):
        """Test creation with no input"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create()
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)
    
    @pytest.mark.online
    def test_data(self):
        """Test presence of data"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_b)
        lc2 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a)
        assert lc1.data.empty == False
        assert lc2.data.empty == False

    @pytest.mark.online
    def test_header(self):
        """Test presence of TELESCOP in header"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_b)
        assert lc1.header['TELESCOP'] == 'HESSI'
    
    @pytest.mark.online
    def test_hsi_url(self):
        """Test creation with url"""
        url = 'http://hesperia.gsfc.nasa.gov/hessidata/metadata/catalog/hsi_obssumm_20030302_146.fits'
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(url)
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)
    
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
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a)
        lc2 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_b)
        #If the dataframes are non-idential it raises an error, if they are
        #identical it returns True
        with pytest.raises((Exception)):
            self.compare(lc1, lc2)
        
    def test_get_url(self):
        """Test the getting of urls"""
        g = sunpy.lightcurve.RHESSISummaryLightCurve
        assert g._get_url_for_date_range(timerange_a) == 'http://hesperia.gsfc.nasa.gov/hessidata/metadata/catalog/hsi_obssumm_20080601_068.fits'
        assert g._get_url_for_date_range(timerange_b) == 'http://hesperia.gsfc.nasa.gov/hessidata/metadata/catalog/hsi_obssumm_20040603_110.fits'