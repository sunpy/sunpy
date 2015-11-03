
"""
RHESSI LightCurve Tests
"""
from __future__ import absolute_import

import pytest
import sunpy.lightcurve
from sunpy.time import TimeRange
from numpy import all


class TestRHESSISummaryLightCurve(object):

    @pytest.fixture
    def timerange_a(self):
        return TimeRange('2008/06/01', '2008/06/02')

    @pytest.fixture
    def timerange_b(self):
        return TimeRange('2004/06/03', '2004/06/04')

    @pytest.mark.online
    def test_hsi_range(self, timerange_a):
        """Test creation with two times"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a.start, timerange_a.end)
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)

    @pytest.mark.online
    def test_hsi_timerange(self, timerange_a):
        """Test creation with a TimeRange"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a)
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)

    @pytest.mark.online
    def test_hsi_default(self):
        """Test creation with no input"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create()
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)

    @pytest.mark.online
    def test_data(self, timerange_a, timerange_b):
        """Test presence of data"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_b)
        lc2 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a)
        assert not lc1.data.empty
        assert not lc2.data.empty

    @pytest.mark.online
    def test_header(self, timerange_b):
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
    def test_filename(self, timerange_a, timerange_b):
        """Compare data from two different time ranges to make
        sure they are not the same"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a)
        lc2 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_b)
        assert not all(lc1.data[lc1.data.columns[0]] == lc2.data[lc2.data.columns[0]])

    @pytest.mark.online
    def test_get_url(self, timerange_a, timerange_b):
        """Test the getting of urls"""
        g = sunpy.lightcurve.RHESSISummaryLightCurve
        assert g._get_url_for_date_range(timerange_a) == 'http://hesperia.gsfc.nasa.gov/hessidata/metadata/catalog/hsi_obssumm_20080601_068.fits'
        assert g._get_url_for_date_range(timerange_b) == 'http://hesperia.gsfc.nasa.gov/hessidata/metadata/catalog/hsi_obssumm_20040603_110.fits'
