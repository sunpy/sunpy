
"""
RHESSI LightCurve Tests
"""
from __future__ import absolute_import

import pytest
import sunpy.lightcurve
from sunpy.time import TimeRange
import numpy as np


class TestRHESSISummaryLightCurve(object):

    @pytest.fixture
    def timerange_a(self):
        return TimeRange('2008/06/01', '2008/06/02')

    @pytest.fixture
    def timerange_b(self):
        return TimeRange('2004/06/03', '2004/06/04')

    @pytest.mark.remote_data
    def test_hsi_range(self, timerange_a):
        """Test creation with two times"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a.start, timerange_a.end)
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)

    @pytest.mark.remote_data
    def test_hsi_timerange(self, timerange_a):
        """Test creation with a TimeRange"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a)
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)

    @pytest.mark.remote_data
    def test_hsi_default(self):
        """Test creation with no input"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create()
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)

    @pytest.mark.remote_data
    def test_data(self, timerange_a, timerange_b):
        """Test presence of data"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_b)
        lc2 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a)
        assert not lc1.data.empty
        assert not lc2.data.empty

    @pytest.mark.remote_data
    def test_header(self, timerange_b):
        """Test presence of TELESCOP in header"""
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_b)
        assert lc1.header['TELESCOP'] == 'HESSI'

    @pytest.mark.remote_data
    def test_hsi_url(self):
        """Test creation with url"""
        url = 'http://soleil.i4ds.ch/hessidata/metadata/catalog/hsi_obssumm_20030302_146.fits'
        lc1 = sunpy.lightcurve.RHESSISummaryLightCurve.create(url)
        assert isinstance(lc1, sunpy.lightcurve.RHESSISummaryLightCurve)

    @pytest.mark.remote_data
    def test_filename(self, timerange_a, timerange_b):
        """Compare data from two different time ranges to make
        sure they are not the same"""
        lc1 = np.array(sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_a).data)
        lc2 = np.array(sunpy.lightcurve.RHESSISummaryLightCurve.create(timerange_b).data)
        with pytest.raises(AssertionError):
            np.testing.assert_allclose(lc1, lc2)

    @pytest.mark.remote_data
    def test_get_url(self, timerange_a, timerange_b):
        """Test the getting of urls"""
        g = sunpy.lightcurve.RHESSISummaryLightCurve
        assert 'hessidata/metadata/catalog/hsi_obssumm_20080601_068.fits' in g._get_url_for_date_range(timerange_a)
        assert 'hessidata/metadata/catalog/hsi_obssumm_20040603_110.fits' in g._get_url_for_date_range(timerange_b)
