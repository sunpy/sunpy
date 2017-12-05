"""
NOAA LightCurve Tests
"""
from __future__ import absolute_import

import pytest
import sunpy.lightcurve
from sunpy.time import TimeRange

timerange_a = TimeRange('2004/01/01', '2007/01/01')

class TestNOAAIndicesLightCurve(object):

    @pytest.mark.remote_data
    def test_create(self):
        lc = sunpy.lightcurve.NOAAIndicesLightCurve.create()
        assert isinstance(lc, sunpy.lightcurve.NOAAIndicesLightCurve)

    @pytest.mark.remote_data
    def test_isempty(self):
        lc = sunpy.lightcurve.NOAAIndicesLightCurve.create()
        assert lc.data.empty == False

    @pytest.mark.remote_data
    def test_url(self):
        """Test creation with url"""
        url = 'ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt'
        lc1 = sunpy.lightcurve.NOAAIndicesLightCurve.create(url)
        assert isinstance(lc1, sunpy.lightcurve.NOAAIndicesLightCurve)

    @pytest.mark.remote_data
    def test_goes_timerange(self):
        """Test creation with a TimeRange"""
        lc1 = sunpy.lightcurve.NOAAIndicesLightCurve.create(timerange_a)
        assert isinstance(lc1, sunpy.lightcurve.NOAAIndicesLightCurve)

    def test_get_url(self):
        """Test the getting of url"""
        g = sunpy.lightcurve.NOAAIndicesLightCurve
        assert g._get_url_for_date_range(timerange_a) == 'ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt'

    @pytest.mark.remote_data
    def test_header(self):
        """Test presence of GOES satellite number in header"""
        lc1 = sunpy.lightcurve.NOAAIndicesLightCurve.create()
        assert 'comments' in lc1.header.keys()


class TestNOAAPredictIndicesLightCurve(object):

    @pytest.mark.remote_data
    def test_create(self):
        """Test creation with no input"""
        lc = sunpy.lightcurve.NOAAPredictIndicesLightCurve.create()
        assert isinstance(lc, sunpy.lightcurve.NOAAPredictIndicesLightCurve)

    @pytest.mark.remote_data
    def test_isempty(self):
        """Test presence of data"""
        lc = sunpy.lightcurve.NOAAPredictIndicesLightCurve.create()
        assert lc.data.empty == False

    @pytest.mark.remote_data
    def test_goes_timerange(self):
        """Test creation with a TimeRange"""
        lc1 = sunpy.lightcurve.NOAAPredictIndicesLightCurve.create(timerange_a)
        assert isinstance(lc1, sunpy.lightcurve.NOAAPredictIndicesLightCurve)

    @pytest.mark.remote_data
    def test_url(self):
        """Test creation with url"""
        url = 'http://services.swpc.noaa.gov/text/predicted-sunspot-radio-flux.txt'
        lc1 = sunpy.lightcurve.NOAAPredictIndicesLightCurve.create(url)
        assert isinstance(lc1, sunpy.lightcurve.NOAAPredictIndicesLightCurve)

    def test_get_url(self):
        """Test the getting of url"""
        g = sunpy.lightcurve.NOAAPredictIndicesLightCurve
        assert g._get_url_for_date_range(timerange_a) == 'http://services.swpc.noaa.gov/text/predicted-sunspot-radio-flux.txt'

    @pytest.mark.remote_data
    def test_header(self):
        """Test presence of GOES satellite number in header"""
        lc1 = sunpy.lightcurve.NOAAPredictIndicesLightCurve.create()
        assert 'comments' in lc1.header.keys()
