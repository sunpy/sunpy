# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:08:21 2016

@author: alex_
"""

import sunpy.data.sample
import sunpy.timeseries
    
#==============================================================================
# Map Factory Tests
#==============================================================================
class TestTimeSeries(object):
    def test_concatenate(self):
        #Test making a TimeSeries that is the concatenation of multiple files.
        ts = sunpy.timeseries.TimeSeries(a_list_of_many, concatenate=True)
        assert isinstance(ts, sunpy)



#==============================================================================
# Sources Tests
#==============================================================================
    def test_eve(self):
        #Test an EVE TimeSeries
        ts_eve = sunpy.timeseries.TimeSeries(sunpy.data.sample.EVE_LIGHTCURVE, source='EVE')
        assert isinstance(ts_eve, sunpy.timeseries.sources.eve.EVELightCurve)

    def test_fermi_gbm(self):
        #Test a GBMSummary TimeSeries
        ts_gbm = sunpy.timeseries.TimeSeries(sunpy.data.sample.GBM_LIGHTCURVE, source='GBMSummary')
        assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryLightCurve)

    def test_norrh(self):
        #Test a NoRH TimeSeries
        ts_norrh = sunpy.timeseries.TimeSeries(sunpy.data.sample.NORH_LIGHTCURVE, source='NoRH')
        assert isinstance(ts_norrh, sunpy.timeseries.sources.norh.NoRHLightCurve)

    def test_goes(self):
        #Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_LIGHTCURVE, source='GOES')
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.GOESLightCurve)

    def test_lyra(self):
        #Test a LYRA TimeSeries
        ts_lyra = sunpy.timeseries.TimeSeries(sunpy.data.sample.LYRA_LEVEL3_LIGHTCURVE, source='LYRA')
        assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRALightCurve)

    def test_rhessi(self):
        #Test a RHESSI TimeSeries
        ts_rhessi = sunpy.timeseries.TimeSeries(sunpy.data.sample.RHESSI_LIGHTCURVE, source='RHESSI')
        assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryLightCurve)

    def test_noaa_ind(self):
        #Test a NOAAPredictIndices TimeSeries
        ts_noaa_ind = sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAINDICES_LIGHTCURVE, source='NOAAPredictIndices')
        assert isinstance(ts_noaa_ind, sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries)
            
    def test_noaa_pre(self):
        #Test a NOAAIndices TimeSeries
        ts_noaa_pre = sunpy.timeseries.TimeSeries(sunpy.data.sample.NOAAPREDICT_LIGHTCURVE, source='NOAAIndices')
        assert isinstance(ts_noaa_pre, sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries)
    
    
    