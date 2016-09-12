# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:08:21 2016

@author: alex_
"""

import os
import glob
import pytest
import sunpy.data.sample
import sunpy.data.test
import sunpy.timeseries
import datetime
import numpy as np
from pandas import DataFrame
from collections import OrderedDict
from sunpy.util.metadata import MetaDict
import astropy.units as u
from astropy.table import Table
from astropy.time import Time


#==============================================================================
# TimeSeries Factory Tests
#==============================================================================

filepath = sunpy.data.test.rootdir
a_list_of_many = glob.glob(os.path.join(filepath, "goes", "*"))

eve_filepath = os.path.join(filepath, '')
fermi_gbm_filepath = os.path.join(filepath, 'gbm.fits')
norrh_filepath = os.path.join(filepath, '')
goes_filepath = os.path.join(filepath, 'goes.fits')
lyra_filepath = os.path.join(filepath, '')
rhessi_filepath = os.path.join(filepath, '')
noaa_ind_filepath = os.path.join(filepath, '')
noaa_pre_filepath = os.path.join(filepath, '')

goes_filepath = sunpy.data.sample.GOES_LIGHTCURVE
eve_filepath = sunpy.data.sample.EVE_LIGHTCURVE
norrh_filepath = sunpy.data.sample.NORH_LIGHTCURVE
lyra_filepath = sunpy.data.sample.LYRA_LEVEL3_LIGHTCURVE
rhessi_filepath = sunpy.data.sample.RHESSI_LIGHTCURVE
noaa_ind_filepath = sunpy.data.sample.NOAAINDICES_LIGHTCURVE
noaa_pre_filepath = sunpy.data.sample.NOAAPREDICT_LIGHTCURVE

a_list_of_many_goes = ['C:\\Users\\alex_\\sunpy\\data\\go1420101102.fits',
 'C:\\Users\\alex_\\sunpy\\data\\go1420101103.fits',
 'C:\\Users\\alex_\\sunpy\\data\\go1420101104.fits',
 'C:\\Users\\alex_\\sunpy\\data\\go1520101105.fits',
 'C:\\Users\\alex_\\sunpy\\data\\go1520101106.fits',
 'C:\\Users\\alex_\\sunpy\\data\\go1520101107.fits']
a_list_of_many = glob.glob(os.path.join(filepath, "eve", "*"))


class TestTimeSeries(object):
    def test_factory_concatenate(self):
        # Test making a TimeSeries that is the concatenation of multiple files
        ts_from_list = sunpy.timeseries.TimeSeries(a_list_of_many, source='EVE', concatenate=True)
        assert isinstance(ts_from_list, sunpy.timeseries.sources.eve.EVELightCurve)
        ts_from_folder = sunpy.timeseries.TimeSeries(os.path.join(filepath, "eve"), source='EVE', concatenate=True)
        assert isinstance(ts_from_folder, sunpy.timeseries.sources.eve.EVELightCurve)
        assert ts_from_list == ts_from_folder

    def test_factory_generate_list_of_ts(self):
        # Test making a list TimeSeries from multiple files
        ts_list = sunpy.timeseries.TimeSeries(a_list_of_many, source='EVE')
        assert isinstance(ts_list, list)
        for ts in ts_list:
          assert isinstance(ts, sunpy.timeseries.sources.eve.EVELightCurve)

    def test_factory_generate_from_glob(self):
        # Test making a TimeSeries from a glob
        ts_from_glob = sunpy.timeseries.TimeSeries(os.path.join(filepath, "eve", "*"), source='EVE', concatenate=True)
        assert isinstance(ts_from_glob, sunpy.timeseries.sources.eve.EVELightCurve)

#==============================================================================
# Sources Tests
#==============================================================================

    def test_eve(self):
        #Test an EVE TimeSeries
        ts_eve = sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')
        assert isinstance(ts_eve, sunpy.timeseries.sources.eve.EVELightCurve)

    def test_fermi_gbm(self):
        #Test a GBMSummary TimeSeries
        ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath)
        assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryLightCurve)

    def test_norrh(self):
        #Test a NoRH TimeSeries
        ts_norrh = sunpy.timeseries.TimeSeries(norrh_filepath, source='NoRH')
        assert isinstance(ts_norrh, sunpy.timeseries.sources.norh.NoRHLightCurve)

    def test_goes(self):
        #Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(goes_filepath, source='GOES')
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.GOESLightCurve)

    def test_lyra(self):
        #Test a LYRA TimeSeries
        ts_lyra = sunpy.timeseries.TimeSeries(lyra_filepath, source='LYRA')
        assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRALightCurve)

    def test_rhessi(self):
        #Test a RHESSI TimeSeries
        ts_rhessi = sunpy.timeseries.TimeSeries(rhessi_filepath, source='RHESSI')
        assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryLightCurve)

    def test_noaa_ind(self):
        #Test a NOAAPredictIndices TimeSeries
        ts_noaa_ind = sunpy.timeseries.TimeSeries(noaa_ind_filepath, source='NOAAIndices')
        assert isinstance(ts_noaa_ind, sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries)

    def test_noaa_pre(self):
        #Test a NOAAIndices TimeSeries
        ts_noaa_pre = sunpy.timeseries.TimeSeries(noaa_pre_filepath, source='NOAAPredictIndices')
        assert isinstance(ts_noaa_pre, sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries)

    def test_generic_construction(self):
        # Generate the data and the corrisponding dates
        base = datetime.datetime.today()
        times = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
        intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity, index=times, columns=['intensity'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        meta = MetaDict({'key':'value'})

        # Create normal TS from dataframe and check
        ts_generic = sunpy.timeseries.TimeSeries(data, meta, units)
        assert isinstance(ts_generic, sunpy.timeseries.timeseriesbase.GenericTimeSeries)

        # Create TS omitting some input arguments
        ts_1 = sunpy.timeseries.TimeSeries(data, meta)
        assert isinstance(ts_1, sunpy.timeseries.timeseriesbase.GenericTimeSeries)
        ts_2 = sunpy.timeseries.TimeSeries(data, units)
        assert isinstance(ts_2, sunpy.timeseries.timeseriesbase.GenericTimeSeries)

        # Create TS using a tuple of values
        ts_3 = sunpy.timeseries.TimeSeries((data, meta, units))
        assert isinstance(ts_3, sunpy.timeseries.timeseriesbase.GenericTimeSeries)

    def test_table_to_ts(self):
        # Generate the data and the corresponding dates
        base = datetime.datetime.today()
        times = Time([base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)])
        intensity = u.Quantity(np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60)))), u.W/u.m**2)

        # Create the units and meta objects
        units = OrderedDict([('intensity', u.W/u.m**2)])
        meta = MetaDict({'key':'value'})
        tbl_meta = MetaDict({'t_key':'t_value'})

        # Create a suitable mixin qtable
        table = Table([times, intensity], names=['time', 'intensity'], meta=tbl_meta)
        table.add_index('time')

        # Create TS from table and check
        ts_table = sunpy.timeseries.TimeSeries(table, meta, units)
        assert isinstance(ts_table, sunpy.timeseries.timeseriesbase.GenericTimeSeries)
        ts_table2 = sunpy.timeseries.TimeSeries(table, units, meta)
        assert (ts_table2 == ts_table)

        # Create TS using a tuple of values
        ts_table3 = sunpy.timeseries.TimeSeries((table, meta, units))
        assert isinstance(ts_table3, sunpy.timeseries.timeseriesbase.GenericTimeSeries)

        # ToDo: Try an incompatible table
        dual_index_table = Table([times, intensity], names=['time', 'intensity'], meta=tbl_meta)
        dual_index_table.add_index(('time', 'intensity'))
        with pytest.raises(ValueError):
            sunpy.timeseries.TimeSeries((dual_index_table, meta, units))


#==============================================================================
# Test some Errors
#==============================================================================


    def test_validate_units(self):
        valid_units = OrderedDict([('Watt Per Meter Squared', u.Unit("W / m2")), ('Meter Cubed', u.Unit("m3"))])
        assert sunpy.timeseries.TimeSeries._validate_units(valid_units)
        # Test for not having only units for values
        invalid_units_1 = OrderedDict([('Watt Per Meter Squared', 'string'), ('Meter Cubed', u.Unit("m3"))])
        assert not sunpy.timeseries.TimeSeries._validate_units(invalid_units_1)
        # Test for being a MetaDict object
        invalid_units_2 = MetaDict(OrderedDict([('Watt Per Meter Squared', u.Unit("W / m2")), ('Meter Cubed', u.Unit("m3"))]))
        assert not sunpy.timeseries.TimeSeries._validate_units(invalid_units_2)

    def test_validate_meta(self):
        valid_meta_1 = MetaDict({'key':'value'})
        assert sunpy.timeseries.TimeSeries._validate_meta(valid_meta_1)
        valid_meta_2 = OrderedDict({'key':'value'})
        assert sunpy.timeseries.TimeSeries._validate_meta(valid_meta_2)
