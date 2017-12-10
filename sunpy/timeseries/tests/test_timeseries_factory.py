# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 12:08:21 2016

@author: alex_
"""

import os
import glob
import pytest
import datetime
import numpy as np
from pandas import DataFrame
from collections import OrderedDict

import sunpy.data.test
import sunpy.timeseries
from sunpy.util.metadata import MetaDict
import sunpy.io
from sunpy.util.datatype_factory_base import NoMatchError

import astropy.units as u
from astropy.table import Table
from astropy.time import Time
from astropy.io import fits

# ==============================================================================
# TimeSeries Factory Tests
# ==============================================================================

filepath = sunpy.data.test.rootdir
eve_filepath = os.path.join(filepath, 'EVE_L0CS_DIODES_1m_truncated.txt')
fermi_gbm_filepath = os.path.join(filepath, 'gbm.fits')
norh_filepath = os.path.join(filepath, 'tca110810_truncated')
lyra_filepath = os.path.join(filepath, 'lyra_20150101-000000_lev3_std_truncated.fits.gz')
rhessi_filepath = os.path.join(filepath, 'hsi_obssumm_20120601_018_truncated.fits.gz')
noaa_ind_filepath = os.path.join(filepath, 'RecentIndices_truncated.txt')
noaa_pre_filepath = os.path.join(filepath, 'predicted-sunspot-radio-flux_truncated.txt')
goes_filepath_com = os.path.join(filepath, 'go1520120601.fits.gz')
goes_filepath = os.path.join(filepath, 'go1520110607.fits')
a_list_of_many = glob.glob(os.path.join(filepath, "eve", "*"))

# ==============================================================================
# Multi file Tests
# ==============================================================================


class TestTimeSeries(object):
    def test_factory_concatenate_same_source(self):
        # Test making a TimeSeries that is the concatenation of multiple files
        ts_from_list = sunpy.timeseries.TimeSeries(a_list_of_many, source='EVE', concatenate=True)
        assert isinstance(ts_from_list, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
        ts_from_folder = sunpy.timeseries.TimeSeries(os.path.join(filepath, "eve"), source='EVE', concatenate=True)
        assert isinstance(ts_from_folder, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
        # text the two methods get identical dataframes
        assert ts_from_list == ts_from_folder
        # test the frames have correct headings/keys (correct concatenation axis)
        ts_from_list.columns == sunpy.timeseries.TimeSeries(a_list_of_many[0], source='EVE', concatenate=True).columns

    def test_factory_concatenate_different_source(self):
        # Test making a TimeSeries that is the concatenation of multiple files
        ts_from_list = sunpy.timeseries.TimeSeries(a_list_of_many, source='EVE', concatenate=True)
        assert isinstance(ts_from_list, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
        ts_from_folder = sunpy.timeseries.TimeSeries(os.path.join(filepath, "eve"), source='EVE', concatenate=True)
        assert isinstance(ts_from_folder, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
        # text the two methods get identical dataframes
        assert ts_from_list == ts_from_folder
        # test the frames have correct headings/keys (correct concatenation axis)
        ts_from_list.columns == sunpy.timeseries.TimeSeries(a_list_of_many[0], source='EVE', concatenate=True).columns

    def test_factory_generate_list_of_ts(self):
        # Test making a list TimeSeries from multiple files
        ts_list = sunpy.timeseries.TimeSeries(a_list_of_many, source='EVE')
        assert isinstance(ts_list, list)
        for ts in ts_list:
          assert isinstance(ts, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)

    def test_factory_generate_from_glob(self):
        # Test making a TimeSeries from a glob
        ts_from_glob = sunpy.timeseries.TimeSeries(os.path.join(filepath, "eve", "*"), source='EVE', concatenate=True)
        assert isinstance(ts_from_glob, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)

#==============================================================================
# Individual Implicit Source Tests
#==============================================================================

    def test_implicit_fermi_gbm(self):
        # Test a GBMSummary TimeSeries
        ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath)
        assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)

    def test_implicit_norh(self):
        # Test a NoRH TimeSeries
        ts_norh = sunpy.timeseries.TimeSeries(norh_filepath)
        assert isinstance(ts_norh, sunpy.timeseries.sources.norh.NoRHTimeSeries)

    def test_implicit_goes(self):
        # Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(goes_filepath)
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)

    def test_implicit_goes_com(self):
        # Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(goes_filepath_com)
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)
        
    def test_implicit_lyra(self):
        # Test a LYRA TimeSeries
        ts_lyra = sunpy.timeseries.TimeSeries(lyra_filepath)
        assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRATimeSeries)

    def test_implicit_rhessi(self):
        # Test a RHESSI TimeSeries
        ts_rhessi = sunpy.timeseries.TimeSeries(rhessi_filepath)
        assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryTimeSeries)

#==============================================================================
# Individual Explicit Sources Tests
#==============================================================================

    def test_eve(self):
        #Test an EVE TimeSeries
        ts_eve = sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')
        assert isinstance(ts_eve, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)

    def test_fermi_gbm(self):
        #Test a GBMSummary TimeSeries
        ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath, source='GBMSummary')
        assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)

    def test_norh(self):
        #Test a NoRH TimeSeries
        ts_norh = sunpy.timeseries.TimeSeries(norh_filepath, source='NoRH')
        assert isinstance(ts_norh, sunpy.timeseries.sources.norh.NoRHTimeSeries)

    def test_goes(self):
        #Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(goes_filepath, source='XRS')
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)

    def test_goes_com(self):
        #Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(goes_filepath_com, source='XRS')
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)
        
    def test_lyra(self):
        #Test a LYRA TimeSeries
        ts_lyra = sunpy.timeseries.TimeSeries(lyra_filepath, source='LYRA')
        assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRATimeSeries)

    def test_rhessi(self):
        #Test a RHESSI TimeSeries
        ts_rhessi = sunpy.timeseries.TimeSeries(rhessi_filepath, source='RHESSI')
        assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryTimeSeries)

    def test_noaa_ind(self):
        #Test a NOAAPredictIndices TimeSeries
        ts_noaa_ind = sunpy.timeseries.TimeSeries(noaa_ind_filepath, source='NOAAIndices')
        assert isinstance(ts_noaa_ind, sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries)

    def test_noaa_pre(self):
        #Test a NOAAIndices TimeSeries
        ts_noaa_pre = sunpy.timeseries.TimeSeries(noaa_pre_filepath, source='NOAAPredictIndices')
        assert isinstance(ts_noaa_pre, sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries)

#==============================================================================
# Manual TimeSeries Tests
#==============================================================================

    def test_meta_from_fits_header(self):
        # Generate the data and the corrisponding dates
        base = datetime.datetime.today()
        times = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
        intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
        data = DataFrame(intensity, index=times, columns=['intensity'])

        # Use a FITS file HDU using sunpy.io
        hdulist = sunpy.io.read_file(goes_filepath)
        meta = hdulist[0].header
        meta_md = MetaDict(OrderedDict(meta))
        ts_hdu_meta = sunpy.timeseries.TimeSeries(data, meta)
        ts_md_meta = sunpy.timeseries.TimeSeries(data, meta_md)
        assert ts_hdu_meta == ts_md_meta

        # Use a FITS file HDU using astropy.io
        hdulist = fits.open(goes_filepath)
        meta = hdulist[0].header
        hdulist.close()
        meta_md = MetaDict(sunpy.io.header.FileHeader(meta))
        ts_hdu_meta = sunpy.timeseries.TimeSeries(data, meta)
        ts_md_meta = sunpy.timeseries.TimeSeries(data, meta_md)
        assert ts_hdu_meta == ts_md_meta

    def test_generic_construction_basic(self):
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
        assert ts_generic.columns == ['intensity']
        assert ts_generic.units == units
        assert ts_generic.meta.metadata[0][2] == meta

        # Create TS using a tuple of values
        ts_tuple = sunpy.timeseries.TimeSeries(((data, meta, units),))
        assert isinstance(ts_tuple, sunpy.timeseries.timeseriesbase.GenericTimeSeries)
        assert ts_generic == ts_tuple


    def test_generic_construction_basic_omitted_details(self):
        # Generate the data and the corrisponding dates
        base = datetime.datetime.today()
        times = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
        intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity, index=times, columns=['intensity'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        meta = MetaDict({'key':'value'})

        # Create TS omitting units input arguments
        ts_1 = sunpy.timeseries.TimeSeries(data, meta)
        assert isinstance(ts_1, sunpy.timeseries.timeseriesbase.GenericTimeSeries)
        assert ts_1.columns == ['intensity']
        assert ts_1.units == OrderedDict([('intensity', u.dimensionless_unscaled)])
        assert ts_1.meta.metadata[0][2] == meta

        ts_2 = sunpy.timeseries.TimeSeries(data, units)
        assert isinstance(ts_2, sunpy.timeseries.timeseriesbase.GenericTimeSeries)
        assert ts_2.columns == ['intensity']
        assert ts_2.units == units
        assert ts_2.meta.metadata[0][2] == MetaDict()

    def test_generic_construction_basic_different_meta_types(self):
        # Generate the data and the corrisponding dates
        base = datetime.datetime.today()
        times = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
        intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity, index=times, columns=['intensity'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        meta_md = MetaDict({'key':'value'})
        meta_di = {'key':'value'}
        meta_od = OrderedDict({'key':'value'})

        # Create TS using different dictionary meta types
        ts_md = sunpy.timeseries.TimeSeries(data, meta_md, units)
        ts_di = sunpy.timeseries.TimeSeries(data, meta_di, units)
        ts_od = sunpy.timeseries.TimeSeries(data, meta_od, units)
        assert ts_md == ts_di == ts_od
        assert ts_md.meta.metadata[0][2] == ts_di.meta.metadata[0][2] == ts_od.meta.metadata[0][2]


    def test_generic_construction_ts_list(self):
        # Generate the data and the corrisponding dates
        base = datetime.datetime.today()
        times = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
        intensity1 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
        intensity2 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity1, index=times, columns=['intensity'])
        data2 = DataFrame(intensity2, index=times, columns=['intensity2'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        units2 = OrderedDict([('intensity', u.W/u.m**2)])
        meta = MetaDict({'key':'value'})
        meta2 = MetaDict({'key2':'value2'})

        # Create TS individually
        ts_1 = sunpy.timeseries.TimeSeries(data, meta, units)
        ts_2 = sunpy.timeseries.TimeSeries(data2, meta2, units2)

        # Create TS list using
        ts_list = sunpy.timeseries.TimeSeries(data, meta, units, data2, meta2, units2)
        assert isinstance(ts_list, list)
        assert len(ts_list) == 2
        assert ts_list[0] == ts_1
        assert ts_list[1] == ts_2

        # Create TS using a tuple
        ts_list2 = sunpy.timeseries.TimeSeries(((data, meta, units),(data2, meta2, units2)))
        assert ts_list == ts_list2

    def test_generic_construction_concatenation(self):
        # Generate the data and the corrisponding dates
        base = datetime.datetime.today()
        times = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
        intensity1 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
        intensity2 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity1, index=times, columns=['intensity'])
        data2 = DataFrame(intensity2, index=times, columns=['intensity2'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        units2 = OrderedDict([('intensity', u.W/u.m**2)])
        meta = MetaDict({'key':'value'})
        meta2 = MetaDict({'key2':'value2'})

        # Create TS individually
        ts_1 = sunpy.timeseries.TimeSeries(data, meta, units)
        ts_2 = sunpy.timeseries.TimeSeries(data2, meta2, units2)
        ts_concat_1 = ts_1.concatenate(ts_2)

        # Concatinate during construction
        ts_concat_2 = sunpy.timeseries.TimeSeries(data, meta, units, data2, meta2, units2, concatenate=True)
        assert isinstance(ts_concat_2, sunpy.timeseries.timeseriesbase.GenericTimeSeries)

        # Create TS using a tuple
        ts_concat_3 = sunpy.timeseries.TimeSeries(((data, meta, units),(data2, meta2, units2)), concatenate=True)
        assert isinstance(ts_concat_3, sunpy.timeseries.timeseriesbase.GenericTimeSeries)
        assert ts_concat_1 == ts_concat_2 == ts_concat_3

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
# Test some other options
#==============================================================================

    def test_passed_ts(self):
        # Test an EVE TimeSeries
        ts_eve = sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')
        ts_from_ts_1 = sunpy.timeseries.TimeSeries(ts_eve, source='EVE')
        ts_from_ts_2 = sunpy.timeseries.TimeSeries(ts_eve)
        assert ts_eve == ts_from_ts_1 == ts_from_ts_2

#==============================================================================
# Test some Errors
#==============================================================================

    def test_invalid_manual_data(self):
        meta = MetaDict({'key':'value'})
        data = []
        with pytest.raises(NoMatchError):
            sunpy.timeseries.TimeSeries(data, meta)

    def test_invalid_filepath(self):
        invalid_filepath = os.path.join(filepath, 'invalid_filepath_here')
        with pytest.raises(NoMatchError):
            sunpy.timeseries.TimeSeries(invalid_filepath)
        # Now with silence_errors kwarg set
        with pytest.raises(NoMatchError):
            sunpy.timeseries.TimeSeries(invalid_filepath, silence_errors=True)

    def test_invalid_file(self):
        invalid_filepath = os.path.join(filepath, 'annotation_ppt.db')
        with pytest.raises(TypeError):
            sunpy.timeseries.TimeSeries(invalid_filepath)
        # Now with silence_errors kwarg set
        with pytest.raises(TypeError):
            sunpy.timeseries.TimeSeries(invalid_filepath, silence_errors=True)

    def test_validate_units(self):
        valid_units = OrderedDict([('Watt Per Meter Squared', u.Unit("W / m2")), ('Meter Cubed', u.Unit("m3"))])
        assert sunpy.timeseries.TimeSeries._validate_units(valid_units)
        # Test for not having only units for values
        invalid_units_1 = OrderedDict([('Watt Per Meter Squared', 'string'), ('Meter Cubed', u.Unit("m3"))])
        assert not sunpy.timeseries.TimeSeries._validate_units(invalid_units_1)
        # Test for being a MetaDict object
        invalid_units_2 = MetaDict(OrderedDict([('Watt Per Meter Squared', u.Unit("W / m2")), ('Meter Cubed', u.Unit("m3"))]))
        assert not sunpy.timeseries.TimeSeries._validate_units(invalid_units_2)

    def test_validate_meta_basic(self):
        valid_meta_1 = MetaDict({'key':'value'})
        assert sunpy.timeseries.TimeSeries._validate_meta(valid_meta_1)
        valid_meta_2 = OrderedDict({'key':'value'})
        assert sunpy.timeseries.TimeSeries._validate_meta(valid_meta_2)
        invalid_meta = []
        assert not sunpy.timeseries.TimeSeries._validate_meta(invalid_meta)

    def test_validate_meta_astropy_header(self):
        # Manually open a goes file for the sunpy.io.header.FileHeader test
        hdus = sunpy.io.read_file(goes_filepath)
        header = hdus[0].header
        assert sunpy.timeseries.TimeSeries._validate_meta(header)
        # Manually open a goes file for the astropy.io.fits.header.Header test
        hdulist = fits.open(goes_filepath)
        header = hdulist[0].header
        hdulist.close()
        assert sunpy.timeseries.TimeSeries._validate_meta(header)

