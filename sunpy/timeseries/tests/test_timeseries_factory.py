import os
import glob
import datetime
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pytest
from pandas import DataFrame

import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.time import TimeDelta

import sunpy.data.test
import sunpy.io
import sunpy.timeseries
from sunpy.time import parse_time
from sunpy.util import SunpyUserWarning
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.exceptions import SunpyDeprecationWarning
from sunpy.util.metadata import MetaDict

# =============================================================================
# TimeSeries Factory Tests
# =============================================================================

filepath = sunpy.data.test.rootdir
eve_filepath = os.path.join(filepath, 'EVE_L0CS_DIODES_1m_truncated.txt')
esp_filepath = os.path.join(filepath, 'eve_l1_esp_2011046_00_truncated.fits')
fermi_gbm_filepath = os.path.join(filepath, 'gbm.fits')
norh_filepath = os.path.join(filepath, 'tca110810_truncated')
lyra_filepath = os.path.join(filepath, 'lyra_20150101-000000_lev3_std_truncated.fits.gz')
rhessi_filepath = os.path.join(filepath, 'hsi_obssumm_20120601_018_truncated.fits.gz')
noaa_ind_json_filepath = os.path.join(filepath, 'observed-solar-cycle-indices-truncated.json')
noaa_pre_json_filepath = os.path.join(filepath, 'predicted-solar-cycle-truncated.json')
noaa_ind_txt_filepath = os.path.join(filepath, 'RecentIndices_truncated.txt')
noaa_pre_txt_filepath = os.path.join(filepath, 'predicted-sunspot-radio-flux_truncated.txt')
goes_filepath_com = os.path.join(filepath, 'go1520120601.fits.gz')
goes_filepath = os.path.join(filepath, 'go1520110607.fits')
new_goes15_filepath = os.path.join(filepath, 'goes_truncated_test_goes15.nc')
new_goes16_filepath = os.path.join(filepath, 'goes_truncated_test_goes16.nc')
a_list_of_many = glob.glob(os.path.join(filepath, "eve", "*"))

# =============================================================================
# Multi file Tests
# =============================================================================


class TestTimeSeries:
    @pytest.mark.filterwarnings('ignore:Unknown units')
    def test_factory_concatenate_same_source(self):
        # Test making a TimeSeries that is the concatenation of multiple files
        ts_from_list = sunpy.timeseries.TimeSeries(a_list_of_many, source='EVE', concatenate=True)
        assert isinstance(ts_from_list, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)

        ts_from_folder = sunpy.timeseries.TimeSeries(
            os.path.join(filepath, "eve"), source='EVE', concatenate=True)
        assert isinstance(ts_from_folder, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
        # text the two methods get identical dataframes
        assert ts_from_list == ts_from_folder
        # test the frames have correct headings/keys (correct concatenation axis)
        ts_from_list.columns == sunpy.timeseries.TimeSeries(
            a_list_of_many[0], source='EVE', concatenate=True).columns

    @pytest.mark.filterwarnings('ignore:Unknown units')
    def test_factory_concatenate_different_source(self):
        # Test making a TimeSeries that is the concatenation of multiple files
        ts_from_list = sunpy.timeseries.TimeSeries(a_list_of_many, source='EVE', concatenate=True)
        assert isinstance(ts_from_list, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
        ts_from_folder = sunpy.timeseries.TimeSeries(
            os.path.join(filepath, "eve"), source='EVE', concatenate=True)
        assert isinstance(ts_from_folder, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
        # text the two methods get identical dataframes
        assert ts_from_list == ts_from_folder
        # test the frames have correct headings/keys (correct concatenation axis)
        ts_from_list.columns == sunpy.timeseries.TimeSeries(
            a_list_of_many[0], source='EVE', concatenate=True).columns

    @pytest.mark.filterwarnings('ignore:Unknown units')
    def test_factory_generate_list_of_ts(self):
        # Test making a list TimeSeries from multiple files
        ts_list = sunpy.timeseries.TimeSeries(a_list_of_many, source='EVE')
        assert isinstance(ts_list, list)
        for ts in ts_list:
            assert isinstance(ts, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)

    @pytest.mark.filterwarnings('ignore:Unknown units')
    def test_factory_generate_from_glob(self):
        # Test making a TimeSeries from a glob
        ts_from_glob = sunpy.timeseries.TimeSeries(os.path.join(
            filepath, "eve", "*"), source='EVE', concatenate=True)
        assert isinstance(ts_from_glob, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)

    @pytest.mark.filterwarnings('ignore:Unknown units')
    def test_factory_generate_from_pathlib(self):
        # Test making a TimeSeries from a : pathlib.PosixPath
        ts_from_pathlib = sunpy.timeseries.TimeSeries(Path(filepath).joinpath("gbm.fits"),
                                                      source="GBMSummary")
        assert isinstance(ts_from_pathlib, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)
# =============================================================================
# Individual Implicit Source Tests
# =============================================================================

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

    def test_implicit_new_goes15(self):
        # Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(new_goes15_filepath)
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)

    def test_implicit_new_goes16(self):
        # Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(new_goes16_filepath)
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)

    def test_implicit_lyra(self):
        # Test a LYRA TimeSeries
        ts_lyra = sunpy.timeseries.TimeSeries(lyra_filepath)
        assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRATimeSeries)

    def test_implicit_rhessi(self):
        # Test a RHESSI TimeSeries
        ts_rhessi = sunpy.timeseries.TimeSeries(rhessi_filepath)
        assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryTimeSeries)

    def test_implicit_esp(self):
        # Test an ESP TimeSeries
        ts_esp = sunpy.timeseries.TimeSeries(esp_filepath)
        assert isinstance(ts_esp, sunpy.timeseries.sources.eve.ESPTimeSeries)

# =============================================================================
# Individual Explicit Sources Tests
# =============================================================================

    def test_eve(self):
        # Test an EVE TimeSeries
        with pytest.warns(SunpyUserWarning, match='Unknown units for x_cool proxy'):
            ts_eve = sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')
        assert isinstance(ts_eve, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)

    def test_esp(self):
        # Test an ESP TimeSeries
        ts_esp = sunpy.timeseries.TimeSeries(esp_filepath, source='ESP')
        assert isinstance(ts_esp, sunpy.timeseries.sources.eve.ESPTimeSeries)

    def test_fermi_gbm(self):
        # Test a GBMSummary TimeSeries
        ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath, source='GBMSummary')
        assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)

    def test_norh(self):
        # Test a NoRH TimeSeries
        ts_norh = sunpy.timeseries.TimeSeries(norh_filepath, source='NoRH')
        assert isinstance(ts_norh, sunpy.timeseries.sources.norh.NoRHTimeSeries)

    def test_goes(self):
        # Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(goes_filepath, source='XRS')
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)

    def test_goes_com(self):
        # Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(goes_filepath_com, source='XRS')
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)

    def test_new_goes15(self):
        # Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(new_goes15_filepath, source='XRS')
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)

    def test_new_goes16(self):
        # Test a GOES TimeSeries
        ts_goes = sunpy.timeseries.TimeSeries(new_goes16_filepath, source='XRS')
        assert isinstance(ts_goes, sunpy.timeseries.sources.goes.XRSTimeSeries)

    def test_lyra(self):
        # Test a LYRA TimeSeries
        ts_lyra = sunpy.timeseries.TimeSeries(lyra_filepath, source='LYRA')
        assert isinstance(ts_lyra, sunpy.timeseries.sources.lyra.LYRATimeSeries)

    def test_rhessi(self):
        # Test a RHESSI TimeSeries
        ts_rhessi = sunpy.timeseries.TimeSeries(rhessi_filepath, source='RHESSI')
        assert isinstance(ts_rhessi, sunpy.timeseries.sources.rhessi.RHESSISummaryTimeSeries)

    def test_noaa_ind_json(self):
        # Test a NOAAPredictIndices TimeSeries json
        ts_noaa_ind = sunpy.timeseries.TimeSeries(noaa_ind_json_filepath, source='NOAAIndices')
        assert isinstance(ts_noaa_ind, sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries)

    def test_noaa_ind_txt(self):
        # Test a NOAAPredictIndices TimeSeries txt
        with pytest.warns(SunpyDeprecationWarning):
            ts_noaa_ind = sunpy.timeseries.TimeSeries(noaa_ind_txt_filepath, source='NOAAIndices')
        assert isinstance(ts_noaa_ind, sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries)

    # The pre- data involves dates long in the future, so ignore an ERFA warning
    # when parsing these dates.
    @pytest.mark.filterwarnings('ignore:ERFA function.*dubious year')
    def test_noaa_pre_json(self):
        # Test a NOAAIndices TimeSeries json
        ts_noaa_pre = sunpy.timeseries.TimeSeries(
            noaa_pre_json_filepath, source='NOAAPredictIndices')
        assert isinstance(ts_noaa_pre, sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries)

    def test_noaa_pre_txt(self):
        # Test a NOAAIndices TimeSeries txt
        with pytest.warns(SunpyDeprecationWarning):
            ts_noaa_pre = sunpy.timeseries.TimeSeries(
                noaa_pre_txt_filepath, source='NOAAPredictIndices')
        assert isinstance(ts_noaa_pre, sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries)

# ==============================================================================
# Remote Sources Tests
# ==============================================================================

    @pytest.mark.remote_data
    def test_goes_remote(self):
        # Older format file
        goes = sunpy.timeseries.TimeSeries(
            'https://umbra.nascom.nasa.gov/goes/fits/1986/go06860129.fits')
        assert isinstance(goes, sunpy.timeseries.sources.goes.XRSTimeSeries)
        # Newer format
        goes = sunpy.timeseries.TimeSeries(
            'https://umbra.nascom.nasa.gov/goes/fits/2018/go1520180626.fits')
        assert isinstance(goes, sunpy.timeseries.sources.goes.XRSTimeSeries)

# =============================================================================
# Manual TimeSeries Tests
# =============================================================================

    def test_meta_from_fits_header(self):
        # Generate the data and the corrisponding dates
        base = parse_time(datetime.datetime.today())
        times = base - TimeDelta(np.arange(24*60)*u.minute)
        intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
        units = {'intensity': u.W/u.m**2}
        data = DataFrame(intensity, index=times, columns=['intensity'])

        # Use a FITS file HDU using sunpy.io
        hdulist = sunpy.io.read_file(goes_filepath)
        meta = hdulist[0].header
        meta_md = MetaDict(OrderedDict(meta))
        ts_hdu_meta = sunpy.timeseries.TimeSeries(data, meta, units)
        ts_md_meta = sunpy.timeseries.TimeSeries(data, meta_md, units)
        assert ts_hdu_meta == ts_md_meta

        # Use a FITS file HDU using astropy.io
        hdulist = fits.open(goes_filepath)
        meta = hdulist[0].header
        hdulist.close()
        meta_md = MetaDict(sunpy.io.header.FileHeader(meta))
        ts_hdu_meta = sunpy.timeseries.TimeSeries(data, meta, units)
        ts_md_meta = sunpy.timeseries.TimeSeries(data, meta_md, units)
        assert ts_hdu_meta == ts_md_meta

    def test_generic_construction_basic(self):
        # Generate the data and the corrisponding dates
        base = parse_time(datetime.datetime.today())
        times = base - TimeDelta(np.arange(24 * 60)*u.minute)
        intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity, index=times, columns=['intensity'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        meta = MetaDict({'key': 'value'})

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
        base = parse_time(datetime.datetime.today())
        times = base - TimeDelta(np.arange(24 * 60)*u.minute)
        intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity, index=times, columns=['intensity'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        meta = MetaDict({'key': 'value'})

        # Create TS omitting units input arguments
        with pytest.warns(SunpyUserWarning, match='Unknown units for intensity'):
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
        base = parse_time(datetime.datetime.today())
        times = base - TimeDelta(np.arange(24 * 60)*u.minute)
        intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
        tr = sunpy.time.TimeRange(times[0], times[-1])

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity, index=times, columns=['intensity'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        meta_md = MetaDict({'key': 'value'})
        meta_di = {'key': 'value'}
        meta_od = OrderedDict({'key': 'value'})
        meta_obj = sunpy.timeseries.TimeSeriesMetaData(timerange=tr, colnames=['GOES'],
                                                       meta=MetaDict({'key': 'value'}))

        # Create TS using different dictionary meta types
        ts_md = sunpy.timeseries.TimeSeries(data, meta_md, units)
        ts_di = sunpy.timeseries.TimeSeries(data, meta_di, units)
        ts_od = sunpy.timeseries.TimeSeries(data, meta_od, units)
        ts_obj = sunpy.timeseries.TimeSeries(data, meta_obj, units)
        assert ts_md == ts_di == ts_od == ts_obj
        assert ts_md.meta.metadata[0][2] == ts_di.meta.metadata[0][2] == ts_od.meta.metadata[0][2] == ts_obj.meta.metadata[0][2]

    def test_generic_construction_ts_list(self):
        # Generate the data and the corrisponding dates
        base = parse_time(datetime.datetime.today())
        times = base - TimeDelta(np.arange(24 * 60)*u.minute)
        intensity1 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
        intensity2 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity1, index=times, columns=['intensity'])
        data2 = DataFrame(intensity2, index=times, columns=['intensity2'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        units2 = OrderedDict([('intensity2', u.W/u.m**2)])
        meta = MetaDict({'key': 'value'})
        meta2 = MetaDict({'key2': 'value2'})

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
        ts_list2 = sunpy.timeseries.TimeSeries(((data, meta, units), (data2, meta2, units2)))
        assert ts_list == ts_list2

    def test_generic_construction_concatenation(self):
        # Generate the data and the corrisponding dates
        base = parse_time(datetime.datetime.today())
        times = base - TimeDelta(np.arange(24 * 60)*u.minute)
        intensity1 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
        intensity2 = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))

        # Create the data DataFrame, header MetaDict and units OrderedDict
        data = DataFrame(intensity1, index=times, columns=['intensity'])
        data2 = DataFrame(intensity2, index=times, columns=['intensity2'])
        units = OrderedDict([('intensity', u.W/u.m**2)])
        units2 = OrderedDict([('intensity2', u.W/u.m**2)])
        meta = MetaDict({'key': 'value'})
        meta2 = MetaDict({'key2': 'value2'})

        # Create TS individually
        ts_1 = sunpy.timeseries.TimeSeries(data, meta, units)
        ts_2 = sunpy.timeseries.TimeSeries(data2, meta2, units2)
        ts_concat_1 = ts_1.concatenate(ts_2)

        # Concatinate during construction
        ts_concat_2 = sunpy.timeseries.TimeSeries(
            data, meta, units, data2, meta2, units2, concatenate=True)
        assert isinstance(ts_concat_2, sunpy.timeseries.timeseriesbase.GenericTimeSeries)

        # Create TS using a tuple
        ts_concat_3 = sunpy.timeseries.TimeSeries(
            ((data, meta, units), (data2, meta2, units2)), concatenate=True)
        assert isinstance(ts_concat_3, sunpy.timeseries.timeseriesbase.GenericTimeSeries)
        assert ts_concat_1 == ts_concat_2 == ts_concat_3

    def test_table_to_ts(self):
        # Generate the data and the corresponding dates
        base = parse_time(datetime.datetime.today())
        times = base - TimeDelta(np.arange(24 * 60)*u.minute)
        intensity = u.Quantity(
            np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60)))), u.W/u.m**2)

        # Create the units and meta objects
        units = OrderedDict([('intensity', u.W/u.m**2)])
        meta = MetaDict({'key': 'value'})
        tbl_meta = MetaDict({'t_key': 't_value'})

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

# =============================================================================
# Test some other options
# =============================================================================

    def test_passed_ts(self):
        # Test an EVE TimeSeries
        with pytest.warns(SunpyUserWarning, match='Unknown units'):
            ts_eve = sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')
        ts_from_ts_1 = sunpy.timeseries.TimeSeries(ts_eve, source='EVE')
        ts_from_ts_2 = sunpy.timeseries.TimeSeries(ts_eve)
        assert ts_eve == ts_from_ts_1 == ts_from_ts_2

# =============================================================================
# Test some Errors
# =============================================================================

    def test_invalid_manual_data(self):
        meta = MetaDict({'key': 'value'})
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
        valid_units = OrderedDict(
            [('Watt Per Meter Squared', u.Unit("W / m2")), ('Meter Cubed', u.Unit("m3"))])
        assert sunpy.timeseries.TimeSeries._validate_units(valid_units)
        # Test for not having only units for values
        invalid_units_1 = OrderedDict(
            [('Watt Per Meter Squared', 'string'), ('Meter Cubed', u.Unit("m3"))])
        assert not sunpy.timeseries.TimeSeries._validate_units(invalid_units_1)
        # Test for being a MetaDict object
        invalid_units_2 = MetaDict(OrderedDict(
            [('Watt Per Meter Squared', u.Unit("W / m2")), ('Meter Cubed', u.Unit("m3"))]))
        assert not sunpy.timeseries.TimeSeries._validate_units(invalid_units_2)

    def test_validate_meta_basic(self):
        valid_meta_1 = MetaDict({'key': 'value'})
        assert sunpy.timeseries.TimeSeries._validate_meta(valid_meta_1)
        valid_meta_2 = OrderedDict({'key': 'value'})
        assert sunpy.timeseries.TimeSeries._validate_meta(valid_meta_2)
        time_range = sunpy.time.TimeRange('2020-01-01 12:00', '2020-01-02 12:00')
        valid_meta_3 = sunpy.timeseries.TimeSeriesMetaData(time_range)
        assert sunpy.timeseries.TimeSeries._validate_meta(valid_meta_3)
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
