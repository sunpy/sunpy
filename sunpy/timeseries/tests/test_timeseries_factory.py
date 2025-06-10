import os
import logging
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

import sunpy.io
import sunpy.net.attrs as a
import sunpy.timeseries
from sunpy.data.test import get_test_data_filenames, get_test_filepath, rootdir
from sunpy.net import Fido
from sunpy.time import parse_time
from sunpy.util import SunpyUserWarning
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.metadata import MetaDict

eve_filepath = get_test_filepath('EVE_L0CS_DIODES_1m_truncated.txt')
eve_many_filepath = [f for f in get_test_data_filenames()
                     if f.parents[0].relative_to(f.parents[1]).name == 'eve']
goes_filepath = get_test_filepath('go1520110607.fits')
psp_filepath = get_test_filepath('psp_fld_l2_mag_rtn_1min_20200104_v02.cdf')
swa_filepath = get_test_filepath('solo_L1_swa-pas-mom_20200706_V01.cdf')
fermi_gbm_filepath = get_test_filepath('gbm.fits')
hsi_filepath = get_test_filepath('hsi_image_20101016_191218.fits')


@pytest.mark.filterwarnings('ignore:Unknown units')
def test_factory_concatenate_same_source():
    # Test making a TimeSeries that is the concatenation of multiple files
    ts_from_list = sunpy.timeseries.TimeSeries(eve_many_filepath, source='EVE', concatenate=True)
    assert isinstance(ts_from_list, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)

    ts_from_folder = sunpy.timeseries.TimeSeries(
        eve_many_filepath[0].parent, source='EVE', concatenate=True)
    assert isinstance(ts_from_folder, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
    # text the two methods get identical dataframes
    assert ts_from_list == ts_from_folder
    # test the frames have correct headings/keys (correct concatenation axis)
    ts_from_list.columns == sunpy.timeseries.TimeSeries(
        eve_many_filepath[0], source='EVE', concatenate=True).columns


@pytest.mark.filterwarnings('ignore:Unknown units')
def test_factory_concatenate_different_source():
    # Test making a TimeSeries that is the concatenation of multiple files
    ts_from_list = sunpy.timeseries.TimeSeries(eve_many_filepath, source='EVE', concatenate=True)
    assert isinstance(ts_from_list, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
    ts_from_folder = sunpy.timeseries.TimeSeries(
        eve_many_filepath[0].parent, source='EVE', concatenate=True)
    assert isinstance(ts_from_folder, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)
    # text the two methods get identical dataframes
    assert ts_from_list == ts_from_folder
    # test the frames have correct headings/keys (correct concatenation axis)
    ts_from_list.columns == sunpy.timeseries.TimeSeries(
        eve_many_filepath[0], source='EVE', concatenate=True).columns


@pytest.mark.filterwarnings('ignore:Unknown units')
def test_factory_generate_list_of_ts():
    # Test making a list TimeSeries from multiple files
    ts_list = sunpy.timeseries.TimeSeries(eve_many_filepath, source='EVE')
    assert isinstance(ts_list, list)
    for ts in ts_list:
        assert isinstance(ts, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)


@pytest.mark.filterwarnings('ignore:Unknown units')
def test_factory_generate_from_glob():
    # Test making a TimeSeries from a glob
    ts_from_glob = sunpy.timeseries.TimeSeries(os.path.join(
        rootdir, "eve", "*"), source='EVE', concatenate=True)
    assert isinstance(ts_from_glob, sunpy.timeseries.sources.eve.EVESpWxTimeSeries)


@pytest.mark.filterwarnings('ignore:Unknown units')
def test_factory_generate_from_pathlib():
    # Test making a TimeSeries from a : pathlib.PosixPath
    ts_from_pathlib = sunpy.timeseries.TimeSeries(Path(fermi_gbm_filepath),
                                                  source="GBMSummary")
    assert isinstance(ts_from_pathlib, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)


@pytest.mark.remote_data
def test_from_url():
    # This is the same PSP file we have in our test data, but accessed from a URL
    url = ('https://spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/mag_rtn_1min/2020/'
           'psp_fld_l2_mag_rtn_1min_20200104_v02.cdf')
    ts = sunpy.timeseries.TimeSeries(url)
    assert isinstance(ts[0], sunpy.timeseries.GenericTimeSeries)
    assert isinstance(ts[1], sunpy.timeseries.GenericTimeSeries)

@pytest.mark.remote_data
def test_from_uri():
    # Test read on PSP file saved on public sumpy s3 repository.
    uri = ('s3://data.sunpy.org/sunpy/v1/psp_fld_l2_mag_rtn_1min_20200104_v02.cdf')
    ts = sunpy.timeseries.TimeSeries(uri, fsspec_kwargs={'anon':True})
    assert isinstance(ts[0], sunpy.timeseries.GenericTimeSeries)
    assert isinstance(ts[1], sunpy.timeseries.GenericTimeSeries)

def test_read_cdf():
    ts_psp = sunpy.timeseries.TimeSeries(psp_filepath)
    assert len(ts_psp) == 2

    ts = ts_psp[0]
    assert ts.columns == ['psp_fld_l2_mag_RTN_1min_0',
                          'psp_fld_l2_mag_RTN_1min_1',
                          'psp_fld_l2_mag_RTN_1min_2']
    assert ts.quantity('psp_fld_l2_mag_RTN_1min_0').unit == u.nT
    assert len(ts.quantity('psp_fld_l2_mag_RTN_1min_0')) == 118

    ts = ts_psp[1]
    assert ts.columns == ['psp_fld_l2_quality_flags']
    assert ts.quantity('psp_fld_l2_quality_flags').unit == u.dimensionless_unscaled
    assert len(ts.quantity('psp_fld_l2_quality_flags')) == 1440


@pytest.mark.remote_data
def test_read_cdf_empty_variable():
    # This tests that:
    # - A CDF file with an empty column can be read
    # - Unknown unit handling works as expected
    result = sunpy.net.Fido.search(a.Time('2020-01-01', '2020-01-02'),
                                   a.cdaweb.Dataset('AC_H6_SWI'))
    filename = Fido.fetch(result[0, 0])

    # Temporarily reset sunpy.io._cdf registry of known unit conversions
    import sunpy.io._cdf as sunpy_cdf
    known_units = sunpy_cdf._known_units
    sunpy_cdf._known_units = {}

    with pytest.warns(SunpyUserWarning, match='Assigning dimensionless units'):
        ts = sunpy.timeseries.TimeSeries(filename)

    assert ts.quantity('nH').unit == u.dimensionless_unscaled

    # Put back known unit registry, and check that units are recognised
    sunpy_cdf._known_units = known_units
    ts = sunpy.timeseries.TimeSeries(filename)
    assert ts.quantity('nH').unit == u.cm**-3

    # Reset again to check that registering units via. astropy works too
    sunpy_cdf._known_units = {}
    u.add_enabled_units([u.def_unit('#/cm^3', represents=u.cm**-3)])
    ts = sunpy.timeseries.TimeSeries(filename)
    assert ts.quantity('nH').unit == u.cm**-3

    sunpy_cdf._known_units = known_units


def test_read_empty_cdf(caplog):
    with caplog.at_level(logging.DEBUG, logger='sunpy'):
        ts_empty = sunpy.timeseries.TimeSeries(swa_filepath)
    assert ts_empty == []
    assert "No data found in file" in caplog.text
    assert "solo_L1_swa-pas-mom_20200706_V01.cdf" in caplog.text


def test_meta_from_fits_header():
    # Generate the data and the corresponding dates
    base = parse_time(datetime.datetime.today())
    times = base - TimeDelta(np.arange(24*60)*u.minute)
    intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
    units = {'intensity': u.W/u.m**2}
    data = DataFrame(intensity, index=times, columns=['intensity'])

    # Use a FITS file HDU using sunpy.io
    hdulist = sunpy.io._file_tools.read_file(goes_filepath)
    meta = hdulist[0].header
    meta_md = MetaDict(OrderedDict(meta))
    ts_hdu_meta = sunpy.timeseries.TimeSeries(data, meta, units)
    ts_md_meta = sunpy.timeseries.TimeSeries(data, meta_md, units)
    assert ts_hdu_meta == ts_md_meta

    # Use a FITS file HDU using astropy.io
    hdulist = fits.open(goes_filepath)
    meta = hdulist[0].header
    hdulist.close()
    meta_md = MetaDict(sunpy.io._header.FileHeader(meta))
    ts_hdu_meta = sunpy.timeseries.TimeSeries(data, meta, units)
    ts_md_meta = sunpy.timeseries.TimeSeries(data, meta_md, units)
    assert ts_hdu_meta == ts_md_meta


def test_generic_construction_basic():
    # Generate the data and the corresponding dates
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


def test_generic_construction_basic_omitted_details():
    # Generate the data and the corresponding dates
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


def test_generic_construction_basic_different_meta_types():
    # Generate the data and the corresponding dates
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


def test_generic_construction_ts_list():
    # Generate the data and the corresponding dates
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


def test_generic_construction_concatenation():
    # Generate the data and the corresponding dates
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

    # Concatenate during construction
    ts_concat_2 = sunpy.timeseries.TimeSeries(
        data, meta, units, data2, meta2, units2, concatenate=True)
    assert isinstance(ts_concat_2, sunpy.timeseries.timeseriesbase.GenericTimeSeries)

    # Create TS using a tuple
    ts_concat_3 = sunpy.timeseries.TimeSeries(
        ((data, meta, units), (data2, meta2, units2)), concatenate=True)
    assert isinstance(ts_concat_3, sunpy.timeseries.timeseriesbase.GenericTimeSeries)
    assert ts_concat_1 == ts_concat_2 == ts_concat_3


def test_table_to_ts():
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
    with pytest.raises(ValueError, match="Invalid input Table, TimeSeries doesn't support conversion of tables with more then one index column."):
        sunpy.timeseries.TimeSeries((dual_index_table, meta, units))


def test_passed_ts():
    # Test an EVE TimeSeries
    with pytest.warns(SunpyUserWarning, match='Unknown units'):
        ts_eve = sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')
    ts_from_ts_1 = sunpy.timeseries.TimeSeries(ts_eve, source='EVE')
    ts_from_ts_2 = sunpy.timeseries.TimeSeries(ts_eve)
    assert ts_eve == ts_from_ts_1 == ts_from_ts_2


def test_invalid_manual_data():
    meta = MetaDict({'key': 'value'})
    data = []
    with pytest.raises(NoMatchError, match="One of the files failed to validate with"):
        sunpy.timeseries.TimeSeries(data, meta)


@pytest.mark.filterwarnings('ignore:"silence_errors" was deprecated in version 5')
def test_invalid_filepath():
    invalid_filepath = os.path.join(rootdir, 'invalid_filepath_here')
    with pytest.raises(ValueError, match='Did not find any files'):
        sunpy.timeseries.TimeSeries(invalid_filepath)
    # Now with silence_errors kwarg set
    with pytest.raises(ValueError, match='Did not find any files'):
        sunpy.timeseries.TimeSeries(invalid_filepath, silence_errors=True)
    # Now with allow_errors kwarg set
    with pytest.raises(ValueError, match='Did not find any files'):
        sunpy.timeseries.TimeSeries(invalid_filepath, allow_errors=True)


@pytest.mark.filterwarnings('ignore:"silence_errors" was deprecated in version 5')
def test_invalid_file():
    invalid_filepath = os.path.join(rootdir, 'annotation_ppt.db')
    with pytest.raises(NoMatchError):
        sunpy.timeseries.TimeSeries(invalid_filepath)
    # Now with silence_errors kwarg set
    with pytest.warns(SunpyUserWarning, match="One of the files failed to validate with: Could not find any timeseries sources to parse"):
        ts = sunpy.timeseries.TimeSeries(invalid_filepath, silence_errors=True)
    assert ts == []
    # Now with allow_errors kwarg set
    with pytest.warns(SunpyUserWarning, match="One of the files failed to validate with: Could not find any timeseries sources to parse"):
        ts = sunpy.timeseries.TimeSeries(invalid_filepath, allow_errors=True)
    assert ts == []


def test_validate_units():
    valid_units = OrderedDict(
        [('Watt Per Meter Squared', u.Unit("W / m2")), ('Meter Cubed', u.Unit("m3"))])
    assert sunpy.timeseries.TimeSeries._is_units(valid_units)
    # Test for not having only units for values
    invalid_units_1 = OrderedDict(
        [('Watt Per Meter Squared', 'string'), ('Meter Cubed', u.Unit("m3"))])
    assert not sunpy.timeseries.TimeSeries._is_units(invalid_units_1)
    # Test for being a MetaDict object
    invalid_units_2 = MetaDict(OrderedDict(
        [('Watt Per Meter Squared', u.Unit("W / m2")), ('Meter Cubed', u.Unit("m3"))]))
    assert not sunpy.timeseries.TimeSeries._is_units(invalid_units_2)


def test_validate_meta_basic():
    valid_meta_1 = MetaDict({'key': 'value'})
    assert sunpy.timeseries.TimeSeries._is_metadata(valid_meta_1)
    valid_meta_2 = OrderedDict({'key': 'value'})
    assert sunpy.timeseries.TimeSeries._is_metadata(valid_meta_2)
    time_range = sunpy.time.TimeRange('2020-01-01 12:00', '2020-01-02 12:00')
    valid_meta_3 = sunpy.timeseries.TimeSeriesMetaData(time_range)
    assert sunpy.timeseries.TimeSeries._is_metadata(valid_meta_3)
    invalid_meta = []
    assert not sunpy.timeseries.TimeSeries._is_metadata(invalid_meta)


def test_validate_meta_astropy_header():
    # Manually open a goes file for the sunpy.io._header.FileHeader test
    hdus = sunpy.io._file_tools.read_file(goes_filepath)
    header = hdus[0].header
    assert sunpy.timeseries.TimeSeries._is_metadata(header)
    # Manually open a goes file for the astropy.io.fits.header.Header test
    hdulist = fits.open(goes_filepath)
    header = hdulist[0].header
    hdulist.close()
    assert sunpy.timeseries.TimeSeries._is_metadata(header)

def test_get_matching_widget():
    with pytest.raises(NoMatchError, match="failed to validate"):
        sunpy.timeseries.TimeSeries(hsi_filepath)
