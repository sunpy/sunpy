import datetime
from collections import OrderedDict

import numpy as np
import pytest
from erfa.core import ErfaWarning
from pandas import DataFrame

import astropy.units as u
from astropy.table import Table
from astropy.time import TimeDelta

import sunpy
import sunpy.timeseries
from sunpy.data.test import get_test_data_filenames, get_test_filepath
from sunpy.time import parse_time
from sunpy.util import SunpyUserWarning
from sunpy.util.metadata import MetaDict

# =============================================================================
# TimeSeries Tests
# =============================================================================

eve_filepath = get_test_filepath('EVE_L0CS_DIODES_1m_truncated.txt')
esp_filepath = get_test_filepath('eve_l1_esp_2011046_00_truncated.fits')
fermi_gbm_filepath = get_test_filepath('gbm.fits')
norh_filepath = get_test_filepath('tca110810_truncated')
lyra_filepath = get_test_filepath('lyra_20150101-000000_lev3_std_truncated.fits.gz')
rhessi_filepath = get_test_filepath('hsi_obssumm_20120601_018_truncated.fits.gz')
noaa_ind_json_filepath = get_test_filepath('observed-solar-cycle-indices-truncated.json')
noaa_pre_json_filepath = get_test_filepath('predicted-solar-cycle-truncated.json')
goes_filepath = get_test_filepath('go1520120601.fits.gz')
a_list_of_many = [f for f in get_test_data_filenames() if f.parents[0].relative_to(f.parents[1]).name == 'eve']


@pytest.fixture
def eve_test_ts():
    with pytest.warns(SunpyUserWarning, match='Unknown units'):
        return sunpy.timeseries.TimeSeries(eve_filepath, source='EVE')


@pytest.fixture
def esp_test_ts():
    return sunpy.timeseries.TimeSeries(esp_filepath, source='ESP')


@pytest.fixture
def fermi_gbm_test_ts():
    return sunpy.timeseries.TimeSeries(fermi_gbm_filepath, source='GBMSummary')


@pytest.fixture
def norh_test_ts():
    return sunpy.timeseries.TimeSeries(norh_filepath, source='NoRH')


@pytest.fixture
def goes_test_ts():
    return sunpy.timeseries.TimeSeries(goes_filepath, source='XRS')


@pytest.fixture
def lyra_test_ts():
    return sunpy.timeseries.TimeSeries(lyra_filepath, source='LYRA')


@pytest.fixture
def rhessi_test_ts():
    return sunpy.timeseries.TimeSeries(rhessi_filepath, source='RHESSI')


@pytest.fixture
def noaa_ind_json_test_ts():
    return sunpy.timeseries.TimeSeries(noaa_ind_json_filepath, source='NOAAIndices')


@pytest.fixture
def noaa_pre_json_test_ts():
    # NOAA pre data contains years long into the future, which ERFA complains about
    with pytest.warns(ErfaWarning, match=r'.*dubious year'):
        return sunpy.timeseries.TimeSeries(
            noaa_pre_json_filepath, source='NOAAPredictIndices')


@pytest.fixture
def generic_ts():
    # Generate the data and the corresponding dates
    base = parse_time("2016/10/01T05:00:00")
    dates = base - TimeDelta(np.arange(24 * 60)*u.minute)
    intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24 * 60))))
    intensity2 = np.cos(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24 * 60))))

    # Create the data DataFrame, header MetaDict and units dict
    data = DataFrame(np.column_stack([intensity, intensity2]),
                     index=dates.isot.astype('datetime64'),
                     columns=['intensity', 'intensity2'])
    units = {'intensity': u.W / u.m**2,
             'intensity2': u.W / u.m**2}
    meta = MetaDict({'key': 'value'})

    # Create the time series
    return sunpy.timeseries.TimeSeries(data, meta, units)


@pytest.fixture
def concatenate_multi_files_ts():
    with pytest.warns(SunpyUserWarning, match='Unknown units'):
        return sunpy.timeseries.TimeSeries(
            a_list_of_many, source='EVE', concatenate=True)

# =============================================================================
# Test Creating TimeSeries From Various Dataformats
# =============================================================================


@pytest.fixture
def table_ts():
    # Generate the data and the corresponding dates
    base = parse_time(datetime.datetime.today())
    times = base - TimeDelta(np.arange(24 * 60)*u.minute)
    intensity = u.Quantity(
        np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24 * 60)))), u.W / u.m ** 2)

    # Create the units and meta objects
    units = OrderedDict([('intensity', u.W / u.m**2)])
    meta = MetaDict({'key': 'value'})
    tbl_meta = MetaDict({'t_key': 't_value'})

    # Create a suitable mixin qtable
    table = Table(
        [times, intensity], names=['time', 'intensity'], meta=tbl_meta)
    table.add_index('time')

    # Create TS from dataframe and check
    return sunpy.timeseries.TimeSeries(table, meta, units)


# =============================================================================
# Test Resulting TimeSeries Parameters
# =============================================================================
@pytest.fixture(params=['eve_test_ts', 'esp_test_ts', 'fermi_gbm_test_ts', 'norh_test_ts', 'goes_test_ts',
                        'lyra_test_ts', 'rhessi_test_ts', 'noaa_ind_json_test_ts',
                        'noaa_pre_json_test_ts', 'generic_ts', 'table_ts'])
def many_ts(request):
    # Fixture to return lots of different timeseries
    return request.getfixturevalue(request.param)
