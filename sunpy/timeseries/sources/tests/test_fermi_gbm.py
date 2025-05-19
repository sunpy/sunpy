import numpy as np
import pytest

import astropy.units as u

import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test
from sunpy.util import SunpyUserWarning

fermi_gbm_filepath = get_test_filepath('gbm.fits')


def test_implicit_fermi_gbm():
    # Test a GBMSummary TimeSeries
    ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath)
    assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)


def test_fermi_gbm():
    # Test a GBMSummary TimeSeries
    ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath, source='GBMSummary')
    assert isinstance(ts_gbm, sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries)


def test_fermi_gbm_plot_column(fermi_gbm_test_ts):
    ax = fermi_gbm_test_ts.plot(columns=['4.0-15.0 keV', '100.0-300.0 keV'])
    assert len(ax.lines) == 2
    assert '4.0-15.0 keV' == ax.lines[0].get_label()
    assert '100.0-300.0 keV' == ax.lines[1].get_label()


@figure_test
def test_fermi_gbm_peek(fermi_gbm_test_ts):
    fermi_gbm_test_ts.peek()


def test_fermi_gbm_default_energy_bins():
    # Checks for default_ebin, custom_ebin, out-for-bound ebin
    custom_bin = [45, 95, 150, 270, 500] * u.keV
    ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath)
    ts_gbm_1 = ts_gbm.energy_bins(custom_bin)
    ts_gbm_2 = sunpy.timeseries.TimeSeries(fermi_gbm_filepath,
                                           energy_bands = custom_bin)
    with pytest.warns(SunpyUserWarning, match = "Data is not available for"):
        ts_gbm_3 = sunpy.timeseries.TimeSeries(fermi_gbm_filepath, energy_bands = [20, 40, 60, 80, 100] * u.eV)
    with pytest.warns(SunpyUserWarning, match = "Data is not available for"):
        ts_gbm_4 = ts_gbm.energy_bins([20, 40, 60, 80, 100] * u.eV)

    cols_list = [
        '4.0-15.0 keV', '15.0-25.0 keV', '25.0-50.0 keV', '50.0-100.0 keV',
        '100.0-300.0 keV', '300.0-800.0 keV', '800.0-2000.0 keV'
    ]
    cols_list_ = [
        f"{ebin.value}-{next_ebin.value} {ebin.unit}"
        for ebin, next_ebin in zip(custom_bin, custom_bin[1:])
    ]
    assert np.all(ts_gbm.to_dataframe().columns == cols_list)
    assert np.all(ts_gbm_1.to_dataframe().columns == cols_list_)
    assert np.all(ts_gbm_2.to_dataframe().columns == cols_list_)
    assert np.all(ts_gbm_3.to_dataframe().columns == cols_list)
    assert np.all(ts_gbm_4.to_dataframe().columns == cols_list)


@pytest.mark.parametrize("energy_bins", [
    [2 * u.eV, 100 * u.eV, 50 * u.keV, 200 * u.keV],
    [50 * u.keV, 200 * u.keV, 2500 * u.keV],
    [2 * u.eV, 100 * u.eV, 50 * u.keV, 200 * u.keV, 2500 * u.keV]
])
def test_fermi_gbm_custom_energy_bins(fermi_gbm_test_ts, energy_bins):
    # Checks for fully overlapped ebin, partial ebin ranges
    e_min, e_max = 4*u.keV, 2000*u.keV
    ts_gbm_base = fermi_gbm_test_ts
    with pytest.warns(SunpyUserWarning, match = "Only Partially Data is available for"):
        ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath, energy_bands=energy_bins)
    with pytest.warns(SunpyUserWarning, match = "Only Partially Data is available for"):
        ts_gbm_energy_bins = ts_gbm_base.energy_bins(energy_bins)
    with pytest.warns(SunpyUserWarning, match = "Only Partially Data is available for"):
        ts_gbm_energy_bins__ = ts_gbm.energy_bins(energy_bins)

    energy_bins = [ebin.to(u.keV) if ebin.unit != u.keV else ebin for ebin in energy_bins]
    energy_bins = sorted(energy_bins)
    # This only checks for fully overlapped ebin and partial ebin ranges
    if min(energy_bins) < e_min and max(energy_bins) > e_max:
        energy_bins = list(filter(lambda x: e_min <= x <= e_max, energy_bins))
        energy_bins.append(e_max * u.keV)
        energy_bins.insert(0, e_min * u.keV)
    elif min(energy_bins) >= e_min and max(energy_bins) > e_max:
        energy_bins = list(filter(lambda x: x <= e_max, energy_bins))
        energy_bins.append(e_max)
    elif min(energy_bins) < e_min and max(energy_bins) <= e_max:
        energy_bins = list(filter(lambda x: x >= e_min, energy_bins))
        energy_bins.insert(0, e_min * u.keV)

    cols_list_ = [
        f"{ebin.value}-{next_ebin.value} keV"
        for ebin, next_ebin in zip(energy_bins, energy_bins[1:])
    ]
    assert np.all(ts_gbm.to_dataframe().columns == cols_list_)
    assert np.all(ts_gbm_energy_bins.to_dataframe().columns == cols_list_)
    assert np.all(ts_gbm_energy_bins__.to_dataframe().columns == cols_list_)


@pytest.mark.parametrize("energy_bins", [
    [3 * u.eV],
    []
])
def test_fermi_gbm_invalid_energy_bins(fermi_gbm_test_ts, energy_bins):
    with pytest.raises(ValueError, match = "'energy_bands' must contain more than one element."):
        sunpy.timeseries.TimeSeries(fermi_gbm_filepath, energy_bands=energy_bins)
    with pytest.raises(ValueError, match = "'energy_bands' must contain more than one element."):
        fermi_gbm_test_ts.energy_bins(energy_bins)


@pytest.mark.parametrize("energy_bins", [
    [4*u.keV, 10*u.m, 20*u.keV],
    [4*u.keV, 10*u.mag, 20*u.keV]
])
def test_fermi_gbm_invalid_unit_energy_bins(fermi_gbm_test_ts, energy_bins):
    with pytest.raises(u.UnitConversionError, match = r"are not convertible.$"):
        sunpy.timeseries.TimeSeries(fermi_gbm_filepath, energy_bands=energy_bins)
    with pytest.raises(u.UnitConversionError, match = r"are not convertible.$"):
        fermi_gbm_test_ts.energy_bins(energy_bins)
