import numpy as np
import pytest

import astropy.units as u

import sunpy.timeseries
from sunpy.data.test import get_test_filepath
from sunpy.tests.helpers import figure_test

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
    ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath)
    cols_list = [
        '4.0-15.0 keV', '15.0-25.0 keV', '25.0-50.0 keV',
        '50.0-100.0 keV', '100.0-300.0 keV', '300.0-800.0 keV',
        '800.0-2000.0 keV'
    ]
    assert np.all(ts_gbm.to_dataframe().columns == cols_list)


@pytest.mark.parametrize("energy_bins", [
    [45, 95, 150, 270, 500] * u.keV,
    [45, 95, 75, 300, 200] * u.keV,
    [50, 100, 300, 400] * u.J,
    [200 * u.keV, 1000 * u.keV, 1 * u.J, 30 * u.J],
    [2 * u.eV, 5 * u.J, 32 * u.keV, 7 * u.J]
])
def test_fermi_gbm_custom_energy_bins(fermi_gbm_test_ts, energy_bins):
    ts_gbm_base = fermi_gbm_test_ts
    ts_gbm = sunpy.timeseries.TimeSeries(fermi_gbm_filepath, energy_bands=energy_bins)

    ts_gbm_energy_bins = ts_gbm_base.energy_bins(energy_bins)
    ts_gbm_energy_bins__ = ts_gbm.energy_bins(energy_bins)
    energy_bins = sorted(energy_bins)
    col_list = [
        f'{(first.to(u.keV)).value}-{(second.to(u.keV)).value} keV'
        for first, second in zip(energy_bins, energy_bins[1:])
    ]
    assert np.all(ts_gbm.to_dataframe().columns == col_list)
    assert np.all(ts_gbm_energy_bins.to_dataframe().columns == col_list)
    assert np.all(ts_gbm_energy_bins__.to_dataframe().columns == col_list)


@pytest.mark.parametrize("energy_bins", [
    [3 * u.eV],
    []
])
def test_fermi_gbm_invalid_energy_bins(fermi_gbm_test_ts, energy_bins):
    with pytest.raises(UserWarning, match = "'energy_bands' must contain more than one element."):
        sunpy.timeseries.TimeSeries(fermi_gbm_filepath, energy_bands=energy_bins)
    with pytest.raises(UserWarning, match = "'energy_bands' must contain more than one element."):
        fermi_gbm_test_ts.energy_bins(energy_bins)
