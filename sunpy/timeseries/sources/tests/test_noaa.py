import pytest

import sunpy.timeseries
from sunpy.data.test import get_test_filepath

noaa_ind_json_filepath = get_test_filepath('observed-solar-cycle-indices-truncated.json')
noaa_pre_json_filepath = get_test_filepath('predicted-solar-cycle-truncated.json')


def test_noaa_ind_json():
    # Test a NOAAPredictIndices TimeSeries json
    ts_noaa_ind = sunpy.timeseries.TimeSeries(noaa_ind_json_filepath, source='NOAAIndices')
    assert isinstance(ts_noaa_ind, sunpy.timeseries.sources.noaa.NOAAIndicesTimeSeries)


# The pre- data involves dates long in the future, so ignore an ERFA warning
# when parsing these dates.
@pytest.mark.filterwarnings('ignore:ERFA function.*dubious year')
def test_noaa_pre_json():
    # Test a NOAAIndices TimeSeries json
    ts_noaa_pre = sunpy.timeseries.TimeSeries(
        noaa_pre_json_filepath, source='NOAAPredictIndices')
    assert isinstance(ts_noaa_pre, sunpy.timeseries.sources.noaa.NOAAPredictIndicesTimeSeries)
