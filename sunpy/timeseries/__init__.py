# Check if user has installed the timeseries extras
from sunpy.util.sysinfo import _warn_missing_deps

_warn_missing_deps('timeseries')

from sunpy.timeseries.metadata import TimeSeriesMetaData
from sunpy.timeseries.sources import *
from sunpy.timeseries.timeseries_factory import TimeSeries
from sunpy.timeseries.timeseriesbase import GenericTimeSeries

try:
    # Register pandas datetime converter with matplotlib
    import pandas.plotting
    pandas.plotting.register_matplotlib_converters()
except ImportError:
    pass

__all__ = ["TimeSeriesMetaData", "TimeSeries", "GenericTimeSeries"]
