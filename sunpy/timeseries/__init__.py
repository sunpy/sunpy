from sunpy.timeseries.metadata import TimeSeriesMetaData
from sunpy.timeseries.timeseries_factory import TimeSeries
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.timeseries.sources.eve import EVESpWxTimeSeries
from sunpy.timeseries.sources.goes import XRSTimeSeries
from sunpy.timeseries.sources.noaa import NOAAIndicesTimeSeries, NOAAPredictIndicesTimeSeries
from sunpy.timeseries.sources.lyra import LYRATimeSeries
from sunpy.timeseries.sources.norh import NoRHTimeSeries
from sunpy.timeseries.sources.rhessi import RHESSISummaryTimeSeries
from sunpy.timeseries.sources.fermi_gbm import GBMSummaryTimeSeries

try:
    # Register pandas datetime converter with matplotlib
    import pandas.plotting
    pandas.plotting.register_matplotlib_converters()
except ImportError:
    pass

__all__ = ["TimeSeriesMetaData", "TimeSeries", "GenericTimeSeries"]
