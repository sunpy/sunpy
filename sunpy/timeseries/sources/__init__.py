"""
This module provides a collection of datasource-specific
`~sunpy.timeseries.TimeSeries` classes.

Each mission should have its own file with one or more classes defined.
Typically, these classes will be subclasses of the
`sunpy.timeseries.TimeSeries`.
"""
from sunpy.timeseries.sources.eve import ESPTimeSeries, EVESpWxTimeSeries
from sunpy.timeseries.sources.fermi_gbm import GBMSummaryTimeSeries
from sunpy.timeseries.sources.goes import XRSTimeSeries
from sunpy.timeseries.sources.lyra import LYRATimeSeries
from sunpy.timeseries.sources.noaa import NOAAIndicesTimeSeries, NOAAPredictIndicesTimeSeries
from sunpy.timeseries.sources.norh import NoRHTimeSeries
from sunpy.timeseries.sources.rhessi import RHESSISummaryTimeSeries

__all__ = ['ESPTimeSeries', 'EVESpWxTimeSeries', 'GBMSummaryTimeSeries', 'XRSTimeSeries', 'LYRATimeSeries',
           'NOAAIndicesTimeSeries', 'NOAAPredictIndicesTimeSeries', 'NoRHTimeSeries', 'RHESSISummaryTimeSeries']

source_names = {globals()[s]._source for s in __all__}
