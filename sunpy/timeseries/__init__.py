"""
SunPy's TimeSeries module provides a datatype for 1D time series data, replacing the SunPy LightCurve module.

Currently the objects can be instansiated from files (such as CSV and FITS) and urls to these files, but don't include data downloaders for their specific instruments as this will become part of the universal downloader.
"""
from __future__ import absolute_import
from sunpy.timeseries.metadata import TimeSeriesMetaData
from sunpy.timeseries.timeseries_factory import TimeSeries
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.timeseries.sources.eve import *
from sunpy.timeseries.sources.goes import *
from sunpy.timeseries.sources.noaa import *
from sunpy.timeseries.sources.lyra import *
from sunpy.timeseries.sources.norh import *
from sunpy.timeseries.sources.rhessi import *
from sunpy.timeseries.sources.fermi_gbm import *