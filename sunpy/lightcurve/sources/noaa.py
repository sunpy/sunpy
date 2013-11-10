# -*- coding: utf-8 -*-
"""Provides programs to process and analyze GOES data."""
from __future__ import absolute_import

import datetime

import matplotlib
from matplotlib import pyplot as plt  
from pandas.io.parsers import read_csv
import numpy as np

import sunpy
from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time, TimeRange

__all__ = ['NOAAIndicesLightCurve']

class NOAAIndicesLightCurve(LightCurve):
    """NOAA Solar Indices weekly indices definition

    Examples
    --------
    from sunpy import lightcurve as lc
    noaa = lc.NOAAIndicesLightCurve.create()
    noaa.peek()

    References
    ----------
    | http://www.swpc.noaa.gov/Data/index.html#indices
    | http://www.swpc.noaa.gov/ftpdir/weekly/README3
    | http://www.swpc.noaa.gov/ftpdir/weekly/RecentIndices.txt
    """

    def plot(self, axes=None, **plot_args):
        """Plots NOAA Indices as a function of time"""
        if axes is None:
            axes = plt.gca()

        dates = matplotlib.dates.date2num(self.data.index)

        self.data['sunspot_swo'].plot()
        self.data['sunspot_smooth_swo'].plot()

        axes.set_ylim(0)
        axes.set_title('Solar Cycle Sunspot Number Progression')
        axes.set_ylabel('Sunspot Number')
        axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')
        axes.legend()

        # @todo: display better tick labels for date range (e.g. 06/01 - 06/05)
        formatter = matplotlib.dates.DateFormatter('%Y/%m')
        axes.xaxis.set_major_formatter(formatter)
        axes.fmt_xdata = matplotlib.dates.DateFormatter('%Y:%m')

        return axes

    @classmethod
    def _get_default_uri(cls):
        """Return the url to download indices"""
        return "http://www.swpc.noaa.gov/ftpdir/weekly/RecentIndices.txt"

    @staticmethod
    def _parse_csv(filepath):
        """Parses an NOAA indices csv"""
        with open(filepath, 'r') as fp:
            fields = ('yyyy', 'mm', 'sunspot_swo', 'sunspot_ri', 'sunspot_ratio', 'sunspot_smooth_swo', 'sunspot_smooth_ri', 'radioflux', 'radioflux_smooth', 'geo_ap', 'geo_ap_smooth')
            data = read_csv(fp, delim_whitespace=True, names = fields, comment='#', skiprows=2, dtype={'yyyy':np.str, 'mm':np.str})
            data = data.dropna(how='any')
            timeindex = [datetime.datetime.strptime(x + '/' + y, '%Y/%m') for x,y in zip(data['yyyy'], data['mm'])]
            data['time']=timeindex
            data = data.set_index('time')
            data = data.drop('mm',1)
            data = data.drop('yyyy',1)
            return "", data
