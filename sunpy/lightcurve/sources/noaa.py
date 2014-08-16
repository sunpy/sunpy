# -*- coding: utf-8 -*-
"""Provides programs to process and analyze NOAA Solar Cycle data."""
from __future__ import absolute_import

import datetime

from matplotlib import pyplot as plt
from pandas.io.parsers import read_csv
import numpy as np

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

from sunpy.lightcurve import LightCurve

__all__ = ['NOAAIndicesLightCurve', 'NOAAPredictIndicesLightCurve']


class NOAAIndicesLightCurve(LightCurve):
    """NOAA Solar Cycle monthly indices.

    Solar activity is measured by a number of different values. The NOAA Solar
    Weather Prediction Center (SWPC) publishes the following indices. All of
    these indices are also provided as a 13-month running smoothed value. The 
    data columns in this object are

    * **sunspot SWO** - The SWO sunspot number is issued by the NOAA Space
    Weather Prediction Center (SWPC).
    * **sunspot SWO smooth** - Smoothed SWO sunspot number.
    * **sunspot RI** - The RI sunspot number is the official International
      Sunspot Number and is issued by the `Solar Influence Data Analysis Center (SDIC)
      <http://sidc.oma.be>`_ in Brussels, Belgium.
    * **sunspot SWO smooth** - Smoothed RI sunspot number.
    * **radio flux** - Radio flux at 10.7 cm is produced by
      `Penticon/Ottawa <http://www.ngdc.noaa.gov/stp/solar/flux.html>`_
      and the units are in sfu.
    * **radio flux smooth** - Smoothed radio flux at 10.7 cm.
    * **geomagnetic ap** - The Ap Geomagnetic Index is produced by the United
    States Air Force (USAF).
    * **geomagnetic smooth** - Smoothed geomagnetic index.
    * **sunspot ratio** - The ratio between the RI and SWO sunspot number.

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> noaa = lc.NOAAIndicesLightCurve.create()
    >>> noaa.peek()

    References
    ----------
    | http://www.swpc.noaa.gov/Data/index.html#indices
    | http://www.swpc.noaa.gov/ftpdir/weekly/README3
    | http://www.swpc.noaa.gov/ftpdir/weekly/RecentIndices.txt
    | http://www.swpc.noaa.gov/SolarCycle/
    """

    def plot(self, title='Solar Cycle Progression', type='sunspot SWO', axes=None, **plot_args):
        """Plots GOES light curve is the usual manner"""

        if axes is None:
            axes = plt.gca()

        if type == 'sunspot SWO':
            axes = self.data['sunspot SWO'].plot()
            self.data['sunspot SWO smooth'].plot()
            ylabel = 'Sunspot Number'
            xlabel = 'Start time: ' + self.data['sunspot SWO'].index[0].strftime(TIME_FORMAT)
        if type == 'sunspot RI':
            axes = self.data['sunspot RI'].plot()
            self.data['sunspot RI smooth'].plot()
            ylabel = 'Sunspot Number'
            xlabel = 'Start time: ' + self.data['sunspot RI'].index[0].strftime(TIME_FORMAT)
        if type == 'sunspot compare':
            axes = self.data['sunspot RI'].plot()
            self.data['sunspot SWO'].plot()
            ylabel = 'Sunspot Number'
            xlabel = 'Start time: ' + self.data['sunspot SWO'].index[0].strftime(TIME_FORMAT)
        if type == 'radio':
            axes = self.data['radio flux'].plot()
            self.data['radio flux smooth'].plot()
            ylabel = 'Radio Flux [sfu]'
        if type == 'geo':
            axes = self.data['geomagnetic ap'].plot()
            self.data['geomagnetic ap smooth'].plot()
            ylabel = 'Geomagnetic AP Index'

        axes.set_ylim(0)
        axes.set_title(title)
        axes.set_ylabel(ylabel)

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.set_xlabel(xlabel)

        axes.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.gcf().autofmt_xdate()
        
        return axes

    @classmethod
    def _get_default_uri(cls):
        """Return the url to download indices"""
        return "http://www.swpc.noaa.gov/ftpdir/weekly/RecentIndices.txt"

    @staticmethod
    def _get_url_for_date_range(*args, **kwargs):
        """Returns a URL for the specified date."""
        return NOAAIndicesLightCurve._get_default_uri()

    @staticmethod
    def _parse_csv(filepath):
        """Parses an NOAA indices csv"""
        header = []
        with open(filepath, 'r') as fp:
            line = fp.readline()
            # Read header at top of file
            while line.startswith((":", "#")):
                header += line
                line = fp.readline()
            fields = ('yyyy', 'mm', 'sunspot SWO', 'sunspot RI', 'sunspot ratio', 'sunspot SWO smooth', 'sunspot RI smooth', 'radio flux', 'radio flux smooth', 'geomagnetic ap', 'geomagnetic smooth')
            data = read_csv(fp, delim_whitespace=True, names = fields, comment='#', dtype={'yyyy':np.str, 'mm':np.str})
            data = data.dropna(how='any')
            timeindex = [datetime.datetime.strptime(x + '/' + y, '%Y/%m') for x,y in zip(data['yyyy'], data['mm'])]
            data['time']=timeindex
            data = data.set_index('time')
            data = data.drop('mm',1)
            data = data.drop('yyyy',1)
            return {'comments': header}, data


class NOAAPredictIndicesLightCurve(LightCurve):
    """NOAA Solar Cycle Predicted Progression

    The predictions are updated monthly and are produced by ISES. Observed
    values are initially the preliminary values which are replaced with the
    final values as they become available. The following data columns are provided.
    
    * **sunspot** - The predicted RI sunspot number is the official International Sunspot
      Number and is issued by the `Solar Influence Data Analysis Center (SDIC)
      <http://sidc.oma.be>`_ in Brussels, Belgium.
    * **sunspot low** - The low bound on the prediction.
    * **sunspot high** - The high bound on the prediction.
    * **radio flux** - The predicted radio flux at 10.7 cm is produced by
      `Penticon/Ottawa <http://www.ngdc.noaa.gov/stp/solar/flux.html>`_ and
      the units are in sfu.
    * **radio flux low** - The low bound on the radio predict.
    * **radio flux high** - The high bound on the radio predict.

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> noaa = lc.NOAAPredictIndicesLightCurve.create()
    >>> noaa.peek()

    References
    ----------
    | http://www.swpc.noaa.gov/Data/index.html#indices
    | http://www.swpc.noaa.gov/ftpdir/weekly/Predict.txt
    | http://www.swpc.noaa.gov/SolarCycle/
    """

    def plot(self, axes=None, title='Solar Cycle Prediction', **plot_args):
        """Plots NOAA Indices as a function of time"""
        if axes is None:
            axes = plt.gca()

        axes = self.data['sunspot'].plot(**plot_args)
        plt.fill_between(self.data.index, self.data['sunspot high'],
                         y2=self.data['sunspot low'], interpolate=True,
                         alpha=0.5)

        axes.set_ylim(0)
        axes.set_title(title)
        axes.set_ylabel('Sunspot Number')
        #axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')

        axes.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.gcf().autofmt_xdate()        

        return axes

    @classmethod
    def _get_default_uri(cls):
        """Return the url to download indices"""
        return "http://www.swpc.noaa.gov/ftpdir/weekly/Predict.txt"

    @staticmethod
    def _get_url_for_date_range(*args, **kwargs):
        """Returns a URL for the specified date."""
        return NOAAPredictIndicesLightCurve._get_default_uri()

    @staticmethod
    def _parse_csv(filepath):
        """Parses an NOAA indices csv"""
        header = ''
        with open(filepath, 'r') as fp:
            line = fp.readline()
            # Read header at top of file
            while line.startswith((":", "#")):
                header += line
                line = fp.readline()
            fields = ('yyyy', 'mm', 'sunspot', 'sunspot low', 'sunspot high', 'radio flux', 'radio flux low', 'radio flux high')
            data = read_csv(filepath, delim_whitespace=True, names = fields, comment='#', skiprows=2, dtype={'yyyy':np.str, 'mm':np.str})
            data = data.dropna(how='any')
            timeindex = [datetime.datetime.strptime(x + '/' + y, '%Y/%m') for x,y in zip(data['yyyy'], data['mm'])]
            data['time']=timeindex
            data = data.set_index('time')
            data = data.drop('mm',1)
            data = data.drop('yyyy',1)
            return {'comments': header}, data
