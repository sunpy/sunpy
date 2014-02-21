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

__all__ = ['NOAAIndicesLightCurve', 'NOAAPredictIndicesLightCurve']

class NOAAIndicesLightCurve(LightCurve):
    """NOAA Solar cycle monthly indices monthly. 
        
    Solar activity is measured by a number of different values. The NOAA Solar Weather
    Prediction Center (SWPC) publishes the following indices. All of these indices are 
    also available as a 13-month running smoothed value.
      
    * The SWO sunspot number is issued by the NOAA Space Weather Prediction Center (SWPC)
        
    * The RI sunspot number is the official International Sunspot Number and is 
    issued by the `Solar Influence Data Analysis Center (SDIC) <http://sidc.oma.be>`_
    in Brussels, Belgium.
    
    * The ratio between the SWO and RI indices.
    
    * Radio flux at 10.7 cm is produced by `Penticon/Ottawa <http://www.ngdc.noaa.gov/stp/solar/flux.html>`_ 
    and the units are in sfu.
    
    * The Ap Geomagnetic Index is produced by the United States Air Force (USAF).
    
    * Predictions provided by ISES

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
    | http://www.swpc.noaa.gov/SolarCycle/
    """

    def plot(self, axes=None, **plot_args):
        """Plots NOAA Indices as a function of time"""
        if axes is None:
            axes = plt.gca()

        dates = matplotlib.dates.date2num(self.data.index)

        self.data['sunspot SWO'].plot()
        self.data['sunspot SWO smooth'].plot()

        axes.set_ylim(0)
        axes.set_title('Solar Cycle Sunspot Number Progression')
        axes.set_ylabel('Sunspot Number')
        #axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.legend()

        return axes

    @classmethod
    def _get_default_uri(cls):
        """Return the url to download indices"""
        return "http://www.swpc.noaa.gov/ftpdir/weekly/RecentIndices.txt"

    @staticmethod
    def _parse_csv(filepath):
        """Parses an NOAA indices csv"""
        with open(filepath, 'r') as fp:
            fields = ('yyyy', 'mm', 'sunspot SWO', 'sunspot RI', 'sunspot ratio', 'sunspot SWO smooth', 'sunspot RI smooth', 'radio flux', 'radio flux smoothed', 'geomagnetic ap', 'geomagnetic smooth')
            data = read_csv(fp, delim_whitespace=True, names = fields, comment='#', skiprows=2, dtype={'yyyy':np.str, 'mm':np.str})
            data = data.dropna(how='any')
            timeindex = [datetime.datetime.strptime(x + '/' + y, '%Y/%m') for x,y in zip(data['yyyy'], data['mm'])]
            data['time']=timeindex
            data = data.set_index('time')
            data = data.drop('mm',1)
            data = data.drop('yyyy',1)
            return "", data

class NOAAPredictIndicesLightCurve(LightCurve):
    """NOAA Solar cycle Predicted Progression

    The predictions are updated monthly and are produced by ISES. Observed values are 
    initially the preliminary values which are replaced with the final values as they 
    become available.
        
    The following predicted values are available.
              
    * The predicted RI sunspot number is the official International Sunspot Number and is 
    issued by the `Solar Influence Data Analysis Center (SDIC) <http://sidc.oma.be>`_
    in Brussels, Belgium.
    
    * The predicted radio flux at 10.7 cm is produced by `Penticon/Ottawa <http://www.ngdc.noaa.gov/stp/solar/flux.html>`_ 
    and the units are in sfu.
    
    Examples
    --------
    from sunpy import lightcurve as lc
    noaa = lc.NOAAPredictIndicesLightCurve.create()
    noaa.peek()

    References
    ----------
    | http://www.swpc.noaa.gov/Data/index.html#indices
    | http://www.swpc.noaa.gov/ftpdir/weekly/Predict.txt
    | http://www.swpc.noaa.gov/SolarCycle/
    """

    def plot(self, axes=None, **plot_args):
        """Plots NOAA Indices as a function of time"""
        if axes is None:
            axes = plt.gca()

        dates = matplotlib.dates.date2num(self.data.index)

        self.data['sunspot'].plot()
        self.data['sunspot low'].plot()
        self.data['sunspot high'].plot()

        axes.set_ylim(0)
        axes.set_title('Solar Cycle Sunspot Number Prediction')
        axes.set_ylabel('Sunspot Number')
        #axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.legend()

        return axes

    @classmethod
    def _get_default_uri(cls):
        """Return the url to download indices"""
        return "http://www.swpc.noaa.gov/ftpdir/weekly/Predict.txt"

    @staticmethod
    def _parse_csv(filepath):
        """Parses an NOAA indices csv"""
        with open(filepath, 'r') as fp:
            fields = ('yyyy', 'mm', 'sunspot', 'sunspot low', 'sunspot high', 'radio flux', 'radio flux low', 'radio flux high')
            data = read_csv(fp, delim_whitespace=True, names = fields, comment='#', skiprows=2, dtype={'yyyy':np.str, 'mm':np.str})
            data = data.dropna(how='any')
            timeindex = [datetime.datetime.strptime(x + '/' + y, '%Y/%m') for x,y in zip(data['yyyy'], data['mm'])]
            data['time']=timeindex
            data = data.set_index('time')
            data = data.drop('mm',1)
            data = data.drop('yyyy',1)
            return "", data

