# -*- coding: utf-8 -*-
"""Provides programs to process and analyze NOAA Solar Cycle data."""
from __future__ import absolute_import

import datetime

from matplotlib import pyplot as plt
from pandas.io.parsers import read_csv
import numpy as np

from sunpy.lightcurve import LightCurve

__all__ = ['NOAAIndicesLightCurve', 'NOAAPredictIndicesLightCurve']


class NOAAIndicesLightCurve(LightCurve):
    """NOAA Solar Cycle monthly indices.

    Solar activity is measured by a number of different values. The NOAA Solar
    Weather Prediction Center (SWPC) publishes the following indices. All of
    these indices are also provided as a 13-month running smoothed value.

    * The SWO sunspot number is issued by the NOAA Space Weather
      Prediction Center (SWPC)
    * The RI sunspot number is the official International Sunspot Number and is
      issued by the `Solar Influence Data Analysis Center (SDIC)
      <http://sidc.oma.be>`_ in Brussels, Belgium.
    * The ratio between the SWO and RI indices.
    * Radio flux at 10.7 cm is produced by
      `Penticon/Ottawa <http://www.ngdc.noaa.gov/stp/solar/flux.html>`_ and the
      units are in sfu.
    * The Ap Geomagnetic Index is produced by the United States Air Force (USAF).

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> noaa = lc.NOAAIndicesLightCurve.create()
    >>> noaa.peek()   # doctest: +SKIP

    References
    ----------
    * `Solar and Geomagnetic Indices Data Archive <http://legacy-www.swpc.noaa.gov/Data/index.html#indices>`_
    * `Recent solar indices <ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt>`_
    * `Indices Descriptions <ftp://ftp.swpc.noaa.gov/pub/weekly/README3>`_
    * `NOAA plots of Solar Cycle Progression <http://www.swpc.noaa.gov/products/solar-cycle-progression>`_
    * `NOAA Product List <http://www.swpc.noaa.gov/products-and-data>`_
    """

    def peek(self, type='sunspot SWO', **plot_args):
        """Plots NOAA Indices as a function of time. An example is shown below.

        .. plot::

            from sunpy import lightcurve as lc
            from sunpy.data.sample import NOAAINDICES_LIGHTCURVE
            noaa = lc.NOAAIndicesLightCurve.create(NOAAINDICES_LIGHTCURVE)
            noaa.peek()

        Parameters
        ----------
        type : str
            The type of plot required.

        **plot_args : dict
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """

        figure = plt.figure()
        axes = plt.gca()

        if type == 'sunspot SWO':
            axes = self.data['sunspot SWO'].plot()
            self.data['sunspot SWO smooth'].plot()
            axes.set_ylabel('Sunspot Number')
        if type == 'sunspot RI':
            axes = self.data['sunspot RI'].plot()
            self.data['sunspot RI smooth'].plot()
            axes.set_ylabel('Sunspot Number')
        if type == 'sunspot compare':
            axes = self.data['sunspot RI'].plot()
            self.data['sunspot SWO'].plot()
            axes.set_ylabel('Sunspot Number')
        if type == 'radio':
            axes = self.data['radio flux'].plot()
            self.data['radio flux smooth'].plot()
            axes.set_ylabel('Radio Flux [sfu]')
        if type == 'geo':
            axes = self.data['geomagnetic ap'].plot()
            self.data['geomagnetic ap smooth'].plot()
            axes.set_ylabel('Geomagnetic AP Index')

        axes.set_ylim(0)
        axes.set_title('Solar Cycle Progression')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.legend()

        figure.show()
        return figure

    @classmethod
    def _get_default_uri(cls):
        """Return the url to download indices"""
        return "ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt"

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
            timeindex = [datetime.datetime.strptime(x + y, '%Y%m') for x,y in zip(data['yyyy'], data['mm'])]
            data['time']=timeindex
            data = data.set_index('time')
            data = data.drop('mm',1)
            data = data.drop('yyyy',1)
            return {'comments': header}, data

class NOAAPredictIndicesLightCurve(LightCurve):
    """NOAA Solar Cycle Predicted Progression

    The predictions are updated monthly and are produced by ISES. Observed
    values are initially the preliminary values which are replaced with the
    final values as they become available.

    The following predicted values are available.

    * The predicted RI sunspot number is the official International Sunspot
      Number and is issued by the `Solar Influence Data Analysis Center (SDIC) <http://sidc.oma.be>`_ in Brussels, Belgium.
    * The predicted radio flux at 10.7 cm is produced by
      `Penticon/Ottawa <http://www.ngdc.noaa.gov/stp/solar/flux.html>`_ and the units are in sfu.

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> noaa = lc.NOAAPredictIndicesLightCurve.create()
    >>> noaa.peek()   # doctest: +SKIP

    References
    ----------
    * `Solar and Geomagnetic Indices Data Archive <http://legacy-www.swpc.noaa.gov/Data/index.html#indices>`_
    * `Predicted solar indices <http://services.swpc.noaa.gov/text/predicted-sunspot-radio-flux.txt>`_
    * `NOAA plots of Solar Cycle Progression <http://www.swpc.noaa.gov/products/solar-cycle-progression>`_
    * `NOAA Product List <http://www.swpc.noaa.gov/products-and-data>`_

    """

    def peek(self, **plot_args):
        """Plots predicted NOAA Indices as a function of time. An example is shown below.

        .. plot::

            from sunpy import lightcurve as lc
            from sunpy.data.sample import NOAAPREDICT_LIGHTCURVE
            noaa = lc.NOAAPredictIndicesLightCurve.create(NOAAPREDICT_LIGHTCURVE)
            noaa.peek()

        Parameters
        ----------
        **plot_args : dict
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """

        figure = plt.figure()
        axes = plt.gca()

        axes = self.data['sunspot'].plot(color='b')
        self.data['sunspot low'].plot(linestyle='--', color='b')
        self.data['sunspot high'].plot(linestyle='--', color='b')

        axes.set_ylim(0)
        axes.set_title('Solar Cycle Sunspot Number Prediction')
        axes.set_ylabel('Sunspot Number')
        # axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.legend()

        figure.show()
        return figure

    @classmethod
    def _get_default_uri(cls):
        """Return the url to download indices."""
        return "http://services.swpc.noaa.gov/text/predicted-sunspot-radio-flux.txt"

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
            timeindex = [datetime.datetime.strptime(x + y, '%Y%m') for x,y in zip(data['yyyy'], data['mm'])]
            data['time']=timeindex
            data = data.set_index('time')
            data = data.drop('mm',1)
            data = data.drop('yyyy',1)
            return {'comments': header}, data
