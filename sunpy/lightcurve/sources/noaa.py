# -*- coding: utf-8 -*-
"""Provides programs to process and analyze NOAA Solar Cycle data."""
from __future__ import absolute_import

import datetime

from matplotlib import pyplot as plt
from pandas.io.parsers import read_csv
import numpy as np

from sunpy.lightcurve import LightCurve

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

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

    def plot(self, title=True, axes=None, plot_type=None, **plot_args):
        """Plots NOAA Indices as a function of time. An example is shown below.

        .. plot::

            from sunpy import lightcurve as lc
            from sunpy.data.sample import NOAAINDICES_LIGHTCURVE
            noaa = lc.NOAAIndicesLightCurve.create(NOAAINDICES_LIGHTCURVE)
            noaa.peek()

        Parameters
        ----------
        plot_type : str
            The type of plot required.

        axes : `matplotlib.axes.Axes`
            The axes on which to plot, defaults to the current axes.

        **plot_args : dict
            Any additional plot arguments that should be used when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """

        if axes is None:
            axes = plt.gca()

        if plot_type == None:
            plot_type = self._get_plot_types()[0]
        ylabel = ''

        if plot_type == self._get_plot_types()[0]:      # sunspot SWO
            self.data['sunspot SWO'].plot(axes=axes)
            self.data['sunspot SWO smooth'].plot(axes=axes)
            ylabel = 'Sunspot Number'
        elif plot_type == self._get_plot_types()[1]:    # sunspot RI
            self.data['sunspot RI'].plot(axes=axes)
            self.data['sunspot RI smooth'].plot(axes=axes)
            ylabel = 'Sunspot Number'
        elif plot_type == self._get_plot_types()[2]:    # sunspot compare
            self.data['sunspot SWO'].plot(axes=axes)
            self.data['sunspot RI'].plot(axes=axes)
            ylabel = 'Sunspot Number'
        elif plot_type == self._get_plot_types()[3]:    # radio
            self.data['radio flux smooth'].plot(axes=axes)
            self.data['radio flux'].plot(axes=axes)
            ylabel = 'Radio Flux [sfu]'
        elif plot_type == self._get_plot_types()[4]:    # ap index
            self.data['geomagnetic smooth'].plot(axes=axes)
            self.data['geomagnetic ap'].plot(axes=axes)
            ylabel = 'Geomagnetic AP Index'
        else:
            raise ValueError('Not a recognized plot type.')

        axes.set_ylim(0)
        axes.set_title(title)
        axes.set_ylabel(ylabel)

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.set_xlabel('Start time: ' + self.data['sunspot SWO'].index[0].strftime(TIME_FORMAT))

        #axes.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
        #axes.legend()
        plt.gcf().autofmt_xdate()

        return axes

    @classmethod
    def _get_plot_types(cls):
        return ['sunspot SWO', 'sunspot RI', 'sunspot compare', 'radio', 'ap index']

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
        comments = []
        with open(filepath, 'r') as fp:
            line = fp.readline()
            # Read header at top of file
            while line.startswith((":", "#")):
                comments += line
                line = fp.readline()
            fields = ('yyyy', 'mm', 'sunspot SWO', 'sunspot RI', 'sunspot ratio', 'sunspot SWO smooth', 'sunspot RI smooth', 'radio flux', 'radio flux smooth', 'geomagnetic ap', 'geomagnetic smooth')
            data = read_csv(fp, delim_whitespace=True, names = fields, comment='#', dtype={'yyyy':np.str, 'mm':np.str})
            data = data.dropna(how='any')
            timeindex = [datetime.datetime.strptime(x + y, '%Y%m') for x,y in zip(data['yyyy'], data['mm'])]
            data['time']=timeindex
            data = data.set_index('time')
            data = data.drop('mm',1)
            data = data.drop('yyyy',1)

            meta = {'comments': comments}
            units = [''] * len(data.columns)
            meta.update({'unit': units})
            meta.update({'instrume': 'NOAA'})
            meta.update({'obsrvtry': 'NOAA'})
            meta.update({'telescope': None})
            meta.update({'wavelnth': ''})
            meta.update({'waveunit': ''})

            return meta, data


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

    def plot(self, title=True, axes=None, plot_type=None, **plot_args):
        """Plots NOAA Indices as a function of time"""
        if axes is None:
            axes = plt.gca()

        if plot_type == None:
            plot_type = self._get_plot_types()[0]

        if plot_type == 'sunspot':
            ylabel = 'Sunspot Number'
            axes = self.data['sunspot'].plot(**plot_args)
            plt.fill_between(self.data.index, self.data['sunspot high'],
                         y2=self.data['sunspot low'], interpolate=True,
                         alpha=0.5)
        elif plot_type == 'radio':
            ylabel = 'Radio flux'
            axes = self.data['radio flux'].plot(**plot_args)
            plt.fill_between(self.data.index, self.data['radio flux high'],
                             y2=self.data['radio flux low'], interpolate=True,
                             alpha=0.5)
        else:
            raise ValueError('Not a recognized plot type.')

        axes = self.data['sunspot'].plot(color='b')
        self.data['sunspot low'].plot(linestyle='--', color='b')
        self.data['sunspot high'].plot(linestyle='--', color='b')

        axes.set_ylim(0)
        axes.set_title(title)
        axes.set_ylabel(ylabel)

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.set_xlabel('Start time: ' + self.data['sunspot'].index[0].strftime(TIME_FORMAT))

        plt.gcf().autofmt_xdate()
        return axes

    @classmethod
    def _get_plot_types(cls):
        return ['sunspot', 'radio']

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
        comments = ''
        with open(filepath, 'r') as fp:
            line = fp.readline()
            # Read header at top of file
            while line.startswith((":", "#")):
                comments += line
                line = fp.readline()
            fields = ('yyyy', 'mm', 'sunspot', 'sunspot low', 'sunspot high', 'radio flux', 'radio flux low', 'radio flux high')
            data = read_csv(filepath, delim_whitespace=True, names = fields, comment='#', skiprows=2, dtype={'yyyy':np.str, 'mm':np.str})
            data = data.dropna(how='any')
            timeindex = [datetime.datetime.strptime(x + y, '%Y%m') for x,y in zip(data['yyyy'], data['mm'])]
            data['time']=timeindex
            data = data.set_index('time')
            data = data.drop('mm',1)
            data = data.drop('yyyy',1)

            meta = {'comments': comments}
            units = [''] * len(data.columns)
            meta.update({'UNIT': units})
            meta.update({'INSTRUME': 'NOAA Prediction'})
            meta.update({'obsrvtry': 'NOAA'})
            meta.update({'telescope': None})
            meta.update({'wavelnth': ''})
            meta.update({'waveunit': ''})

            return meta, data
