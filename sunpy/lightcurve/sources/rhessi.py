# -*- coding: utf-8 -*-
"""Provides programs to process and analyze RHESSI X-ray data."""
from __future__ import absolute_import

import datetime
import matplotlib.dates
import matplotlib.pyplot as plt
from pandas import DataFrame

from sunpy.lightcurve import LightCurve
from sunpy.time import TimeRange, parse_time
from sunpy.instr import rhessi

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

__all__ = ['RHESSISummaryLightCurve']


class RHESSISummaryLightCurve(LightCurve):
    """
    RHESSI X-ray Summary LightCurve. The columns of data in this object are
    
    * **3 - 6 keV**
    * **6 - 12 keV**
    * **12 - 25 keV**
    * **25 - 50 keV**'
    * **50 - 100 keV**
    * **100 - 300 keV**
    * **300 - 800 keV**
    * **800 - 7000 keV**
    * **7000 - 20000 keV**

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> rhessi = lc.RHESSISummaryLightCurve.create()
    >>> rhessi = lc.RHESSISummaryLightCurve.create('2012/06/01', '2012/06/05')
    >>> rhessi.peek()

    References
    ----------
    | http://sprg.ssl.berkeley.edu/~jimm/hessi/hsi_obs_summ_soc.html#hsi_obs_summ_rate
    """

    def plot(self, title="RHESSI Observing Summary Count Rate", axes=None, type='rhessi', **plot_args):
        """Plots RHESSI Count Rate light curve"""

        if axes is None:
            axes = plt.gca()

        dates = matplotlib.dates.date2num(self.data.index)

        for item, frame in self.data.iteritems():
            axes.plot_date(dates, frame.values, '-', label=item, **plot_args)

        axes.set_yscale("log")
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime(TIME_FORMAT))

        axes.set_title(title)
        axes.set_ylabel('Count Rates s$^{-1}$ detector$^{-1}$')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')
        axes.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.gcf().autofmt_xdate()
        
        return axes

    @classmethod
    def _get_plot_types(cls):
        return ['rhessi']
      
    @classmethod
    def _get_default_uri(cls):
        """Retrieve the latest RHESSI data."""
        today = datetime.datetime.today()
        days_back = 3
        time_range = TimeRange(today - datetime.timedelta(days=days_back),
                               today - datetime.timedelta(days=days_back - 1))
        return cls._get_url_for_date_range(time_range)

    @staticmethod
    def _get_url_for_date_range(*args, **kwargs):
        """Returns a URL to the RHESSI data for the specified date range.

        Parameters
        ----------
        args : TimeRange, datetimes, date strings
            Date range should be specified using a TimeRange, or start
            and end dates at datetime instances or date strings.
        """
        if len(args) == 1 and isinstance(args[0], TimeRange):
            time_range = args[0]
        elif len(args) == 2:
            time_range = TimeRange(parse_time(args[0]), parse_time(args[1]))
            if time_range.end() < time_range.start():
                raise ValueError('start time > end time')
        url = rhessi.get_obssum_filename(time_range)
        return url

    @staticmethod
    def _parse_fits(filepath):
        """Parses a RHESSI FITS file"""
        header, d = rhessi.parse_obssumm_file(filepath)
        data = DataFrame(d['data'], columns=d['labels'], index=d['time'])

        return header, data