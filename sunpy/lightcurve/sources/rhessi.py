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

__all__ = ['RHESSISummaryLightCurve']


class RHESSISummaryLightCurve(LightCurve):
    """
    RHESSI X-ray Summary LightCurve.

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

    def peek(self, title="RHESSI Observing Summary Count Rate", **kwargs):
        """Plots RHESSI Count Rate light curve"""
        figure = plt.figure()
        axes = plt.gca()

        dates = matplotlib.dates.date2num(self.data.index)

        # FIXME: not used?!
        lc_linecolors = ('black', 'pink', 'green', 'blue', 'brown', 'red',
                         'navy', 'orange', 'green')

        for item, frame in self.data.iteritems():
            axes.plot_date(dates, frame.values, '-', label=item, lw=2)

        axes.set_yscale("log")
        axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])

        axes.set_title('RHESSI Observing Summary Count Rates, Corrected')
        axes.set_ylabel('Corrected Count Rates s$^{-1}$ detector$^{-1}$')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')
        axes.legend()

        # @todo: display better tick labels for date range (e.g. 06/01 - 06/05)
        formatter = matplotlib.dates.DateFormatter('%H:%M')
        axes.xaxis.set_major_formatter(formatter)

        axes.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
        figure.autofmt_xdate()
        figure.show()


    @staticmethod
    def _parse_fits(filepath):
        """Parses a RHESSI FITS file"""
        header, d = rhessi.parse_obssumm_file(filepath)
        data = DataFrame(d['data'], columns=d['labels'], index=d['time'])

        return header, data
