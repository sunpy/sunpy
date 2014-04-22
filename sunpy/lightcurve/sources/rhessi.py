# -*- coding: utf-8 -*-
"""Provides programs to process and analyze GOES X-ray data."""
from __future__ import absolute_import

import datetime
import matplotlib.dates
from matplotlib import pyplot as plt  
from astropy.io import fits as pyfits
from numpy import nan
from numpy import floor
from pandas import DataFrame

from sunpy.lightcurve import LightCurve
from sunpy.time import TimeRange
from sunpy.instr import rhessi

__all__ = ['RHESSISummaryLightCurve']

class RHESSISummaryLightCurve(LightCurve):
    """
    RHESSI X-ray LightCurve.

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    
    >>> goes = lc.RHESSISummaryLightCurve.create()
    >>> goes = lc.RHESSISummaryLightCurve.create('2012/06/01', '2012/06/05')
    >>> goes.peek()

    References
    ----------
    | http://sprg.ssl.berkeley.edu/~jimm/hessi/hsi_obs_summ_soc.html#hsi_obs_summ_rate
    """

    def peek(self, title="GOES Xray Flux", **kwargs):
        """Plots GOES light curve is the usual manner"""
        figure = plt.figure()
        axes = plt.gca()

        dates = matplotlib.dates.date2num(self.data.index)

        axes.plot_date(dates, self.data['xrsa'], '-', 
                     label='0.5--4.0 $\AA$', color='blue', lw=2)
        axes.plot_date(dates, self.data['xrsb'], '-', 
                     label='1.0--8.0 $\AA$', color='red', lw=2)

        axes.set_yscale("log")
        axes.set_ylim(1e-9, 1e-2)
        axes.set_title(title)
        axes.set_ylabel('Watts m$^{-2}$')
        axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])

        ax2 = axes.twinx()
        ax2.set_yscale("log")
        ax2.set_ylim(1e-9, 1e-2)
        ax2.set_yticks((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
        ax2.set_yticklabels((' ','A','B','C','M','X',' '))

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')
        axes.legend()

        # @todo: display better tick labels for date range (e.g. 06/01 - 06/05)
        formatter = matplotlib.dates.DateFormatter('%H:%M')
        axes.xaxis.set_major_formatter(formatter)

        axes.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
        figure.autofmt_xdate()
        figure.show()

        return figure

    @classmethod
    def _get_default_uri(cls):
        """Retrieve GOES data from yesterday if no other data is specified"""
        today = datetime.datetime.today()
        yesterday = today - datetime.timedelta(days=1)
        return cls._get_url_for_date_range(yesterday, today)
       
    @staticmethod
    def _get_url_for_date_range(*args, **kwargs):
        """Returns a URL to the RHESSI data for the specified date range.

        Parameters
        ----------
        args : TimeRange, datetimes, date strings
            Date range should be specified using a TimeRange, or start
            and end dates at datetime instances or date strings.
        """
        # TimeRange
        if len(args) == 1 and isinstance(args[0], TimeRange):
            time_range = args[0]
            start = args[0].start()
            end = args[0].end()
        elif len(args) == 2:
            # Start & End date
            start = parse_time(args[0])
            end = parse_time(args[1])
            time_range = TimeRange(start, end)
            if end < start:
                print('Warning: start time (argument 1) > end time (argument 2)')
        url = rhessi.get_obssum_filename(time_range)
        print(url)
        return url
        
    @staticmethod
    def _parse_fits(filepath):
        """Parses a RHESSI FITS file"""
        header, d = rhessi.parse_obssumm_file(filepath)
        data = DataFrame(d['data'], columns=d['labels'], index=d['time'])
        
        return header, data
