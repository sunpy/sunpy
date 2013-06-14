# -*- coding: utf-8 -*-
"""Provides programs to process and analyze GOES data."""
from __future__ import absolute_import

import datetime

import matplotlib
from matplotlib import pyplot as plt  
from pandas.io.parsers import read_csv

import sunpy
from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time, TimeRange

__all__ = ['GOESLightCurve']

class GOESLightCurve(LightCurve):
    """GOES light curve definition

    Examples
    --------
    import sunpy
    goes = sunpy.lightcurve.GOESLightCurve.create()
    goes = sunpy.lightcurve.GOESLightCurve.create('2012/06/01', '2012/06/05')
    goes.peek()

    References
    ----------
    | http://www.ngdc.noaa.gov/goes/sem
    | http://www.ngdc.noaa.gov/goes/sem/getData/goes15
    """

    def peek(self, title="GOES Xray Flux", **kwargs):
        """Plots GOES light curve is the usual manner"""
        figure = plt.figure()
        axes = plt.gca()

        dates = matplotlib.dates.date2num(self.data.index)

        axes.plot_date(dates, self.data['A_FLUX'], '-', 
                     label='0.5--4.0 $\AA$', color='blue', lw=2)
        axes.plot_date(dates, self.data['B_FLUX'], '-', 
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
        """Retrieve XRS 2s data from yesterday (most recent data available using
        SEM API) if no other data is specified"""
        today = datetime.datetime.today()
        yesterday = today - datetime.timedelta(days=1)
        return cls._get_url_for_date_range(yesterday, today)

    @staticmethod
    def _get_url_for_date_range(*args, **kwargs):
        """Returns a URL to the GOES data for the specified date.

        Parameters
        ----------
        args : TimeRange, datetimes, date strings
            Date range should be specified using a TimeRange, or start
            and end dates at datetime instances or date strings.
        satellite_number : int
            GOES satellite number (default = 15)
        data_type : string
            Data type to return for the particular GOES satellite. Supported
            types depend on the satellite number specified. (default = xrs_2s) 
        """
        # TimeRange
        if len(args) == 1 and isinstance(args[0], TimeRange):
            start = args[0].start()
            end = args[0].end()
        elif len(args) == 2:
            # Start & End date
            start = parse_time(args[0])
            end = parse_time(args[1])
            if end < start:
                print('Warning: start time (argument 1) > end time (argument 2)')

        # GOES query parameters
        params = {
            "satellite_number": 15,
            "data_type": "xrs_2s"
        }
        params.update(kwargs)

        base_url = 'http://www.ngdc.noaa.gov/goes/sem/getData/goes%d/%s.csv'
        query_str = "?fromDate=%s&toDate=%s&file=true" 

        url = (base_url + query_str) % (params['satellite_number'], 
                                        params['data_type'],
                                        start.strftime("%Y%m%d"), 
                                        end.strftime("%Y%m%d"))

        return url

    @staticmethod
    def _parse_csv(filepath):
        """Parses an GOES CSV"""
        with open(filepath, 'rb') as fp:
            # @todo: check for:
            # "No-Data-Found for the time period requested..." error
            return "", read_csv(fp, sep=",", index_col=0, parse_dates=True)
