# -*- coding: utf-8 -*-
"""Provides programs to process and analyze GOES data."""
from __future__ import absolute_import

import datetime

import matplotlib
from matplotlib import pyplot as plt

from sunpy.lightcurve import GenericLightCurve
from sunpy.time import TimeRange

__all__ = ['GOESLightCurve']

class GOESLightCurve(GenericLightCurve):
    """
    GOES LightCurve.

    Examples
    --------
    >>> from sunpy import lightcurve as lc

    >>> goes = lc.GOESLightCurve.create()
    >>> goes = lc.GOESLightCurve.create('2012/06/01', '2012/06/05')
    >>> goes.peek()

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
    def _get_goes_sat_num_and_data_type(self,start,end):
        """Parses the query time to determine which GOES satellite to use and the correct data type
        The operational times and dataformats are sourced from the GOES SEM API site:
        http://www.ngdc.noaa.gov/goes/sem/getData"""

        goes_operational={
        5:TimeRange('1986-01-01','1987-03-31'),
        6:TimeRange('1986-01-01','1994-12-31'),
        7:TimeRange('1987-03-01','1996-08-31'),
        8:TimeRange('1995-01-01','2003-06-16'),
        9:TimeRange('1996-04-01','1998-07-31'),
        10:TimeRange('1998-07-01','2009-12-31'),
        11:TimeRange('2000-07-01','2011-02-28'),
        12:TimeRange('2003-01-01','2010-08-31'),
        13:TimeRange('2006-08-01','2013-11-19'),
        14:TimeRange('2009-09-01','2012-11-24'),
        15:TimeRange('2010-03-26',datetime.datetime.utcnow())}

        goes_dformat={5:'xrs_1m', 6:'xrs_1m', 7:'xrs_1m', 8:'xrs_1m', 9:'xrs_1m', 10:'xrs_1m',
                      11:'xrs_1m', 12:'xrs_1m',13:'xrs_2s',14:'xrs_2s',15:'xrs_2s'}

        #find out which satellites were available. Start with newest first.
        for sat_num in range(15,5,-1):

            if ((start > goes_operational[sat_num].start() and start < goes_operational[sat_num].end()) and
                (end > goes_operational[sat_num].start() and end < goes_operational[sat_num].end())):
                #if true then the satellite with sat_num is available
                return sat_num,goes_dformat[sat_num]



        #if no satellites were found then raise an exception
        raise Exception, 'No operational GOES satellites found within specified time range'


    @classmethod
    def _get_url_from_timerange(cls, timerange, **kwargs):
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
        start = timerange.start()
        end = timerange.end()

        #find out which satellite and datatype to query from the query times


        sat_num, dtype = cls._get_goes_sat_num_and_data_type(start,end)
        # GOES query parameters
        params = {
            "satellite_number": sat_num,
            "data_type": dtype
        }
        params.update(kwargs)

        base_url = 'http://www.ngdc.noaa.gov/goes/sem/getData/goes%d/%s.csv'
        query_str = "?fromDate=%s&toDate=%s&file=true"

        url = (base_url + query_str) % (params['satellite_number'],
                                        params['data_type'],
                                        start.strftime("%Y%m%d"),
                                        end.strftime("%Y%m%d"))

        return [url]

    @classmethod
    def _is_datasource_for(cls, data, meta, source=None):
        if source is not None:
            return source.lower() == 'goes'
