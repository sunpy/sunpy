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
from sunpy.time import parse_time, TimeRange, is_time_in_given_format

__all__ = ['GOESLightCurve']


class GOESLightCurve(LightCurve):
    """
    GOES X-ray LightCurve. Provides GOES data back to 1981-01-01.
    Most recent data is usually available one or two days late.

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> from sunpy.time import TimeRange
    >>> goes = lc.GOESLightCurve.create(TimeRange('2012/06/01', '2012/06/05'))
    >>> goes.peek()

    References
    ----------
    | http://umbra.nascom.nasa.gov/goes/fits/
    """

    def peek(self, title="GOES Xray Flux"):
        """Plots GOES light curve is the usual manner"""
        figure = plt.figure()
        axes = plt.gca()

        dates = matplotlib.dates.date2num(parse_time(self.data.index))

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
        ax2.set_yticklabels((' ', 'A', 'B', 'C', 'M', 'X', ' '))

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

    @staticmethod
    def _parse_fits(filepath):
        """Parses a GOES FITS file from
        http://umbra.nascom.nasa.gov/goes/fits/"""
        fits = pyfits.open(filepath)
        header = fits[0].header
        if len(fits) == 4:
            if is_time_in_given_format(fits[0].header['DATE-OBS'], '%d/%m/%Y'):
                start_time = datetime.datetime.strptime(fits[0].header['DATE-OBS'], '%d/%m/%Y')
            elif is_time_in_given_format(fits[0].header['DATE-OBS'], '%d/%m/%y'):
                start_time = datetime.datetime.strptime(fits[0].header['DATE-OBS'], '%d/%m/%y')
            else:
                raise ValueError("Date not recognized")
            xrsb = fits[2].data['FLUX'][0][:, 0]
            xrsa = fits[2].data['FLUX'][0][:, 1]
            seconds_from_start = fits[2].data['TIME'][0]
        elif 1 <= len(fits) <= 3:
            start_time = parse_time(header['TIMEZERO'])
            seconds_from_start = fits[0].data[0]
            xrsb = fits[0].data[1]
            xrsa = fits[0].data[2]
        else:
            raise ValueError("Don't know how to parse this file")

        times = [start_time + datetime.timedelta(seconds=int(floor(s)),
                                                 microseconds=int((s - floor(s)) * 1e6)) for s in seconds_from_start]

        # remove bad values as defined in header comments
        xrsb[xrsb == -99999] = nan
        xrsa[xrsa == -99999] = nan

        # fix byte ordering
        newxrsa = xrsa.byteswap().newbyteorder()
        newxrsb = xrsb.byteswap().newbyteorder()

        data = DataFrame({'xrsa': newxrsa, 'xrsb': newxrsb}, index=times)

        return header, data
