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
from sunpy.util import net

__all__ = ['GOESLightCurve']


class GOESLightCurve(LightCurve):
    """
    GOES XRS LightCurve

    Each GOES satellite there are two X-ray Sensors (XRS) which provide solar X
    ray fluxes for the wavelength bands of 0.5 to 4 Å (short channel)
    and 1 to 8 Å (long channel). Most recent data is usually available one or two days late.

    Data is available starting on 1981/01/01.

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> from sunpy.time import TimeRange
    >>> goes = lc.GOESLightCurve.create(TimeRange('2012/06/01', '2012/06/05'))
    >>> goes.peek()   # doctest: +SKIP

    References
    ----------
    * `GOES Mission Homepage <http://www.goes.noaa.gov>`_
    * `GOES XRS Homepage <http://www.swpc.noaa.gov/products/goes-x-ray-flux>`_
    * `GOES XRS Guide <http://ngdc.noaa.gov/stp/satellite/goes/doc/GOES_XRS_readme.pdf>`_
    * `NASCOM Data Archive <http://umbra.nascom.nasa.gov/goes/fits/>`_
    """

    def peek(self, title="GOES Xray Flux"):
        """Plots GOES XRS light curve is the usual manner. An example is shown
        below.

        .. plot::

            from sunpy import lightcurve as lc
            from sunpy.data.sample import GOES_LIGHTCURVE
            goes = lc.GOESLightCurve.create(GOES_LIGHTCURVE)
            goes.peek()

        Parameters
        ----------
        title : str
            The title of the plot.

        **kwargs : dict
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """
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

    @classmethod
    def _get_default_uri(cls):
        """Returns the URL for the latest GOES data."""
        now = datetime.datetime.utcnow()
        time_range = TimeRange(datetime.datetime(now.year, now.month, now.day), now)
        url_does_exist = net.url_exists(cls._get_url_for_date_range(time_range))
        while not url_does_exist:
            time_range = TimeRange(time_range.start-datetime.timedelta(days=1),
                                   time_range.start)
            url_does_exist = net.url_exists(cls._get_url_for_date_range(time_range))
        return cls._get_url_for_date_range(time_range)

    @classmethod
    def _get_goes_sat_num(self, start, end):
        """Parses the query time to determine which GOES satellite to use."""

        goes_operational = {
        2: TimeRange('1980-01-04', '1983-05-01'),
        5: TimeRange('1983-05-02', '1984-08-01'),
        6: TimeRange('1983-06-01', '1994-08-19'),
        7: TimeRange('1994-01-01', '1996-08-14'),
        8: TimeRange('1996-03-21', '2003-06-19'),
        9: TimeRange('1997-01-01', '1998-09-09'),
        10: TimeRange('1998-07-10', '2009-12-02'),
        11: TimeRange('2006-06-20', '2008-02-16'),
        12: TimeRange('2002-12-13', '2007-05-09'),
        13: TimeRange('2006-08-01', '2006-08-01'),
        14: TimeRange('2009-12-02', '2010-11-05'),
        15: TimeRange('2010-09-01', datetime.datetime.utcnow())}

        sat_list = []
        for sat_num in goes_operational:
            if ((start >= goes_operational[sat_num].start and
                 start <= goes_operational[sat_num].end and
                (end >= goes_operational[sat_num].start and
                 end <= goes_operational[sat_num].end))):
                # if true then the satellite with sat_num is available
                sat_list.append(sat_num)

        if not sat_list:
            # if no satellites were found then raise an exception
            raise Exception('No operational GOES satellites within time range')
        else:
            return sat_list

    @staticmethod
    def _get_url_for_date_range(*args):
        """Returns a URL to the GOES data for the specified date.

        Parameters
        ----------
        args : `~sunpy.time.TimeRange`, `datetime.datetime`, str
            Date range should be specified using a TimeRange, or start
            and end dates at datetime instances or date strings.
        satellite_number : int
            GOES satellite number (default = 15)
        data_type : str
            Data type to return for the particular GOES satellite. Supported
            types depend on the satellite number specified. (default = xrs_2s)
        """
        # TimeRange
        if len(args) == 1 and isinstance(args[0], TimeRange):
            start = args[0].start
            end = args[0].end
        elif len(args) == 2:
            start = parse_time(args[0])
            end = parse_time(args[1])
        if end < start:
            raise ValueError('start time > end time')

        # find out which satellite and datatype to query from the query times
        sat_num = GOESLightCurve._get_goes_sat_num(start, end)
        base_url = 'http://umbra.nascom.nasa.gov/goes/fits/'

        if start < parse_time('1999/01/15'):
            url = base_url + "{date:%Y}/go{sat:02d}{date:%y%m%d}.fits".format(
                date=start, sat=sat_num[0])
        else:
            url = base_url + "{date:%Y}/go{sat:02d}{date:%Y%m%d}.fits".format(
                date=start, sat=sat_num[0])
        return url

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
        data.sort_index(inplace=True)
        return header, data
