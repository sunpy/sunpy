# -*- coding: utf-8 -*-
"""GOES XRS TimeSeries subclass definitions."""
from __future__ import absolute_import, print_function, division
# pylint: disable=W0221,W0222,E1101,E1121

from collections import OrderedDict
import datetime
import matplotlib.dates
from matplotlib import pyplot as plt
import numpy as np
from pandas import DataFrame

import sunpy.io
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.time import parse_time, TimeRange, is_time_in_given_format
from sunpy.util.metadata import MetaDict

from astropy import units as u

__author__ = ["Alex Hamilton"]
__email__ = "####"

__all__ = ['XRSTimeSeries']


class XRSTimeSeries(GenericTimeSeries):
    """
    GOES XRS Time Series

    Each GOES satellite there are two X-ray Sensors (XRS) which provide solar X
    ray fluxes for the wavelength bands of 0.5 to 4 Å (short channel)
    and 1 to 8 Å (long channel). Most recent data is usually available one or two days late.

    Data is available starting on 1981/01/01.

    Examples
    --------
    >>> import sunpy.timeseries
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)  # doctest: +REMOTE_DATA
    >>> goes.peek()   # doctest: +SKIP

    References
    ----------
    * `GOES Mission Homepage <http://www.goes.noaa.gov>`_
    * `GOES XRS Homepage <https://www.swpc.noaa.gov/products/goes-x-ray-flux>`_
    * `GOES XRS Guide <https://ngdc.noaa.gov/stp/satellite/goes/doc/GOES_XRS_readme.pdf>`_
    * `NASCOM Data Archive <https://umbra.nascom.nasa.gov/goes/fits/>`_

    Notes:
    https://umbra.nascom.nasa.gov/goes/fits/goes_fits_files_notes.txt
    """

    # Class attribute used to specify the source class of the TimeSeries.
    _source = 'xrs'

    def peek(self, title="GOES Xray Flux"):
        """Plots GOES XRS light curve is the usual manner. An example is shown
        below.

        .. plot::

            import sunpy.timeseries
            import sunpy.data.sample
            ts_goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
            ts_goes.peek()

        Parameters
        ----------
        title : `str`
            The title of the plot.

        **kwargs : `dict`
            Any additional plot arguments that should be used when plotting.
        """
        # Check we have a timeseries valid for plotting
        self._validate_data_for_ploting()

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

    # ToDo: is this part of the DL pipeline? If so delete.
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

    @classmethod
    def _parse_file(cls, filepath):
        """Parses a GOES/XRS FITS file from
        https://umbra.nascom.nasa.gov/goes/fits/"""

        hdus = sunpy.io.read_file(filepath)
        return cls._parse_hdus(hdus)

    @classmethod
    def _parse_hdus(cls, hdulist):
        header = MetaDict(OrderedDict(hdulist[0].header))
        if len(hdulist) == 4:
            if is_time_in_given_format(hdulist[0].header['DATE-OBS'], '%d/%m/%Y'):
                start_time = datetime.datetime.strptime(hdulist[0].header['DATE-OBS'], '%d/%m/%Y')
            elif is_time_in_given_format(hdulist[0].header['DATE-OBS'], '%d/%m/%y'):
                start_time = datetime.datetime.strptime(hdulist[0].header['DATE-OBS'], '%d/%m/%y')
            else:
                raise ValueError("Date not recognized")
            xrsb = hdulist[2].data['FLUX'][0][:, 0]
            xrsa = hdulist[2].data['FLUX'][0][:, 1]
            seconds_from_start = hdulist[2].data['TIME'][0]
        elif 1 <= len(hdulist) <= 3:
            start_time = parse_time(header['TIMEZERO'])
            seconds_from_start = hdulist[0].data[0]
            xrsb = hdulist[0].data[1]
            xrsa = hdulist[0].data[2]
        else:
            raise ValueError("Don't know how to parse this file")

        times = [start_time + datetime.timedelta(seconds=int(np.floor(s)),
                                                    microseconds=int((s - np.floor(s)) * 1e6)) for s in seconds_from_start]

        # remove bad values as defined in header comments
        xrsb[xrsb == -99999] = np.nan
        xrsa[xrsa == -99999] = np.nan

        # fix byte ordering
        newxrsa = xrsa.byteswap().newbyteorder()
        newxrsb = xrsb.byteswap().newbyteorder()

        data = DataFrame({'xrsa': newxrsa, 'xrsb': newxrsb}, index=times)
        data.sort_index(inplace=True)

        # Add the units
        units = OrderedDict([('xrsa', u.W/u.m**2),
                             ('xrsb', u.W/u.m**2)])
        return data, header, units

    @classmethod
    def is_datasource_for(cls, **kwargs):
        """Determines if header corresponds to a GOES lightcurve TimeSeries"""
        if 'source' in kwargs.keys():
            if kwargs.get('source', ''):
                return kwargs.get('source', '').lower().startswith(cls._source)
        if 'meta' in kwargs.keys():
            return kwargs['meta'].get('TELESCOP', '').startswith('GOES')
