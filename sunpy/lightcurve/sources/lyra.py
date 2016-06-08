# -*- coding: utf-8 -*-
"""Provides programs to process and analyze PROBA2/LYRA data."""
from __future__ import absolute_import, division, print_function

import datetime
import sys
from collections import OrderedDict

from matplotlib import pyplot as plt
from astropy.io import fits
import pandas

from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time

from sunpy import config

from sunpy.extern.six.moves import urllib

TIME_FORMAT = config.get("general", "time_format")

__all__ = ['LYRALightCurve']

class LYRALightCurve(LightCurve):
    """
    Proba-2 LYRA LightCurve.

    LYRA (Large Yield RAdiometer) is an ultraviolet irradiance radiometer that
    observes the Sun in four passbands, chosen for their relevance to
    solar physics, aeronomy and space weather.
    LYRA is composed of three (redundant) units, each of them constituted of the
    same four channels:

        * 120-123 nm Lyman-alpha channel
        * 190-222 nm Herzberg continuum channel
        * Aluminium filter channel (17-80 nm + a contribution below 5 nm), including He II at 30.4 nm
        * Zirconium filter channel (6-20 nm + a contribution below 2 nm), rejecting He II

    LYRA can take data with cadences chosen in the 100Hz to 0.1Hz interval.

    PROBA2 was launched on 2 November 2009.

    This class can download and hold either Level 2 data (the default) which
    has sub-second resolution or Level 3 which is the Level 2 data averaged to
    one minute cadence. The level can be specified with the ``level`` keyword
    argument to `~sunpy.lightcurve.LyraLightCurve.create`.

    Examples
    --------
    >>> import sunpy
    >>> lyra = sunpy.lightcurve.LYRALightCurve.create()
    >>> lyra = sunpy.lightcurve.LYRALightCurve.create('~/Data/lyra/lyra_20110810-000000_lev2_std.fits')   # doctest: +SKIP
    >>> lyra = sunpy.lightcurve.LYRALightCurve.create('2011/08/10')
    >>> lyra = sunpy.lightcurve.LYRALightCurve.create('2011/08/10', level=3)
    >>> lyra = sunpy.lightcurve.LYRALightCurve.create("http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev2_std.fits")
    >>> lyra.peek()   # doctest: +SKIP

    References
    ----------
    * `Proba2 SWAP Science Center <http://proba2.sidc.be/about/SWAP/>`_
    * `LYRA Data Homepage <http://proba2.sidc.be/data/LYRA>`_
    * `LYRA Instrument Homepage <http://proba2.sidc.be/about/LYRA>`_
    """

    def plot(self, title=True, axes=None, plot_type=None, **plot_args):
        """Plots the LYRA data. Available plot types are

        .. plot::

            import sunpy.lightcurve
            from sunpy.data.sample import LYRA_LEVEL3_LIGHTCURVE
            lyra = sunpy.lightcurve.LYRALightCurve.create(LYRA_LEVEL3_LIGHTCURVE)
            lyra.peek()

        Parameters
        ----------
        names : int
            The number of columns to plot.

        **kwargs : dict
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """
        lyranames = (('Lyman alpha', 'Herzberg cont.', 'Al filter', 'Zr filter'),
                 ('120-123nm', '190-222nm', '17-80nm + <5nm', '6-20nm + <2nm'))

        if axes is None:
            axes = plt.gca()

        if plot_type == None:
            plot_type = self._get_plot_types()[0]

        if title is True:
            title = self.name

        ylabel = ''

        if plot_type == 'channel 1':
            axes = self.data['CHANNEL1'].plot(ax=axes, **plot_args)
            ylabel = lyranames[0][0] + ' \n (' + lyranames[1][0] + ") [" + str(self.unit) + "]"
        if plot_type == 'channel 2':
            axes = self.data['CHANNEL2'].plot(ax=axes, **plot_args)
            ylabel = lyranames[0][1] + ' \n (' + lyranames[1][1] + ") [" + str(self.unit) + "]"
        if plot_type == 'channel 3':
            axes = self.data['CHANNEL3'].plot(ax=axes, **plot_args)
            ylabel = lyranames[0][2] + ' \n (' + lyranames[1][2] + ") [" + str(self.unit) + "]"
        if plot_type == 'channel 4':
            axes = self.data['CHANNEL4'].plot(ax=axes, **plot_args)
            ylabel = lyranames[0][3] + ' \n (' + lyranames[1][3] + ") [" + str(self.unit) + "]"

        axes.set_xlabel('Start time: ' + self.data['CHANNEL1'].index[0].strftime(TIME_FORMAT))
        axes.set_ylabel(ylabel)
        axes.set_title(title)
        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        plt.gcf().autofmt_xdate()

        return axes

    @classmethod
    def _get_plot_types(cls):
        return ['channel 1', 'channel 2', 'channel 3', 'channel 4']

    @staticmethod
    def _get_url_for_date(date, **kwargs):
        """Returns a URL to the LYRA data for the specified date"""
        dt = parse_time(date or datetime.datetime.utcnow())

        # Filename
        filename = "lyra_{0:%Y%m%d-}000000_lev{1:d}_std.fits".format(
            dt, kwargs.get('level', 2))
        # URL
        base_url = "http://proba2.oma.be/lyra/data/bsd/"
        url_path = urllib.parse.urljoin(dt.strftime('%Y/%m/%d/'), filename)
        return urllib.parse.urljoin(base_url, url_path)

    @classmethod
    def _get_default_uri(cls):
        """Returns URL for latest LYRA data"""
        return cls._get_url_for_date(datetime.datetime.utcnow())

    @staticmethod
    def _parse_fits(filepath):
        """Parses LYRA data from a FITS file"""
        # Open file with PyFITS
        hdulist = fits.open(filepath)
        fits_record = hdulist[1].data
        # secondary_header = hdulist[1].header

        # Start and end dates.  Different LYRA FITS files have
        # different tags for the date obs.
        if 'date-obs' in hdulist[0].header:
            start_str = hdulist[0].header['date-obs']
        elif 'date_obs' in hdulist[0].header:
            start_str = hdulist[0].header['date_obs']
        # end_str = hdulist[0].header['date-end']

        # start = datetime.datetime.strptime(start_str, '%Y-%m-%dT%H:%M:%S.%f')
        start = parse_time(start_str)
        # end = datetime.datetime.strptime(end_str, '%Y-%m-%dT%H:%M:%S.%f')

        # First column are times.  For level 2 data, the units are [s].
        # For level 3 data, the units are [min]
        if hdulist[1].header['TUNIT1'] == 's':
            times = [start + datetime.timedelta(seconds=n)
                     for n in fits_record.field(0)]
        elif hdulist[1].header['TUNIT1'] == 'MIN':
            times = [start + datetime.timedelta(minutes=int(n))
                     for n in fits_record.field(0)]
        else:
            raise ValueError("Time unit in LYRA fits file not recognised.  "
                             "Value = {0}".format(hdulist[1].header['TUNIT1']))

        # Rest of columns are the data
        table = {}

        for i, col in enumerate(fits_record.columns[1:-1]):
            # temporary patch for big-endian data bug on pandas 0.13
            if fits_record.field(i+1).dtype.byteorder == '>' and sys.byteorder =='little':
                table[col.name] = fits_record.field(i + 1).byteswap().newbyteorder()
            else:
                table[col.name] = fits_record.field(i + 1)

        # Return the header and the data
        data = pandas.DataFrame(table, index=times)
        data.sort_index(inplace=True)
        header = OrderedDict(hdulist[0].header)
        header.update({'UNIT': "W m^-2"})
        header.update({'INSTRUME': 'LYRA'})
        header.update({'OBSRVTRY': 'PROBA2'})
        header.update({'TELESCOP': ['Lyman alpha', 'Herzberg cont.', 'Al filter', 'Zr filter']})
        header.update({'DETECTOR': header.get('TELESCOP')})
        header.update({'WAVELNTH': [[120, 123], [190, 222], [17, 80], [6, 20]]})
        header.update({'WAVEUNIT': 'nm'})
        return header, data
