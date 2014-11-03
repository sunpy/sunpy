# -*- coding: utf-8 -*-
"""Provides programs to process and analyze PROBA2/LYRA data."""
from __future__ import absolute_import

import datetime
import urlparse
import sys

from matplotlib import pyplot as plt
from astropy.io import fits
import pandas

from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time
from sunpy.util.odict import OrderedDict

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

__all__ = ['LYRALightCurve']


class LYRALightCurve(LightCurve):
    """
    Proba-2 LYRA LightCurve.

    Examples
    --------
    >>> import sunpy

    >>> lyra = sunpy.lightcurve.LYRALightCurve.create()
    >>> lyra = sunpy.lightcurve.LYRALightCurve.create('~/Data/lyra/lyra_20110810-000000_lev2_std.fits')
    >>> lyra = sunpy.lightcurve.LYRALightCurve.create('2011/08/10')
    >>> lyra = sunpy.lightcurve.LYRALightCurve.create("http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev2_std.fits")
    >>> lyra.peek()

    References
    ----------
    | http://proba2.sidc.be/data/LYRA
    """

    def plot(self, title="LYRA", axes=None, type='channel 1', **plot_args):
        """Plots the LYRA data. Available plot types are
        
        """
        lyranames = (('Lyman alpha', 'Herzberg cont.', 'Al filter', 'Zr filter'),
                 ('120-123nm', '190-222nm', '17-80nm + <5nm', '6-20nm + <2nm'))

        """Shows a plot of all four light curves"""
        if axes is None:
            axes = plt.gca()

        if type == 'channel 1':
            axes = self.data['CHANNEL1'].plot(ax=axes, **plot_args)
            ylabel = lyranames[0][0] + ' \n (' + lyranames[1][0] + ')'
        if type == 'channel 2':
            axes = self.data['CHANNEL2'].plot(ax=axes, **plot_args)
            ylabel = lyranames[0][1] + ' \n (' + lyranames[1][1] + ')'
        if type == 'channel 3':
            axes = self.data['CHANNEL3'].plot(ax=axes, **plot_args)
            ylabel = lyranames[0][2] + ' \n (' + lyranames[1][2] + ')'
        if type == 'channel 4':
            axes = self.data['CHANNEL4'].plot(ax=axes, **plot_args)
            ylabel = lyranames[0][3] + ' \n (' + lyranames[1][3] + ')'
            
        ylabel+= "\n (W/m**2)"

        axes.set_xlabel('Start time: ' + self.data['CHANNEL1'].index[0].strftime(TIME_FORMAT))
        axes.set_ylabel(ylabel)
        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
        plt.gcf().autofmt_xdate()
        
        return axes

    @classmethod
    def _get_plot_types(cls):
        return ['channel 1', 'channel 2', 'channel 3', 'channel 4']
        
    @staticmethod
    def _get_url_for_date(date,**kwargs):
        """Returns a URL to the LYRA data for the specified date
        """
        dt = parse_time(date or datetime.datetime.utcnow())

        # Filename
        filename = "lyra_{0:%Y%m%d-}000000_lev{1:d}_std.fits".format(
            dt, kwargs.get('level',2))
        # URL
        base_url = "http://proba2.oma.be/lyra/data/bsd/"
        url_path = urlparse.urljoin(dt.strftime('%Y/%m/%d/'), filename)
        return urlparse.urljoin(base_url, url_path)

    @classmethod
    def _get_default_uri(cls):
        """Look for and download today's LYRA data"""
        return cls._get_url_for_date(datetime.datetime.utcnow())

    @staticmethod
    def _parse_fits(filepath):
        """Loads LYRA data from a FITS file"""
        # Open file with PyFITS
        hdulist = fits.open(filepath)
        fits_record = hdulist[1].data
        #secondary_header = hdulist[1].header

        # Start and end dates.  Different LYRA FITS files have
        # different tags for the date obs.
        if 'date-obs' in hdulist[0].header:
            start_str = hdulist[0].header['date-obs']
        elif 'date_obs' in hdulist[0].header:
            start_str = hdulist[0].header['date_obs']
        #end_str = hdulist[0].header['date-end']

        #start = datetime.datetime.strptime(start_str, '%Y-%m-%dT%H:%M:%S.%f')
        start = parse_time(start_str)
        #end = datetime.datetime.strptime(end_str, '%Y-%m-%dT%H:%M:%S.%f')

        # First column are times.  For level 2 data, the units are [s].
        # For level 3 data, the units are [min]
        if hdulist[1].header['TUNIT1'] == 's':
            times = [start + datetime.timedelta(seconds=int(n))
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
            #temporary patch for big-endian data bug on pandas 0.13
            if fits_record.field(i+1).dtype.byteorder == '>' and sys.byteorder =='little':
                table[col.name] = fits_record.field(i + 1).byteswap().newbyteorder()
            else:
                table[col.name] = fits_record.field(i + 1)

        # Return the header and the data
        return OrderedDict(hdulist[0].header), pandas.DataFrame(table, index=times)
