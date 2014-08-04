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

    def peek(self, names=3, **kwargs):
        """Plots the LYRA data

        See: http://pandas.sourceforge.net/visualization.html
        """
        lyranames = (('Lyman alpha','Herzberg cont.','Al filter','Zr filter'),
                 ('120-123nm','190-222nm','17-80nm + <5nm','6-20nm + <2nm'))

        # Choose title if none was specified
        #if not kwargs.has_key("title"):
        #    if len(self.data.columns) > 1:
        #        kwargs['title'] = 'LYRA data'
        #    else:
        #        if self._filename is not None:
        #            base = self._filename
        #            kwargs['title'] = os.path.splitext(base)[0]
        #        else:
        #            kwargs['title'] = 'LYRA data'

        """Shows a plot of all four light curves"""
        figure = plt.figure()
        plt.subplots_adjust(left=0.17,top=0.94,right=0.94,bottom=0.15)
        axes = plt.gca()

        axes = self.data.plot(ax=axes, subplots=True, sharex=True, **kwargs)
        #plt.legend(loc='best')

        for i, name in enumerate(self.data.columns):
            if names < 3:
                name = lyranames[names][i]
            else:
                name = lyranames[0][i] + ' \n (' + lyranames[1][i] + ')'
            axes[i].set_ylabel( "%s %s" % (name, "\n (W/m**2)"),fontsize=9.5)

        axes[0].set_title("LYRA ("+ self.data.index[0].strftime('%Y-%m-%d') +")")
        axes[-1].set_xlabel("Time")
        for axe in axes:
            axe.locator_params(axis='y',nbins=6)

        figure.show()

        return figure


    @staticmethod
    def _get_url_for_date(date,**kwargs):
        """Returns a URL to the LYRA data for the specified date
        """
        dt = parse_time(date or datetime.datetime.utcnow())

        # Filename
        filename = "lyra_%s000000_lev%d_%s.fits" % (dt.strftime('%Y%m%d-'),
                                                    kwargs.get('level',2), 'std')
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
