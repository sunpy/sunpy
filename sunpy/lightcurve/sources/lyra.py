# -*- coding: utf-8 -*-
"""Provides programs to process and analyze PROBA2/LYRA data."""
from __future__ import absolute_import

import datetime
import urlparse

from matplotlib import pyplot as plt
from astropy.io import fits
import pandas

from sunpy.lightcurve import LightCurve 
from sunpy.time import parse_time

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
        axes = plt.gca()

        axes = self.data.plot(ax=axes, subplots=True, sharex=True, **kwargs)
        #plt.legend(loc='best')

        for i, name in enumerate(self.data.columns):
            if names < 3:
                name = lyranames[names][i]
            else:
                name = lyranames[0][i] + ' (' + lyranames[1][i] + ')'
            axes[i].set_ylabel("%s (%s)" % (name, "W/m**2"))

        axes[0].set_title("LYRA ("+ self.data.index[0].strftime('%Y-%m-%d') +")")
        axes[-1].set_xlabel("Time")

        figure.show()

        return figure


    @staticmethod
    def _get_url_for_date(date):
        """Returns a URL to the LYRA data for the specified date
        """
        dt = parse_time(date or datetime.datetime.utcnow())

        # Filename
        filename = "lyra_%s000000_lev%d_%s.fits" % (dt.strftime('%Y%m%d-'),
                                                    2, 'std')
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

        # First column are times
        times = [start + datetime.timedelta(0, n) for n in fits_record.field(0)]

        # Rest of columns are the data
        table = {}

        for i, col in enumerate(fits_record.columns[1:-1]):
            table[col.name] = fits_record.field(i + 1)

        # Return the header and the data
        return hdulist[0].header, pandas.DataFrame(table, index=times)


