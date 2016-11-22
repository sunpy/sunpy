"""Provides programs to process and analyse NoRH lightcurve data."""

from __future__ import absolute_import

import datetime
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import pandas

from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time

from sunpy import config

from sunpy.extern.six.moves import urllib

TIME_FORMAT = config.get("general", "time_format")

__all__ = ['NoRHLightCurve']

class NoRHLightCurve(LightCurve):
    """
    Nobeyama Radioheliograph Correlation LightCurve.

    Nobeyama Radioheliograph (NoRH) is a radio telescope dedicated to observing
    the Sun. It consists of 84 parabolic antennas with 80 cm diameter,
    sitting on lines of 490 m long in the east/west and of 220 m long in the north/south.
    It observes the full solar disk at 17 GHz and 34 GHz with a temporal resolution
    down to 0.1 second resolution (typically 1 s). It is located in Japan at
    `35.941667, 138.475833 <https://www.google.com/maps/place/Nobeyama+radio+observatory/@35.9410098,138.470243,14z/data=!4m2!3m1!1s0x0:0xe5a3821a5f6a3c4b>`_.

    Its first observation was in April, 1992 and daily 8-hour observations are
    available starting June, 1992.

    Examples
    --------
    >>> import sunpy.lightcurve
    >>> norh = sunpy.lightcurve.NoRHLightCurve.create('~/Data/norh/tca110607')   # doctest: +SKIP
    >>> norh = sunpy.lightcurve.NoRHLightCurve.create('2011/08/10')
    >>> norh = sunpy.lightcurve.NoRHLightCurve.create('2011/08/10',wavelength='34')
    >>> norh.peek()   # doctest: +SKIP

    References
    ----------
    * `Nobeyama Radioheliograph Homepage <http://solar.nro.nao.ac.jp/norh/>`_
    * `Analysis Manual <http://solar.nro.nao.ac.jp/norh/doc/manuale/index.html>`_
    * `Nobeyama Correlation Plots <http://solar.nro.nao.ac.jp/norh/html/cor_plot/>`_
    """

    def peek(self, **kwargs):
        """Plots the NoRH lightcurve

        .. plot::

            from sunpy import lightcurve as lc
            from sunpy.data.sample import NORH_LIGHTCURVE
            norh = lc.NoRHLightCurve.create(NORH_LIGHTCURVE)
            norh.peek()

        Parameters
        ----------
        **kwargs : dict
            Any additional plot arguments that should be used
            when plotting.
        """

        plt.figure()
        axes = plt.gca()
        data_lab=self.meta['OBS-FREQ'][0:2] + ' ' + self.meta['OBS-FREQ'][2:5]
        axes.plot(self.data.index,self.data,label=data_lab)
        axes.set_yscale("log")
        axes.set_ylim(1e-4,1)
        axes.set_title('Nobeyama Radioheliograph')
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime(TIME_FORMAT))
        axes.set_ylabel('Correlation')
        axes.legend()
        plt.show()

    @classmethod
    def _get_url_for_date(cls, date, **kwargs):
        """
        This method retrieves the url for NoRH correlation data for the given
        date.
        """
        # default urllib password anonymous@ is not accepted by the NoRH FTP
        # server. include an accepted password in base url
        baseurl = 'ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/'
        # date is a datetime object
        if 'wavelength' in kwargs:
            if kwargs['wavelength'] == '34':
                final_url = urllib.parse.urljoin(
                    baseurl, date.strftime('%Y/%m/tcz%y%m%d'))
        else:
            final_url = urllib.parse.urljoin(
                baseurl, date.strftime('%Y/%m/tca%y%m%d'))

        return final_url

    @staticmethod
    def _parse_fits(filepath):
        """This method parses NoRH tca and tcz correlation files."""
        hdulist = fits.open(filepath)
        header = OrderedDict(hdulist[0].header)
        # For these NoRH files, the time series data is recorded in the primary
        # HDU
        data = hdulist[0].data

        # No explicit time array in FITS file, so construct the time array from
        # the FITS header
        obs_start_time=parse_time(header['DATE-OBS'] + 'T' + header['CRVAL1'])
        length = len(data)
        cadence = np.float(header['CDELT1'])
        sec_array = np.linspace(0, length-1, (length/cadence))

        norh_time = []
        for s in sec_array:
            norh_time.append(obs_start_time + datetime.timedelta(0,s))

        return header, pandas.DataFrame(data, index=norh_time)
