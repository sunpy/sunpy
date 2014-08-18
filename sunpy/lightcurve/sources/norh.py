"""Provides programs to process and analyse NoRH lightcurve data."""

from __future__ import absolute_import

import datetime
import urlparse
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import pandas

from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time
from sunpy.util.odict import OrderedDict

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

__all__ = ['NoRHLightCurve']

class NoRHLightCurve(LightCurve):
    """
    Nobeyama Radioheliograph LightCurve. The columns in this data are
    
    * **0** - The radio corrlation index.

    Examples
    --------
    >>> import sunpy

    >>> norh = sunpy.lightcurve.NoRHLightCurve.create('~/Data/norh/tca110607')
    >>> norh = sunpy.lightcurve.NoRHLightCurve.create('2011/08/10')
    >>> norh = sunpy.lightcurve.NoRHLightCurve.create('2011/08/10',wavelength='34')
    >>> norh.peek()

    References
    ----------
    | http://solar.nro.nao.ac.jp/norh/
    """

    def plot(self, title='Nobeyama Radioheliograph', axes=None, **plot_args):
        """Plots the NoRH lightcurve"""
 
        if axes is None:
            axes = plt.gca()

        data_label = self.meta['OBS-FREQ'][0:2] + ' ' + self.meta['OBS-FREQ'][2:5]
        axes.plot(self.data.index, self.data, label=data_label)
        axes.set_yscale("log")
        axes.set_ylim(1e-4,1)
        axes.set_title(title)
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime(TIME_FORMAT))
        axes.set_ylabel('Correlation')
        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.gcf().autofmt_xdate()
        
        return axes

    @classmethod
    def _get_plot_types(cls):
        return ['norh']

    @classmethod
    def _get_url_for_date(cls,date, **kwargs):
        """This method retrieves the url for NoRH correlation data for the given date."""
        #default urllib password anonymous@ is not accepted by the NoRH FTP server.
        #include an accepted password in base url
        baseurl='ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/'
        #date is a datetime object
        if 'wavelength' in kwargs:
            if kwargs['wavelength'] == '34':
                final_url=urlparse.urljoin(baseurl,date.strftime('%Y/%m/' + 'tcz' + '%y%m%d'))
        else:
            final_url=urlparse.urljoin(baseurl, date.strftime('%Y/%m/' + 'tca' + '%y%m%d'))

        return final_url

    @staticmethod
    def _parse_fits(filepath):
        """This method parses NoRH tca and tcz correlation files."""
        hdulist=fits.open(filepath)
        header=OrderedDict(hdulist[0].header)
        #for these NoRH files, the time series data is recorded in the primary HDU
        data=hdulist[0].data

        #No explicit time array in FITS file, so construct the time array from the FITS header
        obs_start_time=parse_time(header['DATE-OBS'] + 'T' + header['CRVAL1'])
        length=len(data)
        cadence=np.float(header['CDELT1'])
        sec_array=np.linspace(0, length-1, (length/cadence))

        norh_time=[]
        for s in sec_array:
            norh_time.append(obs_start_time + datetime.timedelta(0,s))

        return header, pandas.DataFrame(data, index=norh_time)
