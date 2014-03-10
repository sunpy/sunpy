"""Provides programs to process and analyse NoRH lightcurve data."""

from __future__ import absolute_import

import datetime
import urlparse
import numpy as np
import matplotlib.pyplot as plt

import pandas

from sunpy.lightcurve import GenericLightCurve
from sunpy.time import parse_time
from sunpy.map.header import MapMeta


__all__ = ['NoRHLightCurve']

class NoRHLightCurve(GenericLightCurve):
    """
    Nobeyama Radioheliograph LightCurve.

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
    def __init__(self, data, meta, **kwargs):
        obs_start_time=parse_time(meta['DATE-OBS'] + 'T' + meta['CRVAL1'])
        length=len(data)
        cadence=np.float(meta['CDELT1'])
        sec_array=np.linspace(0, length-1, (length/cadence))

        norh_time=[]
        for s in sec_array:
            norh_time.append(obs_start_time + datetime.timedelta(0,s))

        self.data = pandas.DataFrame(data, index=norh_time)
        self.meta = MapMeta(meta)

    def peek(self, **kwargs):
        """Plots the NoRH lightcurve"""
        plt.figure()
        axes = plt.gca()
        data_lab=self.meta['OBS-FREQ'][0:2] + ' ' + self.meta['OBS-FREQ'][2:5]
        axes.plot(self.data.index,self.data,label=data_lab)
        axes.set_yscale("log")
        axes.set_ylim(1e-4,1)
        axes.set_title('Nobeyama Radioheliograph')
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime('%Y-%m-%d %H:%M:%S UT'))
        axes.set_ylabel('Correlation')
        axes.legend()
        plt.show()

    @classmethod
    def _get_url_from_timerange(cls, timerange, **kwargs):
        days = timerange.get_days()
        urls = []
        for day in days:
            urls.append(cls._get_url_for_date(day, **kwargs))
        return urls

    @classmethod
    def _get_url_for_date(cls, date, **kwargs):
        """This method retrieves the url for NoRH correlation data for the given date."""

        # Hack to get around Python 2.x not backporting PEP 3102.
        wavelength = kwargs.pop('wavelength', None)

        #default urllib password anonymous@ is not accepted by the NoRH FTP server.
        #include an accepted password in base url
        baseurl='ftp://anonymous:mozilla@example.com@solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/'

        #date is a datetime.date object
        if wavelength == '34':
            final_url=urlparse.urljoin(baseurl,date.strftime('%Y/%m/' + 'tcz' + '%y%m%d'))
        else:
            final_url=urlparse.urljoin(baseurl, date.strftime('%Y/%m/' + 'tca' + '%y%m%d'))

        return final_url

    @classmethod
    def _is_datasource_for(cls, data, meta, source=None):
        if meta is not None:
            return meta.pop('telescop', '').upper() == 'RADIOHELIOGRAPH'
        if source is not None:
            source = source.lower()
            return source == 'norh'