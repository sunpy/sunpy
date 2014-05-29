"""Provides programs to process and analyse NoRH lightcurve data."""

from __future__ import absolute_import

import datetime
import urlparse
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import pandas

from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

__all__ = ['NoRHLightCurve']

class NoRHLightCurve(LightCurve):
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

    def peek(self, **kwargs):
        """Plots the NoRH lightcurve"""
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

    @classmethod
    def _is_datasource_for(cls, data, meta, source=None):
        if meta is not None:
            return meta.pop('telescop', '').upper() == 'RADIOHELIOGRAPH'
        if source is not None:
            source = source.lower()
            return source == 'norh'
