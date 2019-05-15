"""
This module provies a Nobeyama Radioheliograph `~sunpy.timeseries.TimeSeries`
source.
"""
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas

import astropy.units as u
from astropy.time import TimeDelta

import sunpy.io
from sunpy import config
from sunpy.time import parse_time
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.metadata import MetaDict
from sunpy.visualization import peek_show

TIME_FORMAT = config.get("general", "time_format")

__all__ = ['NoRHTimeSeries']


class NoRHTimeSeries(GenericTimeSeries):
    """
    Nobeyama Radioheliograph Correlation Lightcurve TimeSeries.

    Nobeyama Radioheliograph (NoRH) is a radio telescope dedicated to observing the Sun.
    It consists of 84 parabolic antennas with a 80 cm diameter,
    sitting on lines of 490 m long in the east/west and of 220 m long in the north/south.
    It observes the full solar disk at 17 GHz and 34 GHz with a temporal resolution
    down to 0.1 second resolution (typically 1 second).

    Its first observation was in April, 1992 and daily 8-hour observations are available starting June, 1992.

    Examples
    --------
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> import sunpy.timeseries
    >>> norh = sunpy.timeseries.TimeSeries(sunpy.data.sample.NORH_TIMESERIES, source='NoRH')  # doctest: +REMOTE_DATA
    >>> norh.peek()   # doctest: +SKIP

    References
    ----------
    * `Nobeyama Radioheliograph Homepage <https://solar.nro.nao.ac.jp/norh/>`_
    * `Analysis Manual <https://solar.nro.nao.ac.jp/norh/doc/manuale/index.html>`_
    * `Nobeyama Correlation Plots <https://solar.nro.nao.ac.jp/norh/html/cor_plot/>`_
    """
    # Class attribute used to specify the source class of the TimeSeries.
    _source = 'norh'

    def __init__(self, data, header, units, **kwargs):
        super(NoRHTimeSeries, self).__init__(data, header, units, **kwargs)

    @peek_show
    def peek(self):
        """
        Plot the NoRH lightcurve TimeSeries.

        .. plot::

            import sunpy.data.sample
            import sunpy.timeseries
            norh = sunpy.timeseries.TimeSeries(sunpy.data.sample.NORH_TIMESERIES, source='NoRH')
            norh.peek()
        """
        # Check we have a timeseries valid for plotting
        self._validate_data_for_ploting()

        figure = plt.figure()
        axes = plt.gca()
        data_lab = str(self.meta.get('OBS-FREQ').values()).replace('[', '').replace(
            ']', '').replace('\'', '')
        axes.plot(self.data.index, self.data, label=data_lab)
        axes.set_yscale("log")
        axes.set_ylim(1e-4, 1)
        axes.set_title('Nobeyama Radioheliograph')
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime(TIME_FORMAT))
        axes.set_ylabel('Correlation')
        axes.legend()

        return figure


    @classmethod
    def _parse_file(cls, filepath):
        """
        This method parses NoRH FITS files.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        hdus = sunpy.io.read_file(filepath)
        return cls._parse_hdus(hdus)

    @classmethod
    def _parse_hdus(cls, hdulist):
        """
        This method parses a NoRH `astropy.io.fits.HDUList`.

        Parameters
        ----------
        hdulist : `astropy.io.fits.HDUList`
            A HDU list.
        """
        header = MetaDict(OrderedDict(hdulist[0].header))
        # For these NoRH files, the time series data is recorded in the primary
        # HDU
        data = hdulist[0].data

        # No explicit time array in FITS file, so construct the time array from
        # the FITS header
        obs_start_time = parse_time(header['DATE-OBS'] + 'T' + header['CRVAL1'])
        length = len(data)
        cadence = np.float(header['CDELT1'])
        sec_array = np.linspace(0, length - 1, int(length / cadence))

        norh_time = obs_start_time + TimeDelta(sec_array*u.second)
        norh_time.precision = 9
        norh_time = norh_time.isot.astype('datetime64')

        # Add the units data
        units = OrderedDict([('Correlation Coefficient', u.dimensionless_unscaled)])
        # Todo: check units used.
        return pandas.DataFrame(
            data, index=norh_time, columns=('Correlation Coefficient', )), header, units

    @classmethod
    def is_datasource_for(cls, **kwargs):
        """
        Determines if header corresponds to a Nobeyama Radioheliograph
        Correlation `~sunpy.timeseries.TimeSeries`.
        """
        if 'source' in kwargs.keys():
            if kwargs.get('source', ''):
                return kwargs.get('source', '').lower().startswith(cls._source)
        if 'meta' in kwargs.keys():
            return kwargs['meta'].get('ORIGIN', '').startswith('NOBEYAMA RADIO OBS')
