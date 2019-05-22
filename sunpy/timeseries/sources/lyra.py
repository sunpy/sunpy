"""
This module provies Proba-2 `~sunpy.timeseries.TimeSeries` source.
"""
import sys
from collections import OrderedDict

import pandas
import matplotlib as mpl
from matplotlib import pyplot as plt

import astropy.units as u
from astropy.time import TimeDelta

import sunpy.io
from sunpy import config
from sunpy.time import parse_time
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.metadata import MetaDict
from sunpy.visualization import peek_show

TIME_FORMAT = config.get("general", "time_format")

__all__ = ['LYRATimeSeries']


class LYRATimeSeries(GenericTimeSeries):
    """
    Proba-2 LYRA Lightcurve TimeSeries.

    LYRA (Large Yield RAdiometer) is an ultraviolet irradiance radiometer that observes the Sun in four passbands,
    chosen for their relevance to solar physics and space weather.
    LYRA is composed of three (redundant) units, each of them constituted of the same four channels:

    * 120-123 nm Lyman-alpha channel
    * 190-222 nm Herzberg continuum channel
    * Aluminium filter channel (17-80 nm + a contribution below 5 nm), including He II at 30.4 nm
    * Zirconium filter channel (6-20 nm + a contribution below 2 nm), rejecting He II

    LYRA can take data with cadences chosen in the 100Hz to 0.1Hz interval.

    PROBA2 was launched on 2 November 2009.

    Examples
    --------
    >>> import sunpy.timeseries
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> lyra = sunpy.timeseries.TimeSeries(sunpy.data.sample.LYRA_LEVEL3_TIMESERIES)  # doctest: +REMOTE_DATA
    >>> lyra.peek()   # doctest: +SKIP

    References
    ----------
    * `Proba2 SWAP Science Center <http://proba2.sidc.be/about/SWAP/>`_
    * `LYRA Data Homepage <http://proba2.sidc.be/data/LYRA>`_
    * `LYRA Instrument Homepage <http://proba2.sidc.be/about/LYRA>`_
    """
    # Class attribute used to specify the source class of the TimeSeries.
    _source = 'lyra'

    @peek_show
    def peek(self, names=3, **kwargs):
        """
        Plots the LYRA data. An example is shown below:

        .. plot::

            import sunpy.timeseries
            import sunpy.data.sample
            lyra = sunpy.timeseries.TimeSeries(sunpy.data.sample.LYRA_LEVEL3_TIMESERIES, source='LYRA')
            lyra.peek()

        Parameters
        ----------
        names : `int`, optional
            The number of columns to plot. Defaults to 3.
        **kwargs : `dict`
            Any additional plot arguments that should be used when plotting.
        """
        # Check we have a timeseries valid for plotting
        self._validate_data_for_ploting()

        lyranames = (('Lyman alpha', 'Herzberg cont.', 'Al filter', 'Zr filter'),
                     ('120-123nm', '190-222nm', '17-80nm + <5nm', '6-20nm + <2nm'))
        figure = plt.figure()
        plt.subplots_adjust(left=0.17, top=0.94, right=0.94, bottom=0.15)
        axes = plt.gca()

        axes = self.data.plot(ax=axes, subplots=True, sharex=True, **kwargs)

        for i, name in enumerate(self.data.columns):
            if names < 3:
                name = lyranames[names][i]
            else:
                name = lyranames[0][i] + ' \n (' + lyranames[1][i] + ')'
            axes[i].set_ylabel("{name} \n (W/m**2)".format(name=name), fontsize=9.5)

        axes[0].set_title("LYRA ({0:{1}})".format(self.data.index[0],TIME_FORMAT))
        axes[-1].set_xlabel("Time")
        for axe in axes:
            axe.locator_params(axis='y',nbins=6)

        return figure

    @classmethod
    def _parse_file(cls, filepath):
        """
        Parses Lyra FITS data files to create TimeSeries.

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
        Parses LYRA `astropy.io.fits.HDUList` from a FITS file.

        Parameters
        ----------
        hdulist : `astropy.io.fits.HDUList`
            A HDU list.
        """
        # Open file with PyFITS
        fits_record = hdulist[1].data

        metadata = MetaDict(OrderedDict(hdulist[0].header))
        start_str = metadata.get('date-obs', metadata.get('date_obs', ''))
        start = parse_time(start_str)

        # First column are times.  For level 2 data, the units are [s].
        # For level 3 data, the units are [min]
        if hdulist[1].header['TUNIT1'] == 's':
            times = start + TimeDelta(fits_record.field(0)*u.second)
        elif hdulist[1].header['TUNIT1'] == 'MIN':
            td = [int(n) for n in fits_record.field(0)]
            times = start + TimeDelta(td*u.minute)
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
        times.precision = 9
        data = pandas.DataFrame(table, index=times.isot.astype('datetime64'))
        data.sort_index(inplace=True)

        # Add the units data
        units = OrderedDict([('CHANNEL1', u.W/u.m**2),
                             ('CHANNEL2', u.W/u.m**2),
                             ('CHANNEL3', u.W/u.m**2),
                             ('CHANNEL4', u.W/u.m**2)])
        # TODO: check: http://www.wmo-sat.info/oscar/instruments/view/733
        return data, metadata, units

    @classmethod
    def is_datasource_for(cls, **kwargs):
        """
        Determines if the file corresponds to a LYRA LightCurve
        `~sunpy.timeseries.TimeSeries`.
        """
        # Check if source is explicitly assigned
        if 'source' in kwargs.keys():
            if kwargs.get('source', ''):
                return kwargs.get('source', '').lower().startswith(cls._source)
        # Check if HDU defines the source instrument
        if 'meta' in kwargs.keys():
            return kwargs['meta'].get('INSTRUME', '').startswith('LYRA')
