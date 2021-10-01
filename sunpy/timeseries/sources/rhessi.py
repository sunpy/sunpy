"""
This module provies a RHESSI `~sunpy.timeseries.TimeSeries` source.
"""
import datetime
import itertools
from collections import OrderedDict

import matplotlib.dates
import matplotlib.pyplot as plt
import numpy as np
from pandas import DataFrame

import astropy.units as u
from astropy.time import TimeDelta

import sunpy.io
from sunpy.time import parse_time
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.metadata import MetaDict
from sunpy.visualization import peek_show

__all__ = ['RHESSISummaryTimeSeries']


def uncompress_countrate(compressed_countrate):
    """
    Convert the compressed count rate inside of observing summary file from a
    compressed byte to a true count rate.

    Parameters
    ----------
    compressed_countrate : `byte` array
        A compressed count rate returned from an observing summary file.

    References
    ----------
    `Hsi_obs_summ_decompress.pro <https://hesperia.gsfc.nasa.gov/ssw/hessi/idl/qlook_archive/hsi_obs_summ_decompress.pro>`_
    """

    # Ensure uncompressed counts are between 0 and 255
    if (compressed_countrate.min() < 0) or (compressed_countrate.max() > 255):
        raise ValueError(
            f'Exepected uncompressed counts {compressed_countrate} to in range 0-255')

    # TODO Must be a better way than creating entire lookup table on each call
    ll = np.arange(0, 16, 1)
    lkup = np.zeros(256, dtype='int')
    _sum = 0
    for i in range(0, 16):
        lkup[16 * i:16 * (i + 1)] = ll * 2 ** i + _sum
        if i < 15:
            _sum = lkup[16 * (i + 1) - 1] + 2 ** i

    return lkup[compressed_countrate]


def parse_observing_summary_hdulist(hdulist):
    """
    Parse a RHESSI observation summary file.

    Parameters
    ----------
    hdulist : `list`
        The HDU list from the fits file.

    Returns
    -------
    out : `dict`
        Returns a dictionary.
    """
    header = hdulist[0].header

    reference_time_ut = parse_time(hdulist[5].data.field('UT_REF')[0],
                                   format='utime')
    time_interval_sec = hdulist[5].data.field('TIME_INTV')[0]
    # label_unit = fits[5].data.field('DIM1_UNIT')[0]
    # labels = fits[5].data.field('DIM1_IDS')
    labels = ['3 - 6 keV', '6 - 12 keV', '12 - 25 keV', '25 - 50 keV',
              '50 - 100 keV', '100 - 300 keV', '300 - 800 keV',
              '800 - 7000 keV', '7000 - 20000 keV']

    # The data stored in the fits file are "compressed" countrates stored as
    # one byte
    compressed_countrate = np.array(hdulist[6].data.field('countrate'))

    countrate = uncompress_countrate(compressed_countrate)
    dim = np.array(countrate[:, 0]).size

    time_array = parse_time(reference_time_ut) + \
        TimeDelta(time_interval_sec * np.arange(dim) * u.second)

    #  TODO generate the labels for the dict automatically from labels
    data = {'time': time_array, 'data': countrate, 'labels': labels}

    return header, data


class RHESSISummaryTimeSeries(GenericTimeSeries):
    """
    RHESSI X-ray Summary lightcurve TimeSeries.

    The RHESSI mission consists of a single spin-stabilized spacecraft in a low-altitude orbit
    inclined 38 degrees to the Earth's equator.
    The only instrument on board is a set of 9 Germanium spectrometers with the ability to
    obtain high fidelity solar spectra from X rays (down to 3 keV) to gamma rays (1 MeV).
    Each spectrometer is coupled to a set of grids with different pitches which enable
    fourier-style imaging as the spacecraft spins.

    RHESSI provides summary lightcurves in the following passbands:

    * 3 - 6 keV
    * 6 - 12 keV
    * 12 - 25 keV
    * 25 - 50 keV
    * 50 - 100 keV
    * 100 - 300 keV
    * 300 - 800 keV
    * 800 - 7000 keV
    * 7000 - 20000 keV

    RHESSI was launched on 5th February 2002.

    Examples
    --------
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> import sunpy.timeseries
    >>> rhessi = sunpy.timeseries.TimeSeries(sunpy.data.sample.RHESSI_TIMESERIES)  # doctest: +REMOTE_DATA
    >>> rhessi.peek()  # doctest: +SKIP

    References
    ----------
    * `RHESSI Homepage. <https://hesperia.gsfc.nasa.gov/rhessi3/index.html>`_
    * `Mission Paper. <https://doi.org/10.1023/A:1022428818870>`_
    """

    # Class attribute used to specify the source class of the TimeSeries.
    _source = 'rhessi'

    def plot(self, axes=None, **kwargs):
        """
        Plots RHESSI Count Rate light curve from a pandas dataframe.

        Parameters
        ----------
        axes : `matplotlib.axes.Axes`, optional
            The axes on which to plot the TimeSeries. Defaults to current axes.
        **kwargs : `dict`
            Additional plot keyword arguments that are handed to `~matplotlib.axes.Axes.plot`
            functions.

        Returns
        -------
        `~matplotlib.axes.Axes`
            The plot axes.
        """
        self._validate_data_for_plotting()
        if axes is None:
            axes = plt.gca()
        # These are a matplotlib version of the default RHESSI color cycle
        default_colors = ('black', 'tab:pink', 'tab:green', 'tab:cyan',
                          'tab:olive', 'tab:red', 'tab:blue', 'tab:orange',
                          'tab:brown')
        colors = kwargs.pop('colors', default_colors)

        for color, (item, frame) in zip(itertools.cycle(colors),
                                        self.to_dataframe().items()):

            axes.plot(self.to_dataframe().index, frame.values,
                      color=color, label=item, **kwargs)

        axes.set_yscale("log")
        axes.set_xlabel(datetime.datetime.isoformat(self.to_dataframe().index[0])[0:10])
        axes.set_ylabel('Count Rate s$^{-1}$ detector$^{-1}$')
        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')
        axes.legend()
        # TODO: display better tick labels for date range (e.g. 06/01 - 06/05)
        formatter = matplotlib.dates.DateFormatter('%H:%M:%S')
        axes.xaxis.set_major_formatter(formatter)
        axes.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M:%S')
        return axes

    @peek_show
    def peek(self, title="RHESSI Observing Summary Count Rate", **kwargs):
        """
        Displays the RHESSI Count Rate light curve by calling
        `~sunpy.timeseries.sources.rhessi.RHESSISummaryTimeSeries.plot`.

        .. plot::

            import sunpy.data.sample
            import sunpy.timeseries
            rhessi = sunpy.timeseries.TimeSeries(sunpy.data.sample.RHESSI_TIMESERIES, source='RHESSI')
            rhessi.peek()

        Parameters
        ----------
        title : `str`
            The title of the plot. Defaults to "RHESSI Observing Summary Count Rate".
        **kwargs : `dict`
            Additional plot keyword arguments that are handed to `~matplotlib.axes.Axes.plot`
            functions.
        """
        fig, ax = plt.subplots()
        axes = self.plot(axes=ax, **kwargs)
        axes.set_title(title)
        fig.autofmt_xdate()
        return fig

    @classmethod
    def _parse_file(cls, filepath):
        """
        Parses rhessi FITS data files to create TimeSeries.

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
        Parses a RHESSI `astropy.io.fits.HDUList` from a FITS file.

        Parameters
        ----------
        hdulist : `astropy.io.fits.HDUList`
            A HDU list.
        """
        header, d = parse_observing_summary_hdulist(hdulist)
        # The time of dict `d` is astropy.time, but dataframe can only take datetime
        d['time'] = d['time'].datetime
        header = MetaDict(OrderedDict(header))
        data = DataFrame(d['data'], columns=d['labels'], index=d['time'])
        # Add the units data
        units = OrderedDict([('3 - 6 keV', u.ct / u.s / u.Unit('detector')),
                             ('6 - 12 keV', u.ct / u.s / u.Unit('detector')),
                             ('12 - 25 keV', u.ct / u.s / u.Unit('detector')),
                             ('25 - 50 keV', u.ct / u.s / u.Unit('detector')),
                             ('50 - 100 keV', u.ct / u.s / u.Unit('detector')),
                             ('100 - 300 keV', u.ct / u.s / u.Unit('detector')),
                             ('300 - 800 keV', u.ct / u.s / u.Unit('detector')),
                             ('800 - 7000 keV', u.ct / u.s / u.Unit('detector')),
                             ('7000 - 20000 keV', u.ct / u.s / u.Unit('detector'))])
        # Todo: check units used. https://hesperia.gsfc.nasa.gov/ssw/hessi/doc/guides/hessi_data_access.htm
        return data, header, units

    @classmethod
    def is_datasource_for(cls, **kwargs):
        """
        Determines if the file corresponds to a RHESSI X-ray Summary
        `~sunpy.timeseries.TimeSeries`.
        """
        # Check if source is explicitly assigned
        if 'source' in kwargs.keys():
            if kwargs.get('source', ''):
                return kwargs.get('source', '').lower().startswith(cls._source)
        # Check if HDU defines the source instrument
        if 'meta' in kwargs.keys():
            return kwargs['meta'].get('telescop', '').startswith('HESSI')
