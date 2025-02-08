"""
This module FERMI GBM `~sunpy.timeseries.TimeSeries` source.
"""
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import astropy.units as u
from astropy.time import TimeDelta

import sunpy.io
import sunpy.io._file_tools
from sunpy.time import parse_time
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.metadata import MetaDict
from sunpy.visualization import peek_show

__all__ = ['GBMSummaryTimeSeries']


class GBMSummaryTimeSeries(GenericTimeSeries):
    """
    Fermi/GBM Summary lightcurve TimeSeries.

    The Gamma-ray Burst Monitor (GBM) is an instrument on board Fermi.
    It is meant to detect gamma-ray bursts but also detects solar flares.
    It consists of 12 Sodium Iodide (NaI) scintillation detectors and 2 Bismuth Germanate (BGO) scintillation detectors.
    The NaI detectors cover from a few keV to about 1 MeV and provide burst triggers and locations.
    The BGO detectors cover the energy range from about 150 keV to about 30 MeV.

    This summary lightcurve makes use of the CSPEC (daily version) data set which consists of the counts
    accumulated every 4.096 seconds in 128 energy channels for each of the 14 detectors.

    Note that the data is re-binned from the original 128 into the following 8 pre-determined energy channels.
    The rebinning method treats the counts in each of the original 128 channels as
    all having the energy of the average energy of that channel. For example, the
    counts in an 14.5--15.6 keV original channel would all be accumulated into the
    15--25 keV rebinned channel.

    * 4-15 keV
    * 15-25 keV
    * 25-50 keV
    * 50-100 keV
    * 100-300 keV
    * 300-800 keV
    * 800-2000 keV

    Examples
    --------
    >>> import sunpy.timeseries
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> gbm = sunpy.timeseries.TimeSeries(sunpy.data.sample.GBM_TIMESERIES, source='GBMSummary')  # doctest: +REMOTE_DATA
    >>> gbm.peek()   # doctest: +SKIP

    References
    ----------
    * `Fermi Mission Homepage <https://fermi.gsfc.nasa.gov>`__
    * `Fermi GBM Homepage <https://fermi.gsfc.nasa.gov/science/instruments/gbm.html>`__
    * `Fermi Science Support Center <https://fermi.gsfc.nasa.gov/ssc/>`__
    * `Fermi Data Product <https://fermi.gsfc.nasa.gov/ssc/data/access/>`__
    * `GBM Instrument Papers <https://gammaray.nsstc.nasa.gov/gbm/publications/instrument_journal_gbm.html>`__
    """
    # Class attributes used to specify the source class of the TimeSeries
    # and a URL to the mission website.
    _source = 'gbmsummary'
    _url = "https://gammaray.nsstc.nasa.gov/gbm/#"

    def plot(self, axes=None, columns=None, **kwargs):
        """
        Plots the GBM timeseries.

        Parameters
        ----------
        axes : `matplotlib.axes.Axes`, optional
            The axes on which to plot the TimeSeries. Defaults to current axes.
        columns : list[str], optional
            If provided, only plot the specified columns.
        **kwargs : `dict`
            Additional plot keyword arguments that are handed to `~matplotlib.axes.Axes.plot`
            functions.

        Returns
        -------
        `~matplotlib.axes.Axes`
            The plot axes.
        """
        axes, columns = self._setup_axes_columns(axes, columns)
        for d in columns:
            axes.plot(self._data.index, self._data[d], label=d, **kwargs)
        axes.set_yscale("log")
        axes.set_ylabel('Counts/s/keV')
        axes.legend()
        self._setup_x_axis(axes)
        return axes

    @peek_show
    def peek(self, *, title=None, columns=None, **kwargs):
        """
        Displays the GBM timeseries by calling
        `~sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries.plot`.

        .. plot::

            import sunpy.timeseries
            import sunpy.data.sample
            gbm = sunpy.timeseries.TimeSeries(sunpy.data.sample.GBM_TIMESERIES, source='GBMSummary')
            gbm.peek()

        Parameters
        ----------
        title : `str`, optional
            The title of the plot.
        columns : list[str], optional
            If provided, only plot the specified columns.
        **kwargs : `dict`
            Additional plot keyword arguments that are handed to `~matplotlib.axes.Axes.plot`
            functions.
        """
        if title is None:
            title = 'Fermi GBM Summary data ' + str(self.meta.get('DETNAM').values())
        fig, ax = plt.subplots()
        axes = self.plot(axes=ax, columns=columns, **kwargs)
        axes.set_title(title)
        return fig

    @classmethod
    def _parse_file(cls, filepath):
        """
        Parses a GBM CSPEC FITS file.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        hdus = sunpy.io._file_tools.read_file(filepath)
        return cls._parse_hdus(hdus)

    @classmethod
    def _parse_hdus(cls, hdulist):
        """
        Parses a GBM CSPEC `astropy.io.fits.HDUList`.

        Parameters
        ----------
        hdulist : `str`
            The path to the file you want to parse.
        """
        header = MetaDict(OrderedDict(hdulist[0].header))
        # these GBM files have three FITS extensions.
        # extn1 - this gives the energy range for each of the 128 energy bins
        # extn2 - this contains the data, e.g. counts, exposure time, time of observation
        # extn3 - eclipse times?
        energy_bins = hdulist[1].data
        count_data = hdulist[2].data

        # rebin the 128 energy channels into some summary ranges
        # 4-15 keV, 15 - 25 keV, 25-50 keV, 50-100 keV, 100-300 keV, 300-800 keV, 800 - 2000 keV
        # put the data in the units of counts/s/keV
        summary_counts = _bin_data_for_summary(energy_bins, count_data)

        # get the time information in datetime format with the correct MET adjustment
        met_ref_time = parse_time('2001-01-01 00:00')  # Mission elapsed time
        gbm_times = met_ref_time + TimeDelta(count_data['time'], format='sec')
        gbm_times.precision = 9
        gbm_times = gbm_times.isot.astype('datetime64')

        column_labels = ['4-15 keV', '15-25 keV', '25-50 keV', '50-100 keV',
                         '100-300 keV', '300-800 keV', '800-2000 keV']

        # Add the units data
        units = OrderedDict([('4-15 keV', u.ct / u.s / u.keV), ('15-25 keV', u.ct / u.s / u.keV),
                             ('25-50 keV', u.ct / u.s / u.keV), ('50-100 keV', u.ct / u.s / u.keV),
                             ('100-300 keV', u.ct / u.s / u.keV), ('300-800 keV', u.ct / u.s / u.keV),
                             ('800-2000 keV', u.ct / u.s / u.keV)])
        return pd.DataFrame(summary_counts, columns=column_labels, index=gbm_times), header, units

    @classmethod
    def is_datasource_for(cls, **kwargs):
        """
        Determines if the file corresponds to a GBM summary lightcurve
        `~sunpy.timeseries.TimeSeries`.
        """
        # Check if source is explicitly assigned
        if 'source' in kwargs.keys():
            if kwargs.get('source', ''):
                return kwargs.get('source', '').lower().startswith(cls._source)
        # Check if HDU defines the source instrument
        if 'meta' in kwargs.keys():
            return kwargs['meta'].get('INSTRUME', '').startswith('GBM')


def _bin_data_for_summary(energy_bins, count_data):
    """
    Rebin the 128 energy channels into some summary ranges and put the data in
    the units of counts/s/keV.

    Bin ranges used:
    * 4-15 keV
    * 15-25 keV
    * 25-50 keV
    * 50-100 keV
    * 100-300 keV
    * 300-800 keV
    * 800-2000 keV

    Parameters
    ----------
    energy_bins : `numpy.ndarray`
        The array of energy bins to rebin.
    count_data : `numpy.ndarray`
        The array of count data to rebin.
    """

    # list of energy bands to sum between
    ebands = [4, 15, 25, 50, 100, 300, 800, 2000]
    e_center = (energy_bins['e_min'] + energy_bins['e_max']) / 2
    indices = [np.searchsorted(e_center, e) for e in ebands]

    summary_counts = []
    for ind_start, ind_end in zip(indices[:-1], indices[1:]):
        # sum the counts in the energy bands, and find counts/s/keV
        summed_counts = np.sum(count_data["counts"][:, ind_start:ind_end], axis=1)
        energy_width = (energy_bins["e_max"][ind_end - 1] - energy_bins["e_min"][ind_start])
        summary_counts.append(summed_counts/energy_width/count_data["exposure"])

    return np.array(summary_counts).T
