"""
This module FERMI GBM `~sunpy.timeseries.TimeSeries` source.
"""
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import astropy.units as u
from astropy.time import Time

import sunpy.io
from sunpy.instr import fermi
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.metadata import MetaDict
from sunpy.visualization import peek_show

__all__ = ['GBMSummaryTimeSeries']


class GBMSummaryTimeSeries(GenericTimeSeries):
    """
    Fermi/GBM Summary Lightcurve TimeSeries.

    The Gamma-ray Burst Monitor (GBM) is an instrument aboard Fermi.
    It is meant to detect gamma-ray bursts but also detects solar flares.
    It consists of 12 Sodium Iodide (NaI) scintillation detectors and 2 Bismuth Germanate (BGO) scintillation detectors.
    The NaI detectors cover from a few keV to about 1 MeV and provide burst triggers and locations.
    The BGO detectors cover the energy range from about 150 keV to about 30 MeV.

    This summary lightcurve makes use of the CSPEC (daily version) data set which consists of the counts
    accumulated every 4.096 seconds in 128 energy channels for each of the 14 detectors.
    Note that the data is re-binned from the original 128 into the following 8 pre-determined energy channels.

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
    * `Fermi Mission Homepage <https://fermi.gsfc.nasa.gov>`_
    * `Fermi GBM Homepage <https://fermi.gsfc.nasa.gov/science/instruments/gbm.html>`_
    * `Fermi Science Support Center <https://fermi.gsfc.nasa.gov/ssc/>`_
    * `Fermi Data Product <https://fermi.gsfc.nasa.gov/ssc/data/access/>`_
    * `GBM Instrument Papers <https://gammaray.nsstc.nasa.gov/gbm/publications/instrument_journal_gbm.html>`_
    """
    # Class attribute used to specify the source class of the TimeSeries.
    _source = 'gbmsummary'

    @peek_show
    def peek(self):
        """
        Plots the GBM lightcurve TimeSeries. An example can be seen below:

        .. plot::

            import sunpy.timeseries
            import sunpy.data.sample
            gbm = sunpy.timeseries.TimeSeries(sunpy.data.sample.GBM_TIMESERIES, source='GBMSummary')
            gbm.peek()
        """
        # Check we have a timeseries valid for plotting
        self._validate_data_for_ploting()

        figure = plt.figure()
        axes = plt.gca()
        data_lab = self.data.columns.values

        for d in data_lab:
            axes.plot(self.data.index, self.data[d], label=d)

        axes.set_yscale("log")
        axes.set_title('Fermi GBM Summary data ' + str(self.meta.get(
            'DETNAM').values()))
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime(
            '%Y-%m-%d %H:%M:%S UT'))
        axes.set_ylabel('Counts/s/keV')
        axes.legend()
        figure.autofmt_xdate()

        return figure

    @classmethod
    def _parse_file(cls, filepath):
        """
        Parses a GBM CSPEC FITS file.

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
        Parses a GBM CSPEC `astropy.io.fits.HDUList`.

        Parameters
        ----------
        filepath : `str`
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
        gbm_times = Time([fermi.met_to_utc(t) for t in count_data['time']])
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
    # find the indices corresponding to some standard summary energy bins
    ebands = [4, 15, 25, 50, 100, 300, 800, 2000]
    indices = []
    for e in ebands:
        indices.append(np.searchsorted(energy_bins['e_max'], e))

    summary_counts = []
    for i in range(0, len(count_data['counts'])):
        counts_in_bands = []
        for j in range(1, len(ebands)):
            counts_in_bands.append(
                np.sum(count_data['counts'][i][indices[j - 1]:indices[j]]) /
                (count_data['exposure'][i] *
                 (energy_bins['e_max'][indices[j]] -
                  energy_bins['e_min'][indices[j - 1]])))

        summary_counts.append(counts_in_bands)

    return summary_counts


def _parse_detector(detector):
    """
    Check and fix detector name strings.

    Parameters
    ----------
    detector : `str`
        The detector name to check.
    """
    oklist = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9',
              'n10', 'n11']
    altlist = [str(i) for i in range(12)]
    if detector in oklist:
        return detector
    elif detector in altlist:
        return 'n' + detector
    else:
        raise ValueError('Detector string could not be interpreted')
