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
from sunpy.util.exceptions import warn_user
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

    Note that by default the data is re-binned from the original 128 into the following 8 pre-determined energy channels.
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

    Additionally, GBM timeseries with different energy bins ranges can be calculated by using the ``energy_bands`` keyword.

    Parameters
    ----------
    energy_bands : `list`, `tuple`
        The energy ranges to be calculated.

    Returns
    -------
    `sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries`
        with custom energy bins.

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

    def __init__(self, data, meta=None, units=None, **kwargs):
        super().__init__(data, meta, units, **kwargs)
        self.__header = None
        self.__energy_bins = None
        self.__count_data = None

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
    def _parse_file(cls, filepath, **kwargs):
        """
        Parses a GBM CSPEC FITS file.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        hdus = sunpy.io._file_tools.read_file(filepath)
        return cls._parse_hdus(hdus, **kwargs)

    @classmethod
    def _parse_hdus(cls, hdulist, **kwargs):
        """
        Parses a GBM CSPEC `astropy.io.fits.HDUList`.

        Parameters
        ----------
        hdulist : `str`
            The path to the file you want to parse.
        """
        if hdulist:
            cls.__header = MetaDict(OrderedDict(hdulist[0].header))
            # These GBM files have three FITS extensions.
            # extn1 - Gives the energy range for each of the 128 energy bins
            # extn2 - Contains the data, e.g. counts, exposure time, time of observation
            # extn3 - Eclipse times?
            cls.__energy_bins = hdulist[1].data
            cls.__count_data = hdulist[2].data
        else:
            cls.__header = kwargs.pop('header')
            cls.__energy_bins = kwargs.pop('energy_bins')
            cls.__count_data = kwargs.pop('count_data')

        ebands = kwargs.get('energy_bands', [4, 15, 25, 50, 100, 300, 800, 2000]*u.keV)
        for eband in ebands:
            if not eband.unit.is_equivalent(u.keV):
                unit_str, other_str = get_err_str(eband.unit), get_err_str(u.keV)
                raise u.UnitConversionError(f"{unit_str} and {other_str} are not convertible.")
        ebands = [eband.to(u.keV) if eband.unit != u.keV else eband for eband in ebands]
        ebands = sorted(ebands)

        # rebin the 128 energy channels into some summary ranges
        # some of default bin ranges are
        # 4-15 keV, 15 - 25 keV, 25-50 keV, 50-100 keV, 100-300 keV, 300-800 keV, 800 - 2000 keV
        # put the data in the units of counts/s/keV
        summary_counts, column_labels, ebands = _bin_data_for_summary(
            cls.__energy_bins, cls.__count_data, ebands
        )
        # get the time information in datetime format with the correct MET adjustment
        met_ref_time = parse_time('2001-01-01 00:00')  # Mission elapsed time
        gbm_times = met_ref_time + TimeDelta(cls.__count_data['time'], format='sec')
        gbm_times.precision = 9
        gbm_times = gbm_times.isot.astype('datetime64')

        # Add the units data
        units = OrderedDict((col, u.ct / u.s / eband.unit) for col, eband in zip(column_labels, ebands))

        return pd.DataFrame(summary_counts, columns=column_labels, index=gbm_times), cls.__header, units

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

    @classmethod
    def energy_bins(cls, energy_bands):
        """
        Determine GBM timeseries with different energy bins ranges.

        Parameters
        ----------
        energy_bands : `list`, `tuple`
            The energy ranges to be calculated.

        Returns
        -------
        `sunpy.timeseries.sources.fermi_gbm.GBMSummaryTimeSeries`
            with custom energy bins.

        Examples
        --------
        >>> import astropy.units as u
        >>> import sunpy.timeseries
        >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
        >>> gbm = sunpy.timeseries.TimeSeries(sunpy.data.sample.GBM_TIMESERIES, source='GBMSummary')  # doctest: +REMOTE_DATA
        >>> gbm_ts = gbm.energy_bins([50, 100, 200, 400, 1000] * u.keV)  # doctest: +SKIP
        """
        from sunpy.timeseries.timeseries_factory import TimeSeries
        data, meta, units = cls._parse_hdus(
            None,
            header = cls.__header,
            energy_bins = cls.__energy_bins,
            count_data = cls.__count_data,
            energy_bands = energy_bands
        )
        return TimeSeries(data, meta, units)


def _bin_data_for_summary(energy_bins, count_data, ebands):
    """
    Rebin the 128 energy channels into some summary ranges and put the data in
    the units of counts/s/keV.

    Default bin ranges used:
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
    ebands = ebands
    if len(ebands) <= 1:
        raise ValueError("'energy_bands' must contain more than one element.")
    e_center = (energy_bins['e_min'] + energy_bins['e_max']) / 2
    e_min, e_max = round(energy_bins['e_min'].min()), np.ceil(energy_bins['e_max'].max())

    if min(ebands).value >= e_min and max(ebands).value <= e_max:
        pass
    # ebands doesn't overlap with energy bins
    elif max(ebands).value <= e_min or min(ebands).value >= e_max:
        warn_user(
            "Data is not available for"
            f" {min(ebands)}-{max(ebands)}, "
            "updating energy bins to default value")
        ebands = [4, 15, 25, 50, 100, 300, 800, 2000] * u.keV
    else:
        warn_user("Only Partially Data is available for"
                f" {e_min} keV-{e_max} keV")
        # ebands fully overlap with energy bins
        if min(ebands).value < e_min and max(ebands).value > e_max:
            ebands = list(filter(lambda x: e_min <= x.value <= e_max, ebands))
            ebands.append(e_max * u.keV)
            ebands.insert(0, e_min * u.keV)

        # ebands partially overlap with energy bins
        elif min(ebands).value >= e_min and max(ebands).value > e_max:
            ebands = list(filter(lambda x: x.value <= e_max, ebands))
            ebands.append(e_max * u.keV)
        elif min(ebands).value < e_min and max(ebands).value <= e_max:
            ebands = list(filter(lambda x: x.value >= e_min, ebands))
            ebands.insert(0, e_min * u.keV)

    column_labels = [
        f"{eband.value}-{next_eband.value} {eband.unit}"
        for eband, next_eband in zip(ebands, ebands[1:])
    ]
    indices = [np.searchsorted(e_center, e.value) for e in ebands]
    summary_counts = []
    for ind_start, ind_end in zip(indices[:-1], indices[1:]):
        # sum the counts in the energy bands, and find counts/s/keV
        summed_counts = np.sum(count_data["counts"][:, ind_start:ind_end], axis=1)
        energy_width = (energy_bins["e_max"][ind_end - 1] - energy_bins["e_min"][ind_start])
        summary_counts.append(summed_counts/energy_width/count_data["exposure"])
    return np.array(summary_counts).T, column_labels, ebands


def get_err_str(unit: u.UnitBase) -> str:
    unit_str = unit.to_string("generic")
    physical_type = unit.physical_type
    if physical_type != "unknown":
        unit_str = f"'{unit_str}' ({physical_type})"
    else:
        unit_str = f"'{unit_str}'"
    return unit_str
