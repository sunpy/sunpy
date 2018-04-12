# -*- coding: utf-8 -*-
"""RHESSI TimeSeries subclass definitions."""
from __future__ import absolute_import

from collections import OrderedDict
import datetime
import matplotlib.dates
import matplotlib.pyplot as plt
from pandas import DataFrame

from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.metadata import MetaDict
from sunpy.instr import rhessi
import sunpy.io
from astropy import units as u

__all__ = ['RHESSISummaryTimeSeries']


class RHESSISummaryTimeSeries(GenericTimeSeries):
    """
    RHESSI X-ray Summary Lightcurve TimeSeries.

    The RHESSI mission consists of a single spin-stabilized
    spacecraft in a low-altitude orbit inclined 38 degrees to
    the Earth's equator. The only instrument on board is a set of 9
    Germanium spectrometers with the ability to obtain high
    fidelity solar spectra from X rays (down to 3 keV) to gamma rays (1 MeV).
    Each spectrometer is coupled to a set of grids with different pitches
    which enable fourier-style imaging as the spacecraft spins.

    RHESSI provides summary lightcurves in the following passbands

    * 3 - 6 keV
    * 6 - 12 keV
    * 12 - 25 keV
    * 25 - 50 keV
    * 50 - 100 keV
    * 100 - 300 keV
    * 300 - 800 keV
    * 800 - 7000 keV
    * 7000 - 20000 keV

    RHESSI was launched on 5 February 2002.

    Examples
    --------
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> import sunpy.timeseries
    >>> rhessi = sunpy.timeseries.TimeSeries(sunpy.data.sample.RHESSI_TIMESERIES)  # doctest: +REMOTE_DATA
    >>> rhessi.peek()   # doctest: +SKIP

    References
    ----------
    * RHESSI Homepage `<https://hesperia.gsfc.nasa.gov/rhessi3/index.html>`_
    * Mission Paper `<https://doi.org/10.1023/A:1022428818870>`_
    """

    # Class attribute used to specify the source class of the TimeSeries.
    _source = 'rhessi'

    def peek(self, title="RHESSI Observing Summary Count Rate", **kwargs):
        """Plots RHESSI Count Rate light curve. An example is shown below.

        .. plot::

            import sunpy.data.sample
            import sunpy.timeseries
            rhessi = sunpy.timeseries.TimeSeries(sunpy.data.sample.RHESSI_TIMESERIES, source='RHESSI')
            rhessi.peek()

        Parameters
        ----------
        title : `str`
            The title of the plot.

        **kwargs : `dict`
            Any additional plot arguments that should be used
            when plotting.
        """
        # Check we have a timeseries valid for plotting
        self._validate_data_for_ploting()

        figure = plt.figure()
        axes = plt.gca()

        lc_linecolors = rhessi.hsi_linecolors()

        for lc_color, (item, frame) in zip(lc_linecolors, self.data.iteritems()):
            axes.plot_date(self.data.index, frame.values, '-', label=item, lw=2, color=lc_color)

        axes.set_yscale("log")
        axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])

        axes.set_title(title)
        axes.set_ylabel('Count Rate s$^{-1}$ detector$^{-1}$')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')
        axes.legend()

        # @todo: display better tick labels for date range (e.g. 06/01 - 06/05)
        formatter = matplotlib.dates.DateFormatter('%H:%M')
        axes.xaxis.set_major_formatter(formatter)

        axes.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
        figure.autofmt_xdate()
        figure.show()

    @classmethod
    def _parse_file(cls, filepath):
        """Parses rhessi FITS data files to create TimeSeries."""
        #header, d, hdus = rhessi.parse_obssumm_file(filepath)
        hdus = sunpy.io.read_file(filepath)
        return cls._parse_hdus(hdus)

    @classmethod
    def _parse_hdus(cls, hdulist):
        """Parses a RHESSI FITS HDU list form a FITS file."""
        header, d = rhessi.parse_obssumm_hdulist(hdulist)
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
        """Determines if the file corresponds to a RHESSI X-ray Summary lightcurve"""
        # Check if source is explicitly assigned
        if 'source' in kwargs.keys():
            if kwargs.get('source', ''):
                return kwargs.get('source', '').lower().startswith(cls._source)
        # Check if HDU defines the source instrument
        if 'meta' in kwargs.keys():
            return kwargs['meta'].get('telescop', '').startswith('HESSI')
