# -*- coding: utf-8 -*-
"""Provides programs to process and analyze RHESSI X-ray data."""
from __future__ import absolute_import

import datetime
import matplotlib.pyplot as plt
from pandas import DataFrame
from collections import OrderedDict

from sunpy.lightcurve import LightCurve
from sunpy.time import TimeRange, parse_time
from sunpy.instr import rhessi

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

__all__ = ['RHESSISummaryLightCurve']


class RHESSISummaryLightCurve(LightCurve):
    """
    RHESSI X-ray Summary LightCurve.

    The RHESSI mission consists of a single spin-stabilized
    spacecraft in a low-altitude orbit inclined 38 degrees to
    the Earth's equator. The only instrument on board is an
    Germaniun imaging spectrometer with the ability to obtain high
    fidelity solar images in X rays (down to 3 keV) to gamma rays (1 MeV).

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
    >>> from sunpy import lightcurve as lc
    >>> rhessi = lc.RHESSISummaryLightCurve.create()
    >>> rhessi = lc.RHESSISummaryLightCurve.create('2012/06/01', '2012/06/05')
    >>> rhessi.peek()   # doctest: +SKIP

    References
    ----------
    * RHESSI Homepage `<http://hesperia.gsfc.nasa.gov/rhessi3/index.html>`_
    * Mission Paper `<http://link.springer.com/article/10.1023%2FA%3A1022428818870>`_
    """

    def plot(self, title=True, axes=None, plot_type=None, **plot_args):
        """Plots RHESSI Count Rate light curve"""

        if axes is None:
            axes = plt.gca()

        if plot_type == None:
            plot_type = self._get_plot_types()[0]

        if title is True:
            title = self.name

        if plot_type == 'rhessi':
            for item, frame in self.data.iteritems():
                axes.plot_date(self.data.index, frame.values, '-', label=item, **plot_args)
            axes.set_yscale("log")
            axes.yaxis.grid(True, 'major')
            axes.xaxis.grid(True, 'major')
        else:
            raise ValueError('Not a recognized plot type.')

        axes.set_title(title)
        axes.set_ylabel(self.unit)
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime(TIME_FORMAT))
        plt.gcf().autofmt_xdate()

        return axes

    @classmethod
    def _get_plot_types(cls):
        return ['rhessi']

    @classmethod
    def _get_default_uri(cls):
        """Retrieves the latest RHESSI data."""
        today = datetime.datetime.today()
        days_back = 3
        time_range = TimeRange(today - datetime.timedelta(days=days_back),
                               today - datetime.timedelta(days=days_back - 1))
        return cls._get_url_for_date_range(time_range)

    @staticmethod
    def _get_url_for_date_range(*args, **kwargs):
        """Returns a URL to the RHESSI data for the specified date range.

        Parameters
        ----------
        args : `~sunpy.time.TimeRange`, `datetime.datetime, str
            Date range should be specified using a TimeRange, or start
            and end dates at datetime instances or date strings.
        """
        if len(args) == 1 and isinstance(args[0], TimeRange):
            time_range = args[0]
        elif len(args) == 2:
            time_range = TimeRange(parse_time(args[0]), parse_time(args[1]))
        url = rhessi.get_obssum_filename(time_range)
        return url

    @staticmethod
    def _parse_fits(filepath):
        """Parses a RHESSI FITS file"""
        header, d = rhessi.parse_obssumm_file(filepath)
        header = OrderedDict(header)
        data = DataFrame(d['data'], columns=d['labels'], index=d['time'])
        header.update({'UNIT': 'count s^-1 detector^-1'})
        header.update({'INSTRUME': 'RHESSI'})
        header.update({'OBSRVTRY': 'RHESSI'})
        header.update({'TELESCOP': 'RHESSI'})
        header.update({'WAVELNTH': [[3, 6], [6, 12], [12, 25], [25, 50], [50, 100], [100, 300],
                                    [300, 800], [800, 7000], [7000, 20000]]})
        header.update({'WAVEUNIT': 'keV'})
        return header, data
