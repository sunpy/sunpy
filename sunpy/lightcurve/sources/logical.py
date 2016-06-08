# -*- coding: utf-8 -*-
"""Provides a logical lightcurve.  Only two values are allowed - True or False.
Useful for keeping track of when an event occurred, usually labeled as
"True"."""
from __future__ import absolute_import

import numpy as np

from sunpy.lightcurve import LightCurve
from scipy.ndimage import label
from sunpy.time import TimeRange
import matplotlib.pyplot as plt
from sunpy.extern.six.moves import range

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

__all__ = ['LogicalLightCurve']


# Logical Lightcurve
# TODO
# Change the init to accept a list of TimeRange objects.  Durations between the
# start and end time of each TimeRange object are labeled 'True'.
class LogicalLightCurve(LightCurve):
    """
    Logical LightCurve with only True and False values.

    Examples
    --------
    >>> import sunpy.lightcurve as lightcurve
    >>> import datetime
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> z = [True for x in range(0, 24 * 60)]
    >>> light_curve = lightcurve.LogicalLightCurve.create({"param1": z}, index=dates)
    """

    @classmethod
    def _get_plot_types(cls):
        return ['logical']

    def plot(self, axes=None, title="Logical", plot_type=None, **plot_args):
        """Plots Logical lightcurve"""
        if axes is None:
            axes = plt.gca()

        if title is True:
            title = self.name

        if plot_type == None:
            plot_type = self._get_plot_types()[0]

        if plot_type == self._get_plot_types()[0]:      # logical
            self.data.plot(ax=axes, **plot_args)
            plt.fill_between(self.data.index, self.data['param1'], alpha=0.5)
            axes.set_title(title)
            axes.yaxis.grid(True, 'major')
            axes.xaxis.grid(True, 'major')
            axes.set_xlabel('Start time: ' + self.data.index[0].strftime(TIME_FORMAT))
        else:
            raise ValueError('Not a recognized plot type.')
        plt.gcf().autofmt_xdate()

        return axes

    def complement(self):
        """Return the logical complement of the original lightcurve."""
        return LogicalLightCurve.create(np.invert(self.data),
                                        header=self.header)

    def times(self):
        """Returns a list of time ranges where values are True.

        Returns
        -------
        outtr : `~sunpy.time.TimeRange` array
            An array of time ranges
        """

        labeling = label(self.data)
        timeranges = []
        for i in range(1, labeling[1] + 1):
            eventindices = (labeling[0] == i).nonzero()
            timeranges.append(TimeRange(self.data.index[eventindices[0][0]],
                                        self.data.index[eventindices[0][-1]]))
        return timeranges
