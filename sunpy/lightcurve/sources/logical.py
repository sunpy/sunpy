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

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

__all__ = ['LogicalLightCurve']


# Logical Lightcurve
# TODO
# Change the init to accept a list of TimeRange objects.  Durations between the
# start and end time of each TimeRange object are labeled 'True'.
class LogicalLightCurve(LightCurve):
    """
    Logical LightCurve.

    Originated from a need to analyze the times of HEK
    results, where 'True' indicates an event was observed, and 'False'
    indicates an event was not observed.

    Examples
    --------
    >>> import sunpy.lightcurve as lightcurve
    >>> import datetime

    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> z = [True for x in range(0, 24 * 60)]
    >>> light_curve = lightcurve.LogicalLightCurve.create({"param1": z}, index=dates)
    """

    def plot(self, axes=None, title="Logical", **plot_args):
        """Plots RHESSI Count Rate light curve"""
        if axes is None:
            axes = plt.gca()
        self.data.plot(ax=axes, **plot_args)
        plt.fill_between(self.data.index, self.data['param1'], alpha=0.5)
        axes.set_title(title)
        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime(TIME_FORMAT))

        axes.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.gcf().autofmt_xdate()
        
        return axes

    def complement(self):
        """ Define the complement of the passed lightcurve """
        return LogicalLightCurve.create(np.invert(self.data),
                                        header = self.header)

    def times(self):
        """Label all the periods of time that have the value 'True'. Return
        a list of TimeRange objects """

        labeling = label(self.data)
        timeranges = []
        for i in xrange(1, labeling[1]+1):
            eventindices = (labeling[0] == i).nonzero()
            timeranges.append( TimeRange(self.data.index[ eventindices[0][0] ],
                                         self.data.index[ eventindices[0][-1] ]) )
        return timeranges