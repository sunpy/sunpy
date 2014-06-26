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

__all__ = ['LogicalLightCurve']

#TODO Change the init to accept a list of TimeRange objects.  Durations
# between the start and end time of each TimeRange object are labeled 'True'.

class LogicalLightCurve(LightCurve):
    """
    Logical LightCurve.

    A lightcurve which contains only True or False values.

    Examples
    --------
    >>> import sunpy.lightcurve as lightcurve
    >>> import datetime
    # Create a logical lightcurve which is true every minute divisible by 5
    # in the last hour
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(seconds=x) for x in range(0, 1 * 60 * 60)]
    >>> z = [((d.minute % 5) == 0) for d in dates]
    >>> light_curve = lightcurve.LogicalLightCurve.create({"param1": z}, index=dates)
    """

    def plot(self, axes=None, title="Logical", **plot_args):
        """Plot the logical lightcurve"""
        if axes is None:
            axes = plt.gca()

        axes = self.data.plot(ax=axes, **plot_args)
        axes.fill_between(self.data.index,
                          self.data[self.data.columns[0]].values, alpha=0.5)
        axes.set_title(title)

        return axes

    def complement(self):
        """ Define the complement of the passed lightcurve """
        return LogicalLightCurve.create(np.invert(self.data),
                                        header=self.header)

    def times(self):
        """Label all the periods of time that have the value 'True'. Return
        a list of TimeRange objects """

        labeling = label(self.data)
        timeranges = []
        for i in xrange(1, labeling[1] + 1):
            eventindices = (labeling[0] == i).nonzero()
            timeranges.append(TimeRange(self.data.index[eventindices[0][0]],
                                         self.data.index[eventindices[0][-1]]))
        return timeranges
