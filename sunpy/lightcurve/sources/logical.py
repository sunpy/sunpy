# -*- coding: utf-8 -*-
"""Provides a logical lightcurve.  Only two values are allowed - True or False.
Useful for keeping track of when an event occurred, usually labeled as 
"True"."""
from __future__ import absolute_import

import numpy as np

from sunpy.lightcurve import LightCurve
from scipy.ndimage import label
from sunpy.time import TimeRange

__all__ = ['LogicalLightCurve']

#
#
# Logical Lightcurve
# TODO
# Change the init to accept a list of TimeRange objects.  Durations between the
# start and end time of each TimeRange object are labeled 'True'.
class LogicalLightCurve(LightCurve):
    """
    Logical light curve.  Originated from a need to analyze the times of HEK
    results, where 'True' indicates an event was observed, and 'False' indicates 
    an event was not observed.

    Examples
    --------
    import sunpy
    light_curve = sunpy.lightcurve.LogicalLightCurve.create()
    """

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
