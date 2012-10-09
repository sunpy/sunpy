# -*- coding: utf-8 -*-
"""Provides a logical lightcurve. Useful for keeping track of the duration
of events."""
from __future__ import absolute_import
from sunpy.lightcurve import LightCurve
from scipy.ndimage import label

#
#
# Logical Lightcurve
#
class LogicalLightCurve(LightCurve):
    """
    Logical light curve.  Originated from a need to analyze the times of HEK
    results, where 'True' indicates an event was observed, and 'False' indicates 
    an event was not observed."""
        
    def show(self, **kwargs):
        """Shows a plot of the light curve"""
        fig = self.plot(**kwargs)
        fig.show()
        
        return fig

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
