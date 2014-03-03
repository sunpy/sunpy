"""
LightCurve is a generic LightCurve class from which all other LightCurve classes
inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404,W0612,W0613
__authors__ = ["Keith Hughitt"]
__email__ = "keith.hughitt@nasa.gov"

import os
import shutil
import urllib2
import warnings
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import pandas

import sunpy
from sunpy.time import is_time, TimeRange, parse_time

__all__ = ['GenericLightCurve']

class GenericLightCurve(object):
    """
    LightCurve(filepath)

    A generic light curve object.

    Parameters
    ----------
    args : filepath, url, or start and end dates
        The input for a LightCurve object should either be a filepath, a URL,
        or a date range to be queried for the particular instrument.

    Attributes
    ----------
    meta : string, dict
        The comment string or header associated with the light curve input
    data : pandas.DataFrame
        An pandas DataFrame prepresenting one or more fields as they vary with
        respect to time.

    Examples
    --------
    >>> import sunpy
    >>> import datetime
    >>> import numpy as np

    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]

    >>> intensity = np.sin(np.arange(0, 12 * np.pi, step=(12 * np.pi) / 24 * 60))

    >>> light_curve = sunpy.lightcurve.LightCurve.create(
    ...    {"param1": intensity}, index=dates
    ... )

    >>> light_curve.peek()

    References
    ----------
    | http://pandas.pydata.org/pandas-docs/dev/dsintro.html

    """

    def __init__(self, data, meta=None):
        self.data = pandas.DataFrame(data)
        self.meta = meta

    def __add__(self, other):
        """
        List like concatenation for lightcurves
        """
        if not isinstance(other, GenericLightCurve):
            return NotImplemented

        self.data.append(other.data)
        self.meta.update(other.meta) #TODO: make this MUCH smarter?

    @property
    def header(self):
        """
        Return the lightcurves metadata

        .. deprecated:: 0.4.0
            Use .meta instead
        """
        warnings.warn("""lightcurve.header has been renamed to lightcurve.meta
for compatability with map, please use meta instead""", Warning)
        return self.meta

    def plot(self, axes=None, **plot_args):
        """Plot a plot of the light curve

        Parameters
        ----------
        axes: matplotlib.axes object or None
            If provided the image will be plotted on the given axes. Else the
            current matplotlib axes will be used.

        **plot_args : dict
            Any additional plot arguments that should be used
            when plotting the image.

        """

        #Get current axes
        if axes is None:
            axes = plt.gca()

        axes = self.data.plot(ax=axes, **plot_args)

        return axes

    def peek(self, **kwargs):
        """Displays the light curve in a new figure"""

        figure = plt.figure()

        self.plot(**kwargs)

        figure.show()

        return figure

    @classmethod
    def _get_url_for_date_range(cls, *args, **kwargs):
        """Returns a URL to the data for the specified date range"""
        msg = "Date-range based downloads not supported for for %s"
        raise NotImplementedError(msg % cls.__name__)

    def truncate(self, a, b=None):
        """Returns a truncated version of the timeseries object"""
        if isinstance(a, TimeRange):
            time_range = a
        else:
            time_range = TimeRange(a,b)

        truncated = self.data.truncate(time_range.start(), time_range.end())
        return GenericLightCurve(truncated, self.meta.copy())

    def extract(self, a):
        """Extract a set of particular columns from the DataFrame"""
        # TODO allow the extract function to pick more than one column
        if isinstance(self, pandas.Series):
            return self
        else:
            return GenericLightCurve(self.data[a], self.meta.copy())

    def time_range(self):
        """Returns the start and end times of the LightCurve as a TimeRange
        object"""
        return TimeRange(self.data.index[0], self.data.index[-1])