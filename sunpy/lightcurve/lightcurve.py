"""
LightCurve is a generic LightCurve class from which all other LightCurve classes
inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404,W0612,W0613
__authors__ = ["Keith Hughitt"]
__email__ = "keith.hughitt@nasa.gov"

import os.path
import shutil
import urllib2
import warnings
from datetime import datetime
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
import pandas

from sunpy import config
from sunpy.time import is_time, TimeRange, parse_time
from sunpy.util.cond_dispatch import ConditionalDispatch, run_cls

__all__ = ['LightCurve']


class LightCurve(object):
    """
    LightCurve(filepath)

    A generic light curve object.

    Attributes
    ----------

    meta : `str` or `dict`
        The comment string or header associated with the data.
    data : `~pandas.DataFrame`
        An pandas DataFrame prepresenting one or more fields as a function of time.

    Examples
    --------
    >>> import sunpy
    >>> import datetime
    >>> import numpy as np
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, step=(12 * np.pi) / 24 * 60))
    >>> light_curve = sunpy.lightcurve.LightCurve.create({"param1": intensity}, index=dates)
    >>> light_curve.peek()

    References
    ----------
    * `Pandas Documentation <http://pandas.pydata.org/pandas-docs/dev/dsintro.html>`_

    """
    _cond_dispatch = ConditionalDispatch()
    create = classmethod(_cond_dispatch.wrapper())

    def __init__(self, data, meta=None):
        self.data = pandas.DataFrame(data)
        if meta == '' or meta is None:
            self.meta = OrderedDict()
        else:
            self.meta = OrderedDict(meta)

    @property
    def header(self):
        """
        Return the lightcurves metadata

        .. deprecated:: 0.4.0
            Use .meta instead
        """
        warnings.warn("""lightcurve.header has been renamed to lightcurve.meta
for compatibility with map, please use meta instead""", Warning)
        return self.meta

    @classmethod
    def from_time(cls, time, **kwargs):
        """
        Called by Conditional Dispatch object when valid time is passed as
        input to create method.
        """
        date = parse_time(time)
        url = cls._get_url_for_date(date, **kwargs)
        filepath = cls._download(
            url, kwargs, err="Unable to download data for specified date"
        )
        return cls.from_file(filepath)

    @classmethod
    def from_range(cls, start, end, **kwargs):
        """Called by Conditional Dispatch object when start and end time are
        passed as input to create method.

        :param start:
        :param end:
        :param kwargs:
        :return:
        """
        url = cls._get_url_for_date_range(parse_time(start), parse_time(end), **kwargs)
        filepath = cls._download(
            url, kwargs,
            err="Unable to download data for specified date range"
        )
        result = cls.from_file(filepath)
        result.data = result.data.truncate(start,end)
        return result

    @classmethod
    def from_timerange(cls, timerange, **kwargs):
        """
        Called by Conditional Dispatch object when time range is passed as
        input to create method.
        """
        url = cls._get_url_for_date_range(timerange, **kwargs)
        filepath = cls._download(
            url, kwargs,
            err = "Unable to download data for specified date range"
        )
        result = cls.from_file(filepath)
        result.data = result.data.truncate(timerange.start, timerange.end)
        return result

    @classmethod
    def from_file(cls, filename):
        """Used to return Light Curve object by reading the given filename.

        Parameters
        ----------
        filename: `str`
            Path of the file to be read.

        Returns
        -------
        Lightcurve object.
        """
        filename = os.path.expanduser(filename)
        meta, data = cls._parse_filepath(filename)
        if data.empty:
            raise ValueError("No data found!")
        else:
            return cls(data, meta)

    @classmethod
    def from_url(cls, url, **kwargs):
        """
        Called by Conditional Dispatch object to create Light Curve object when
        given a url. Downloads a file from the given url, attemps to read it
        and returns a Light Curve object.

        Parameters
        ----------
        url : str
            A url given as a string.
        """
        try:
            filepath = cls._download(url, kwargs)
        except (urllib2.HTTPError, urllib2.URLError, ValueError):
            err = "Unable to read location {!s}.".format(url)
            raise ValueError(err)
        return cls.from_file(filepath)

    @classmethod
    def from_data(cls, data, index=None, meta=None):
        """
        Called by Conditional Dispatch object to create Light Curve object when
        corresponding data is passed to create method.

        Parameters
        ----------
        data : `~numpy.ndarray`
            The data array
        index : `~datetime.datetime` array
            The time values
        """

        return cls(
            pandas.DataFrame(data, index=index),
            meta
        )

    @classmethod
    def from_yesterday(cls):
        """
        Called by Conditional Dispatch object if no input if given
        """
        return cls.from_url(cls._get_default_uri())

    @classmethod
    def from_dataframe(cls, dataframe, meta=None):
        """
        Called by Conditional Dispatch object to create Light Curve object when
        Pandas DataFrame is passed to create method.

        Parameters
        ----------
        dataframe : `~pandas.DataFrame`
            The data.
        meta : `str` or `dict`
            The metadata.
        """

        return cls(dataframe, meta)

    def plot(self, axes=None, **plot_args):
        """Plot a plot of the light curve

        Parameters
        ----------
        axes : `~matplotlib.axes.Axes` or None
            If provided the image will be plotted on the given axes. Otherwise
            the current axes will be used.

        **plot_args : `dict`
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        axes : `~matplotlib.axes.Axes`
            The plot axes.
        """

        # Get current axes
        if axes is None:
            axes = plt.gca()

        axes = self.data.plot(ax=axes, **plot_args)

        return axes

    def peek(self, **kwargs):
        """Displays the light curve in a new figure.

        Parameters
        ----------
        **kwargs : `dict`
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """

        figure = plt.figure()
        self.plot(**kwargs)
        figure.show()

        return figure

    @staticmethod
    def _download(uri, kwargs,
                  err='Unable to download data at specified URL'):
        """Attempts to download data at the specified URI.

        Parameters
        ----------
        **kwargs : uri
            A url
        """

        _filename = os.path.basename(uri).split("?")[0]

        # user specifies a download directory
        if "directory" in kwargs:
            download_dir = os.path.expanduser(kwargs["directory"])
        else:
            download_dir = config.get("downloads", "download_dir")

        # overwrite the existing file if the keyword is present
        if "overwrite" in kwargs:
            overwrite = kwargs["overwrite"]
        else:
            overwrite = False

        # If the file is not already there, download it
        filepath = os.path.join(download_dir, _filename)

        if not(os.path.isfile(filepath)) or (overwrite and
                                             os.path.isfile(filepath)):
            try:
                response = urllib2.urlopen(uri)
            except (urllib2.HTTPError, urllib2.URLError):
                raise urllib2.URLError(err)
            with open(filepath, 'wb') as fp:
                shutil.copyfileobj(response, fp)
        else:
            warnings.warn("Using existing file rather than downloading, use overwrite=True to override.", RuntimeWarning)

        return filepath

    @classmethod
    def _get_default_uri(cls):
        """Default data to load when none is specified."""
        msg = "No default action set for {}"
        raise NotImplementedError(msg.format(cls.__name__))

    @classmethod
    def _get_url_for_date(cls, date, **kwargs):
        """Returns a URL to the data for the specified date."""
        msg = "Date-based downloads not supported for for {}"
        raise NotImplementedError(msg.format(cls.__name__))

    @classmethod
    def _get_url_for_date_range(cls, *args, **kwargs):
        """Returns a URL to the data for the specified date range."""
        msg = "Date-range based downloads not supported for for {}"
        raise NotImplementedError(msg.format(cls.__name__))

    @staticmethod
    def _parse_csv(filepath):
        """Place holder method to parse CSV files."""
        msg = "Generic CSV parsing not yet implemented for LightCurve"
        raise NotImplementedError(msg)

    @staticmethod
    def _parse_fits(filepath):
        """Place holder method to parse FITS files."""
        msg = "Generic FITS parsing not yet implemented for LightCurve"
        raise NotImplementedError(msg)

    @classmethod
    def _parse_filepath(cls, filepath):
        """Check the file extension to see how to parse the file."""
        filename, extension = os.path.splitext(filepath)

        if extension.lower() in (".csv", ".txt"):
            return cls._parse_csv(filepath)
        else:
            return cls._parse_fits(filepath)

    def truncate(self, a, b=None):
        """Returns a truncated version of the timeseries object.

        Parameters
        ----------
        a : `sunpy.time.TimeRange`
            A time range to truncate to.

        Returns
        -------
        newlc : `~sunpy.lightcurve.LightCurve`
            A new lightcurve with only the selected times.
        """
        if isinstance(a, TimeRange):
            time_range = a
        else:
            time_range = TimeRange(a,b)

        truncated = self.data.truncate(time_range.start, time_range.end)
        return self.__class__.create(truncated, self.meta.copy())

    def extract(self, column_name):
        """Returns a new lightcurve with the chosen column.

        Parameters
        ----------
        column_name : `str`
            A valid column name

        Returns
        -------
        newlc : `~sunpy.lightcurve.LightCurve`
            A new lightcurve with only the selected column.
        """
        # TODO allow the extract function to pick more than one column
        if isinstance(self, pandas.Series):
            return self
        else:
            return LightCurve(self.data[column_name], self.meta.copy())

    def time_range(self):
        """Returns the start and end times of the LightCurve as a `~sunpy.time.TimeRange`
        object"""
        return TimeRange(self.data.index[0], self.data.index[-1])

# What's happening here is the following: The ConditionalDispatch is just an
# unbound callable object, that is, it does not know which class it is attached
# to. What we do against that is return a wrapper and make that a classmethod -
# thus we get the class as whose member it is called as as the first argument,
# this is why in the type signature we always have type as the first type.

# We then use run_cls, which just returns a wrapper that interprets the first
# argument as the class the function should be called of. So,
# x = run_cls("foo") returns something that turns x(cls, 1) into cls.foo(1).
# Because this has *args, **kwargs as its signature we need to disable the
# check of ConditionalDispatch that makes sure the function and the
# conditional need to have the same signature - but they still do have to.

LightCurve._cond_dispatch.add(
    run_cls("from_time"),
    lambda cls, time, **kwargs: is_time(time),
    # type is here because the class parameter is a class,
    # i.e. an instance of type (which is the base meta-class).
    [type, (basestring, datetime, tuple)],
    False
)

LightCurve._cond_dispatch.add(
    run_cls("from_range"),
    lambda cls, time1, time2, **kwargs: is_time(time1) and is_time(time2),
    [type, (basestring, datetime, tuple), (basestring, datetime, tuple)],
    False
)

LightCurve._cond_dispatch.add(
    run_cls("from_timerange"),
    lambda cls, timerange, **kwargs: True,
    [type, TimeRange],
    False
)

LightCurve._cond_dispatch.add(
    run_cls("from_file"),
    lambda cls, filename: os.path.exists(os.path.expanduser(filename)),
    [type, basestring],
    False
)

LightCurve._cond_dispatch.add(
    run_cls("from_url"),
    lambda cls, url, **kwargs: True,
    [type, basestring],
    False
)

LightCurve._cond_dispatch.add(
    run_cls("from_data"),
    lambda cls, data, index=None, meta=None: True,
    [type, (list, dict, np.ndarray, pandas.Series), object, object],
    False
)

LightCurve._cond_dispatch.add(
    run_cls("from_dataframe"),
    lambda cls, dataframe, meta=None: True,
    [type, pandas.DataFrame, object],
    False
)

LightCurve._cond_dispatch.add(
    run_cls("from_yesterday"),
    lambda cls: True,
    [type],
    False
)
