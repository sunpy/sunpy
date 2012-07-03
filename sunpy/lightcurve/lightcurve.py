"""
LightCurve is a generic LightCurve class from which all other LightCurve classes 
inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404,W0612,W0613
__authors__ = ["Keith Hughitt"]
__email__ = "keith.hughitt@nasa.gov"

import os
import datetime
import pandas
import sunpy
import urllib2
import numpy as np
import matplotlib.pyplot as plt

class LightCurve:
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
    data : pandas.DataFrame
        An pandas DataFrame prepresenting one or more fields as they vary with respect to time.
    header : string, dict
        The comment string or header associated with the light curve input

    Examples
    --------
    >>> import sunpy
    >>> import datetime
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> light_curve = sunpy.lightcurve.LightCurve({"param1": range(24 * 60)}, index=dates)
    >>> light_curve.show()

    References
    ----------
    | http://pandas.pydata.org/pandas-docs/dev/dsintro.html

    """
    def __init__(self, *args, **kwargs):
        self._filename = ""

        # If no arguments specified, perform default action
        if len(args) == 0:
            args = (self._get_default_uri(),)

        # Single argument
        if len(args) == 1:
            # If single string, could be a filepath, URL, date, or TimeRange
            if sunpy.time.is_time(args[0]):
                date = sunpy.time.parse_time(args[0])
                url = self._get_url_for_date(date)
                filepath = self._download(url, kwargs, err = "Unable to download data for specified date")
            elif isinstance(args[0], basestring):
                # Filepath
                if os.path.isfile(os.path.expanduser(args[0])):
                    filepath = os.path.expanduser(args[0])
                else:
                    # Otherwise, assume string input is a URL
                    try:
                        filepath = self._download(args[0], kwargs)
                    except (urllib2.HTTPError, urllib2.URLError, ValueError):
                        err = ("Unable to read location. Did you "
                               "specify a valid filepath or URL?")
                        raise ValueError(err)
            elif isinstance(args[0], sunpy.time.TimeRange):
                # TimeRange
                url = self._get_url_for_date_range(args[0])
                err = "Unable to download data for specified date range"
                filepath = self._download(url, err, kwargs)   

            # Parse resulting file
            header, data = self._parse_filepath(filepath)
        
        # Two arguments
        elif len(args) == 2:
            # Date range
            if (sunpy.time.is_time(args[0]) and sunpy.time.is_time(args[1])):
                url = self._get_url_for_date_range(args[0], args[1])
                filepath = self._download(url, kwargs, err = "Unable to download data for specified date range")
                header, data = self._parse_filepath(filepath)  
                
            # Other light curve creation options (DataFrame, ndarray, etc)
            elif isinstance(args[0], pandas.DataFrame):
                data = args[0]
            elif (isinstance(args[0], list) or 
                  isinstance(args[0], dict) or 
                  isinstance(args[0], np.ndarray) or 
                  isinstance(args[0], pandas.Series)):
            # DataFrame Index
                if "index" in kwargs:
                    index = kwargs["index"]
                else:
                    index = None
            else:
                raise TypeError("Both arguments passed are unrecognized.")

            data = pandas.DataFrame(args[0], index=index)
                
        # @NOTE: should we also support inputting start and end dates or a
        # date range?

        # Check for header
        #header = ""
        
        #if len(args) > 1:
        #    if (isinstance(args[1], basestring) or isinstance(args[1], dict)):
        #        header = args[1]
        #    else:
        #        raise TypeError("Unrecognized input for argument 2")

        if len(args) >= 2:
            raise TypeError("Lightcurve takes a maximum of two arguments.")

        self.data = data
        self.header = header
        
    def show(self, **kwargs):
        """Shows a plot of the light curve"""
        self.data.plot(**kwargs)
        plt.show()
        
    def _download(self, uri, kwargs, err='Unable to download data at specified URL'):
        """Attempts to download data at the specified URI"""
        self._filename = os.path.basename(uri).split("?")[0]
        
        # user specifies a download directory
        if "directory" in kwargs:
            download_dir = os.path.expanduser(kwargs["directory"])
        else:
            download_dir = sunpy.config.get("downloads", "download_dir")
        
        # overwrite the existing file if the keyword is present
        if "overwrite" in kwargs:
            overwrite = kwargs["overwrite"]
        else:
            overwrite = False
        
        # If the file is not already there, download it
        filepath = os.path.join(download_dir, self._filename)

        if not(os.path.isfile(filepath)) or (overwrite and os.path.isfile(filepath)): 
            try:
                response = urllib2.urlopen(uri)
            except (urllib2.HTTPError, urllib2.URLError):
                raise urllib2.URLError(err)
            fp = open(filepath, 'wb')
            fp.write(response.read())
        
        return filepath

    def _get_default_uri(self):
        """Default data to load when none is specified"""
        msg = "No default action set for %s"
        raise NotImplementedError(msg % self.__class__.__name__)
        
    def _get_url_for_date(self, date):
        """Returns a URL to the data for the specified date"""
        msg = "Date-based downloads not supported for for %s"
        raise NotImplementedError(msg % self.__class__.__name__)
    
    def _get_url_for_date_range(self, *args, **kwargs):
        """Returns a URL to the data for the specified date range"""
        msg = "Date-range based downloads not supported for for %s"
        raise NotImplementedError(msg % self.__class__.__name__)
    
    def _parse_csv(self, filepath):
        """Place holder method to parse CSV files."""
        pass
    
    def _parse_fits(self, filepath):
        """Place holder method to parse FITS files."""
        pass
    
    def _parse_filepath(self, filepath):
        filename, extension = os.path.splitext(filepath)
        
        if extension.lower() in (".csv", ".txt"):
            return self._parse_csv(filepath)
        else:
            return self._parse_fits(filepath)


    def discrete_boxcar_average(self, seconds=1):
        """Computes a discrete boxcar average for the DataFrame"""
        date_range = pandas.DateRange(self.data.index[0], self.data.index[-1], 
                                      offset=pandas.datetools.Second(seconds))
        grouped = self.data.groupby(date_range.asof)
        subsampled = grouped.mean()
        
        return LightCurve(subsampled, self.header.copy())
    
    def truncate(self, start=None, end=None):
        """Returns a truncated version of the Lyra object"""
        if start is None:
            start = self.data.index[0]
        if end is None:
            end = self.data.index[-1]
        
        truncated = self.data.truncate(sunpy.time.parse_time(start),
                                       sunpy.time.parse_time(end))
        return LightCurve(truncated, self.header.copy())
