"""
LightCurve is a generic LightCurve class from which all other LightCurve classes 
inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404,W0612,W0613
__authors__ = ["Keith Hughitt"]
__email__ = "keith.hughitt@nasa.gov"

import os
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
        if len(args) is 0:
            args = (self._get_default_uri(),)

        # Single argument
        if len(args) is 1:
            # If single string, could be a filepath, URL, date, or TimeRange
            if sunpy.time.is_time(args[0]):
                date = sunpy.time.parse_time(args[0])
                url = self._get_url_for_date(date)
                err = "Unable to download data for specified date"
                filepath = self._download(url, err)
            elif isinstance(args[0], basestring):
                # Filepath
                if os.path.isfile(os.path.expanduser(args[0])):
                    filepath = os.path.expanduser(args[0])
                else:
                    # Otherwise, assume string input is a URL
                    try:
                        filepath = self._download(args[0])
                    except (urllib2.HTTPError, urllib2.URLError, ValueError):
                        err = ("Unable to read location. Did you "
                               "specify a valid filepath or URL?")
                        raise ValueError(err)
            elif isinstance(args[0], sunpy.time.TimeRange):
                # TimeRange
                url = self._get_url_for_date_range(args[0])
                err = "Unable to download data for specified date range"
                filepath = self._download(url, err)   

            # Parse resulting file
            header, data = self._parse_filepath(filepath)
            
        # Date range
        if (sunpy.time.is_time(args[0]) and sunpy.time.is_time(args[1])):
            url = self._get_url_for_date_range(args[0], args[1])
            err = "Unable to download data for specified date range"
            filepath = self._download(url, err)
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
                
            data = pandas.DataFrame(args[0], index=index)
        else:
            raise TypeError("Unrecognized input for argument 1")
                
        # @NOTE: should we also support inputting start and end dates or a
        # date range?

        # Check for header
        header = ""
        
        if len(args) > 1:
            if (isinstance(args[1], basestring) or isinstance(args[1], dict)):
                header = args[1]
            else:
                raise TypeError("Unrecognized input for argument 2")
        
        self.data = data
        self.header = header
        
    def show(self, **kwargs):
        """Shows a plot of the light curve"""
        self.data.plot(**kwargs)
        plt.show()
        
    def _download(self, uri, err='Unable to download data at specified URL'):
        """Attempts to download data at the specified URI"""
        self._filename = os.path.basename(uri).split("?")[0]
        
        download_dir = sunpy.config.get("downloads", "download_dir")
        
        try:
            response = urllib2.urlopen(uri)
        except (urllib2.HTTPError, urllib2.URLError):
            raise urllib2.URLError(err)
            
        filepath = os.path.join(download_dir, self._filename)
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
        pass
    
    def _parse_fits(self, filepath):
        pass
    
    def _parse_filepath(self, filepath):
        filename, extension = os.path.splitext(filepath)
        
        if extension.lower() in (".csv", ".txt"):
            return self._parse_csv(filepath)
        else:
            return self._parse_fits(filepath)
