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
        An pandas DataFrame prepresenting one or more fields as they vary with
        respect to time.
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
    | http://pandas.pydata.org/pandas-docs/dev/dsintro.html#dataframe
    
    """
    def __init__(self, *args, **kwargs):
        self._filename = ""

        # If no arguments specified, perform default action
        if len(args) is 0:
            args = (self._get_default_uri(),)

        # If single string, could be a filepath, URL, or date string
        if len(args) is 1 and isinstance(args[0], basestring):
            # Filepath
            if os.path.isfile(os.path.expanduser(args[0])):
                filepath = os.path.expanduser(args[0])
            else:
                # Check to see if string is a valid date
                try:
                    date = sunpy.time.parse_time(args[0])
                    url = self._get_url_for_date(date)
                    try:
                        filepath = self._download(url)
                    except (urllib2.HTTPError, urllib2.URLError, ValueError):
                        err = "Unable to download data for specified date"
                        raise urllib2.URLError(err)

                except ValueError:
                    # Otherwise, assume string input is a URL
                    try:
                        filepath = self._download(args[0])
                    except (urllib2.HTTPError, urllib2.URLError, ValueError):
                        raise Exception("Unable to read location. Did you "
                                        "specify a valid filepath or URL?")

            # Parse file
            filename, extension = os.path.splitext(filepath)
            
            if extension.lower() in ("csv", "txt"):
                header, data = self._parse_csv(filepath)
            else:
                header, data = self._parse_fits(filepath)
                
        # @NOTE: should we also support inputting start and end dates or a
        # date range?

        # Other light curve creation options (DataFrame, ndarray, etc)
        header = ""
        
        # DataFrame Index
        if "index" in kwargs:
            index = kwargs["index"]
        else:
            index = None
            
        # DataFrame input
        if isinstance(args[0], pandas.DataFrame):
            data = args[0]
        elif (isinstance(args[0], list) or 
              isinstance(args[0], dict) or 
              isinstance(args[0], np.ndarray) or 
              isinstance(args[0], pandas.Series)):
            # list, dict, ndarray, or Series
            data = pandas.DataFrame(args[0], index=index)
        else:
            raise TypeError("Unrecognized input for argument 1")
        
        # Check for header
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
        
    def _download(self, uri):
        """Attempts to download data at the specified URI"""
        self._filename = os.path.basename(uri)
        
        response = urllib2.urlopen(uri)
        filepath = os.path.join(sunpy.config['data.directory'], self._filename)
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
    
#    @classmethod
#    def parse_file(cls, filepath):
#        """Reads in a map file and returns a header and data array"""
#        data, dict_header = read_file(filepath)
#
#        return dict_header, data
#
#    @classmethod
#    def read(cls, filepath):
#        """LightCurve class factory
#
#        Attempts to determine the type of data associated with input and
#        returns a LightCurve subclass instance. If the file format is not
#        recognized a warning will be displayed.
#
#        Parameters
#        ----------
#        filepath : string
#            Path to input file (FITS, CSV, etc.)
#
#        Returns
#        -------
#        out : LightCurve
#            Returns a LightCurve instance.
#        """
#        header, data = cls.parse_file(filepath)
#
#        if cls.__name__ is not "LightCurve":
#            return cls(filepath)
#
#        for cls in LightCurve.__subclasses__():
#            if cls.is_datasource_for(header):
#                return cls(data, header)


