"""
LightCurve is a generic LightCurve class from which all other LightCurve classes 
inherit from.
"""
from __future__ import absolute_import

#pylint: disable=E1101,E1121,W0404,W0613
__authors__ = ["Keith Hughitt"]
__email__ = "keith.hughitt@nasa.gov"

import os
import pandas
import sunpy
import urllib2
from warnings import warn
from sunpy.io import read_file

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

    References
    ----------
    | http://pandas.pydata.org/pandas-docs/dev/dsintro.html#dataframe
    
    """
    def __init__(self, *args, **kwargs):
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
                        warn("Unable to download data for specified date")
                        return
                except ValueError:
                    # Otherwise, assume string input is a URL
                    try:
                        filepath = self._download(args[0])
                    except (urllib2.HTTPError, urllib2.URLError, ValueError):
                        warn("Unable to read location. Did you specify "
                                      "a valid filepath or URL?")
                        return

            # Parse file
            header, data = self._parse_csv(filepath)


        # Start and end dates?
        # Date range?
        
        # DataFrame
        
        # Python list
        
        # ndarray
        
        self.data = data
        self.header = header
        
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
        warn("No default action set for %s" % self.__class__.__name__)
        
    def _get_url_for_date(self, date):
        """Returns a URL to the data for the specified date"""
        warn("Date-based downloads not supported for for %s" % 
                      self.__class__.__name__)
        
    def parse_csv(self, filepath):
        """CSV parsing support"""
        warn("CSV support not yet implemented for %s" % 
                      self.__class__.__name__)
    
    @classmethod
    def parse_file(cls, filepath):
        """Reads in a map file and returns a header and data array"""
        data, dict_header = read_file(filepath)

        return dict_header, data

    @classmethod
    def read(cls, filepath):
        """LightCurve class factory

        Attempts to determine the type of data associated with input and
        returns a LightCurve subclass instance. If the file format is not
        recognized a warning will be displayed.

        Parameters
        ----------
        filepath : string
            Path to input file (FITS, CSV, etc.)

        Returns
        -------
        out : LightCurve
            Returns a LightCurve instance.
        """
        header, data = cls.parse_file(filepath)

        if cls.__name__ is not "LightCurve":
            return cls(filepath)

        for cls in LightCurve.__subclasses__():
            if cls.is_datasource_for(header):
                return cls(data, header)
        
        # DISPLAY WARNING..

