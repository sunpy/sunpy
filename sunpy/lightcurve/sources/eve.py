# -*- coding: utf-8 -*-
"""Provides programs to process and analyze EVE data."""
from __future__ import absolute_import

import os
from datetime import datetime  

import matplotlib.pyplot as plt
from pandas.io.parsers import read_csv
from os.path import basename

from sunpy.lightcurve import LightCurve

__all__ = ['EVELightCurve']

class EVELightCurve(LightCurve):
    """
    SDO EVE LightCurve.
    
    Examples
    --------
    >>> import sunpy
    
    >>> eve = sunpy.lightcurve.EVELightCurve.create()
    >>> eve = sunpy.lightcurve.EVELightCurve.create('2012/06/20')
    >>> eve = sunpy.lightcurve.EVELightCurve.create(sunpy.data.test.EVE_AVERAGES_CSV)
    >>> eve = sunpy.lightcurve.EVELightCurve.create("http://lasp.colorado.edu/eve/data_access/quicklook/quicklook_data/L0CS/LATEST_EVE_L0CS_DIODES_1m.txt")
    >>> eve.peek(subplots=True)
    
    References
    ----------
    | http://lasp.colorado.edu/home/eve/data/data-access/
    """

    def peek(self, column = None, **kwargs):
        figure = plt.figure()
        # Choose title if none was specified
        if not kwargs.has_key("title") and column is None:
            if len(self.data.columns) > 1:
                kwargs['title'] = 'EVE (1 minute data)'
            else:
                if self._filename is not None:
                    base = self._filename.replace('_', ' ')
                    kwargs['title'] = os.path.splitext(base)[0]
                else:
                    kwargs['title'] = 'EVE Averages'

        if column is None:
            self.plot(**kwargs)
        else:
            data = self.data[column]
            if not kwargs.has_key("title"):
                kwargs['title'] = 'EVE ' + column.replace('_', ' ')
            data.plot(**kwargs)
        figure.show()
        return figure
        
    @staticmethod
    def _get_default_uri():
        """Load latest level 0CS if no other data is specified"""
        return "http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/LATEST_EVE_L0CS_DIODES_1m.txt"
    
    @staticmethod
    def _get_url_for_date(date):
        """Returns a URL to the EVE data for the specified date
        
            @NOTE: currently only supports downloading level 0 data
            .TODO: No data available prior to 2010/03/01!
        """
        base_url = 'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/'
        return base_url + date.strftime('%Y/%Y%m%d') + '_EVE_L0CS_DIODES_1m.txt'
    
    @classmethod
    def _parse_csv(cls, filepath):
        """Parses an EVE CSV file"""
        cls._filename = basename(filepath)
        with open(filepath, 'rb') as fp:
            # Determine type of EVE CSV file and parse
            line1 = fp.readline()
            fp.seek(0)

            if line1.startswith("Date"):
                return cls._parse_average_csv(fp)
            elif line1.startswith(";"):
                return cls._parse_level_0cs(fp)
    
    @staticmethod
    def _parse_average_csv(fp):
        """Parses an EVE Averages file"""
        return "", read_csv(fp, sep=",", index_col=0, parse_dates=True)
    
    @staticmethod
    def _parse_level_0cs(fp):
        """Parses and EVE Level 0CS file"""
        header = []
        fields = []
        line = fp.readline()
        header.append(line)
        # Read header at top of file
        while line.startswith(";"):
            header.append(line)
            line = fp.readline()

        fieldnames_start = False
        for l in header:
            if l.startswith("; Format:"):
                fieldnames_start = False
            if fieldnames_start:
                fields.append(l.split(":")[0].replace(';', ' ').strip())        
            if l.startswith("; Column descriptions:"):
                fieldnames_start = True

        # Next line is YYYY DOY MM DD        
        date_parts = line.split(" ")
                
        year = int(date_parts[0])
        month = int(date_parts[2])
        day = int(date_parts[3])
        #last_pos = fp.tell()
        #line = fp.readline()
        #el = line.split()
        #len
        
        # function to parse date column (HHMM)
        parser = lambda x: datetime(year, month, day, int(x[0:2]), int(x[2:4]))

        data = read_csv(fp, sep="\s*", names=fields, index_col=0, date_parser=parser, header = None)
        #data.columns = fields
        return header, data
