# -*- coding: utf-8 -*-
"""Provides programs to process and analyze EVE data."""
from __future__ import absolute_import

from sunpy.lightcurve import LightCurve

import os
import matplotlib.pyplot as plt
from pandas.io.parsers import read_csv
from datetime import datetime  
from matplotlib import pyplot as plt

class EVELightCurve(LightCurve):
    """SDO EVE light curve definition
    
    References
    ----------
    | http://lasp.colorado.edu/home/eve/data/data-access/
    """
    def __init__(self, *args):
        self._filepath = None

        # Filepath
        if len(args) is 1 and isinstance(args[0], basestring):
            self._filepath = args[0]
            header, data = self._parse_csv(args[0])

        # Start and end dates
            
        # Date range
        
        LightCurve.__init__(self, data, header)
        
    def show(self, **kwargs):
        # Choose title if none was specified
        if not kwargs.has_key("title"):
            if len(self.data.columns) > 1:
                kwargs['title'] = 'EVE GOES Proxy Xray Flux (1 minute data)'
            else:
                if self._filepath is not None:
                    base = os.path.basename(self._filepath).replace('_', ' ')
                    kwargs['title'] = os.path.splitext(base)[0]
                else:
                    kwargs['title'] = 'EVE Averages'

        self.data.plot(**kwargs)
        plt.show()
    
    def _parse_fits(self, filepath):
        """Parses an EVE FITS file"""
        pass
    
    def _parse_csv(self, filepath):
        """Parses an EVE CSV file"""
        fp = open(filepath, 'rb')
        
        # Determine type of EVE CSV file and parse
        line1 = fp.readline()
        fp.seek(0)

        if line1.startswith("Date"):
            return self._parse_average_csv(fp)
        elif line1.startswith(";"):
            return self._parse_level_0cs(fp)
    
    def _parse_average_csv(self, fp):
        """Parses an EVE Averages file"""
        return "", read_csv(fp, sep=",", index_col=0, parse_dates=True)
    
    def _parse_level_0cs(self, fp):
        """Parses and EVE Level 0CS file"""
        header = ""
        
        fields = ('xrs-b', 'xrs-a', 'sem', 'ESPquad', 'esp171', 
                  'esp257', 'esp304', 'esp366', 'espdark', 'megsp', 'megsdark', 
                  'q0esp', 'q1esp', 'q2esp', 'q3esp', 'cmlat', 'cmlon')
        
        line = fp.readline()
        
        # Read comment at top of file
        while line.startswith(";"):
            header += line
            line = fp.readline()
            
        # Next line is YYYY DOY MM DD
        parts = line.split(" ")
                
        year = int(parts[0])
        month = int(parts[2])
        day = int(parts[3])
        
        # function to parse date column (HHMM)
        parser = lambda x: datetime(year, month, day, int(x[0:2]), int(x[2:4]))

        data = read_csv(fp, sep="\s*", names=fields, index_col=0, date_parser=parser)
        
        return header, data
