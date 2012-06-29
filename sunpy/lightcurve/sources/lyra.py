# -*- coding: utf-8 -*-
"""Provides programs to process and analyze PROBA2/LYRA data."""
from __future__ import absolute_import

from sunpy.lightcurve import LightCurve 

import os
from pandas.io.parsers import read_csv
from datetime import datetime  
from matplotlib import pyplot as plt
import urlparse

class LYRALightCurve(LightCurve):
    """SDO EVE light curve definition
    
    Examples
    --------
    >>> import sunpy
    >>> lyra = sunpy.lightcurve.LYRALightCurve()
    >>> lyra = sunpy.lightcurve.LYRALightCurve('~/Downloads/EVE_Fe_IX_171_averages.csv')
    >>> lyra = sunpy.lightcurve.LYRALightCurve('2012/06/20')
    >>> lya = sunpy.lightcurve.LYRALightCurve("http://lasp.colorado.edu/eve/data_access/quicklook/quicklook_data/L0CS/LATEST_EVE_L0CS_DIODES_1m.txt")
    >>> 
    >>> lyra.show()
    
    References
    ----------
    | http://lasp.colorado.edu/home/eve/data/data-access/
    """
    def __init__(self, *args, **kwargs):
        LightCurve.__init__(self, *args, **kwargs)

        
    def show(self, **kwargs):
        """Plots the LYRA data
        
        See: http://pandas.sourceforge.net/visualization.html
        """
        
        # Choose title if none was specified
        if not kwargs.has_key("title"):
            if len(self.data.columns) > 1:
                kwargs['title'] = 'LYRA data'
            else:
                if self._filename is not None:
                    base = self._filename
                    kwargs['title'] = os.path.splitext(base)[0]
                else:
                    kwargs['title'] = 'LYRA data'

        """Shows a plot of the light curve"""
        axes = self.data.plot(subplots=True, sharex=True)       
        plt.legend(loc='best')
        
        for i, name in enumerate(self.data.columns):
            axes[i].set_ylabel("%s (%s" % (name, "UNITS"))
            
        axes[0].set_title("LYRA")
        axes[-1].set_xlabel("Time")
        self.data.plot(**kwargs)
        plt.show()
    
    def _get_url_for_date(self, date):
        """Returns a URL to the LYRA data for the specified date
        """
        dt = sunpy.time.parse_time(date or datetime.datetime.utcnow())

        # Filename
        filename = "lyra_%s000000_lev%d_%s.fits" % (dt.strftime('%Y%m%d-'),
                                                    2, 'std')
        # URL
        base_url = "http://proba2.oma.be/lyra/data/bsd/"
        url_path = urlparse.urljoin(dt.strftime('%Y/%m/%d/'), filename)
        return urlparse.urljoin(base_url, url_path)
        
    def _get_default_uri(self):
        """Look for and download today's LYRA data"""
        return _get_url_for_date(self,datetime.utcnow())

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
