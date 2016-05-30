#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import urllib2
import re
from bs4 import BeautifulSoup

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange


_all__ = ['VSMClient']

class VSMClient(GenericClient):

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        wave = int(kwargs['Wavelength'].min.value)
        table = {6302:'v72', 8542:'v82', 10830: 'v22'}
        if (wave == 6302):
            START_DATE = datetime.datetime(2003, 8, 21)
        elif (wave == 8542):
            START_DATE = datetime.datetime(2003, 8, 26)
        elif (wave == 10830):
            START_DATE = datetime.datetime(2004, 11, 4)

        if timerange.start < START_DATE:
            raise ValueError('Earliest date for which Kanzelhohe data is available is {:%Y-%m-%d}'.format(START_DATE))

        prefix = 'http://gong2.nso.edu/pubkeep/{wave_type}/%Y%m'
        suffix = 'k4{wave_type}%y%m%dt%H%M%S.fts.gz'
        total_days = (timerange.end - timerange.start).days + 1
        all_days = timerange.split(total_days)

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'SOLIS'
        self.map_['instrument'] = 'vsm'

    @classmethod
    def _can_handle_query(cls, *query):
        
        """
        Answers whether client can service the query.
        
        Parameters
        ----------
        query : list of query objects
        
        Returns
        -------
        boolean: answer as to whether client can service the query
        
        """
        chkattr = ['Time', 'Instrument', 'Wavelength']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        values = [6302, 8542, 10830]
        for x in query:
            if (x.__class__.__name__ == 'Instrument' and type(x.value) is str and x.value.lower() == 'vsm'):
                chk_var += 1
            if (x.__class__.__name__ == 'Wavelength' and int(x.min.value) in values and int(x.max.value) in values and (x.unit.name).lower()=='angstrom'):
                chk_var += 1
        if (chk_var==2):
            return True
        return False
        
