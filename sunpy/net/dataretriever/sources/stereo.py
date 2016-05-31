#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import re

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ['SECCHIClient']

class SECCHIClient(GenericClient):

    def _get_url_for_timerange(self, timerange, kwargs):

        regex = re.compile('[^a-zA-Z]')
        source = kwargs.get('source', 'ahead').lower() 
        detector = regex.sub('', kwargs.get('detector')).lower() #default to euvi ?
        table = { 'euvi':'euvi', 'hi1' : 'hi_1', 'hi2': 'hi_2', 'cor2':'cor2'}
        prefix = 'http://stereo-ssc.nascom.nasa.gov/data/beacon/{source}/secchi/img/{det}/%Y%m%d/%Y%m%d_%H%M%S_'
        suffix_table = { 'euvi':'s7eu{char}.fts', 'cor2': 'd7c2{char}.fts', 'hi1':'sHh1{char}.fts', 'hi2':'sHh2B{char}.fts'}
        url_pattern = prefix + suffix_table[detector].format(char = source[0].upper())
        crawler = Scraper(url_pattern, source = source, det = detector)
        if not timerange:
            return []
        result = crawler.filelist(timerange)
        return result

    def _makeimap(self):
        self.map_['instrument'] = 'secchi'


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
        regex = re.compile('[^a-zA-Z]')
        chkattr = ['Time', 'Instrument', 'Source', 'Detector']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        sources = ['ahead', 'behind']
        detectors = ['euvi', 'cor2', 'hi1', 'hi2']
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'secchi':
                chk_var += 1
            if x.__class__.__name__ == 'Source' and x.value.lower() in sources:
                chk_var += 1
            if (x.__class__.__name__ == 'Detector' and regex.sub('',x.value).lower() in detectors):
                chk_var += 1
        if (chk_var == 3):
            return True
        return False
                
                
