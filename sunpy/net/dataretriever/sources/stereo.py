#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import re

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ['SECCHIClient', 'PLASTICClient', 'IMPACTClient', 'SWAVESClient']

class SECCHIClient(GenericClient):
    """
    Returns a list of URLS to SECCHI files corresponding to value of input timerange.
    URL source: `http://stereo-ssc.nascom.nasa.gov/data/beacon/`.
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'secchi'
    
    Source: Two arguments, case-insensitive - 'ahead' or 'behind'

    Detector: Can take any four arguments, case-insensitive - 'euvi', 'cor2', 'hi_1', 'hi_2'
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a    
    >>> results = Fido.search(a.Time('2007/3/2 17:00:00', '2007/3/3 17:15:00'), a.Instrument('secchi'), a.Source('ahead'), a.Detector('euvi'))
    >>> print(results)
        [<Table length=4>
        Start Time           End Time      Source Instrument
           str19               str19         str5     str6   
    ------------------- ------------------- ------ ----------
    2007-03-02 00:00:00 2007-03-03 00:00:00  ahead     secchi
    2007-03-03 00:00:00 2007-03-04 00:00:00  ahead     secchi
    2007-03-04 00:00:00 2007-03-05 00:00:00  ahead     secchi
    2007-03-05 00:00:00 2007-03-06 00:00:00  ahead     secchi]
    >>> response = Fido.fetch(results)
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        regex = re.compile('[^a-zA-Z0-9]')
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
        """
        Helper Function: used to hold information about source.
        """       
        self.map_['instrument'] = 'secchi'

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
        regex = re.compile('[^a-zA-Z0-9]')
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
            if x.__class__.__name__ == 'Detector' and regex.sub('',x.value).lower() in detectors:
                chk_var += 1
        if (chk_var == 3):
            return True
        return False
                

class PLASTICClient(GenericClient):
    """
    Returns a list of URLS to PLASTIC files corresponding to value of input timerange.
    URL source: `http://stereo-ssc.nascom.nasa.gov/data/beacon/`.
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument, case-insensitive = 'plastic'
    
    Source: Two arguments, case-insensitive - 'ahead' or 'behind'
    
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a    
    >>> results = Fido.search(a.Time('2007/3/2 17:00:00', '2007/3/4 17:15:00'), a.Instrument('plastic'), a.Source('ahead'))
    >>> print(results)
       [<Table length=2>
         Start Time           End Time      Source Instrument
           str19               str19         str5     str7   
        ------------------- ------------------- ------ ----------
        2007-03-02 00:00:00 2007-03-03 00:00:00  ahead    plastic
        2007-03-03 00:00:00 2007-03-04 00:00:00  ahead    plastic]
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        source = kwargs.get('source').lower()
        if (source == 'ahead'):
            START_DATE = datetime.datetime(2006, 10, 1)
        else:
            START_DATE = datetime.datetime(2006, 10, 12)
        if (timerange.start < START_DATE):
            raise ValueError('Earliest date for which PLASTIC data is available is {:%Y-%m-%d}'.format(START_DATE))
        prefix = 'http://stereo-ssc.nascom.nasa.gov/data/beacon/{source}/plastic/'
        if (source == 'ahead'):
            suffix = '%Y/%m/ST{char}_LB_PLASTIC_%Y%m%d_V06.cdf'
        else:
            suffix = '%y/%m/ST{char}_LB_PLASTIC_%Y%m%d_V06.cdf'
        url_pattern = prefix + suffix
        crawler = Scraper(url_pattern, source = source, char = source[0].upper())
        if not timerange:
            return []
        result = crawler.filelist(timerange)
        return result

    def _makeimap(self):
        """
        Helper Function: used to hold information about source.
        """       
        self.map_['instrument'] = 'plastic'

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
        chkattr = ['Time', 'Instrument', 'Source']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        sources = ['ahead', 'behind']
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'plastic':
                chk_var += 1
            if x.__class__.__name__ == 'Source' and x.value.lower() in sources:
                chk_var += 1
        if (chk_var == 2):
            return True
        return False

class IMPACTClient(GenericClient):
    """
    Returns a list of URLS to IMPACT files corresponding to value of input timerange.
    URL source: `http://stereo-ssc.nascom.nasa.gov/data/beacon/`.
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument, case-insensitive = 'impact'
    
    Source: Two arguments, case-insensitive - 'ahead' or 'behind'
    
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a    
    >>> results = Fido.search(a.Time('2007/3/20 17:00:00', '2007/3/25 17:15:00'), a.Instrument('impact'), a.Source('ahead'))
    [<Table length=5>
        Start Time           End Time      Source Instrument
           str19               str19         str5     str6   
    ------------------- ------------------- ------ ----------
    2007-03-20 00:00:00 2007-03-21 00:00:00  ahead     impact
    2007-03-21 00:00:00 2007-03-22 00:00:00  ahead     impact
    2007-03-22 00:00:00 2007-03-23 00:00:00  ahead     impact
    2007-03-23 00:00:00 2007-03-24 00:00:00  ahead     impact
    2007-03-24 00:00:00 2007-03-25 00:00:00  ahead     impact]
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        source = kwargs.get('source').lower()
        START_DATE = datetime.datetime(2006, 10, 1)
        if (timerange.start < START_DATE):
            raise ValueError('Earliest date for which IMPACT data is available is {:%Y-%m-%d}'.format(START_DATE))
        prefix = 'http://stereo-ssc.nascom.nasa.gov/data/beacon/{source}/'
        suffix = 'impact/%Y/%m/ST{char}_LB_IMPACT_%Y%m%d_V01.cdf'
        url_pattern = prefix + suffix
        crawler = Scraper(url_pattern, source = source, char = source[0].upper())
        if not timerange:
            return []
        result = crawler.filelist(timerange)
        return result
    def _makeimap(self):
        """
        Helper Function: used to hold information about source.
        """      
        self.map_['instrument'] = 'impact'

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
        chkattr = ['Time', 'Instrument', 'Source']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        sources = ['ahead', 'behind']
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'impact':
                chk_var += 1
            if x.__class__.__name__ == 'Source' and x.value.lower() in sources:
                chk_var += 1
        if (chk_var == 2):
            return True
        return False


class SWAVESClient(GenericClient):
    """
    Returns a list of URLS to SWAVES files corresponding to value of input timerange.
    URL source: `http://stereo-ssc.nascom.nasa.gov/data/beacon/`.
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument, case-insensitive = 'swaves'
    
    Source: Two arguments, case-insensitive - 'ahead' or 'behind'
    
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a    
    >>> results = Fido.search(a.Time('2008/3/20 17:00:00', '2008/3/25 17:15:00'), a.Instrument('swaves'), a.Source('ahead'))
    >>> print(results)
    [<Table length=5>
     Start Time           End Time      Source Instrument
       str19               str19         str5     str6   
    ------------------- ------------------- ------ ----------
    2008-03-20 00:00:00 2008-03-21 00:00:00  ahead     swaves
    2008-03-21 00:00:00 2008-03-22 00:00:00  ahead     swaves
    2008-03-22 00:00:00 2008-03-23 00:00:00  ahead     swaves
    2008-03-23 00:00:00 2008-03-24 00:00:00  ahead     swaves
    2008-03-24 00:00:00 2008-03-25 00:00:00  ahead     swaves]
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        source = kwargs.get('source')
        START_DATE = datetime.datetime(2006, 10, 27)
        if (timerange.start < START_DATE):
            raise ValueError('Earliest date for which SWAVES data is available is {:%Y-%m-%d}'.format(START_DATE))        
        prefix = 'http://stereo-ssc.nascom.nasa.gov/data/beacon/{source}/'
        suffix = 'swaves/%Y/%m/ST{char}_LB_SWAVES_%Y%m%d.idlsave'
        url_pattern = prefix + suffix
        crawler = Scraper(url_pattern, source = source, char = source[0].upper())
        if not timerange:
            return []
        result = crawler.filelist(timerange)
        return result

    def _makeimap(self):
        """
        Helper Function: used to hold information about source.
        """
        self.map_['instrument'] = 'swaves'

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
        chkattr = ['Time', 'Instrument', 'Source']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        sources = ['ahead', 'behind']
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'swaves':
                chk_var += 1
            if x.__class__.__name__ == 'Source' and x.value.lower() in sources:
                chk_var += 1
        if (chk_var == 2):
            return True
        return False
