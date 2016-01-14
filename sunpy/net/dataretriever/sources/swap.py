''' SWAPClient definition '''

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import urllib2


from astropy import units as u
from bs4 import BeautifulSoup
from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper

from sunpy.time import TimeRange

__all__ = ['SWAPClient']

class SWAPClient(GenericClient):
    
    """ Returns a list of URLS to Proba2 SWAP files corresponding to value of input timerange.
    URL source: http://proba2.oma.be/swap/data/bsd/

    The earliest date available is from 24-Nov-2009

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'swap'

    Returns
    -------
    urls: list
        list of urls corresponding to requested time range.

    Examples
    --------
    >>> from sunpy.time.timerange import TimeRange
    >>> from sunpy.net.vso.attrs import Time,Instrument,Source,Level
    >>> from sunpy.net.dataretriever.client import QueryResponse
    >>> from sunpy.net.dataretriever.sources import swap

    >>> LCClient = swap.SWAPClient()

    >>> qr = LCClient.query(Time('2015-12-30 00:00:00','2015-12-31 00:05:00'),Instrument('swap'))
    >>> res = LCClient.get(qr)
    >>> dl = res.wait()
    """
    
    def _get_url_for_timerange(cls, timerange, **kwargs):
        """ returns list of urls corresponding
        to given TimeRange. """

        if timerange.start < datetime.datetime(2009,11,24):
            raise ValueError('Earliest date for which SWAP data is available is 2009-11-24')
        
        url_pattern = url_pattern = ('http://proba2.oma.be/swap/data/bsd/'
               '%Y/%m/%d/'
               '{instrument}_lv1_%Y%m%d_%H%M%S.fits')
        
        crawler = Scraper(url_pattern, instrument = 'swap')
        if not timerange:
            return []
        result = crawler.filelist(timerange)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Proba2'
        self.map_['instrument'] = 'swap'
        self.map_['phyobs'] = 'irradiance'
        self.map_['provider'] = 'esa'

    @classmethod
    def _can_handle_query(cls, *query):
         """
        Answers whether client can service the query.
        
        Parameters
        ----------
        query : list of query objects
        
        Returns
        -------
        boolean
            answer as to whether client can service the query
        
        """
         
         chkattr = ['Time','Instrument','Level']
         chklist = [x.__class__.__name__ in chkattr for x in query]
         for x in query:
             if x.__class__.__name__ == 'Instrument' and x.value=='lyra':
                 return all(chklist)
         return False

    
        
        
        
        
    
