#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import urllib2


from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper

from sunpy.time import TimeRange

__all__ = ['KanzelhoheClient']

class KanzelhoheClient(GenericClient):
    """
    Returns a list of URLS to Kanzelhohe H-alpha files corresponding to value of input timerange.
    URL source: `http://cesar.kso.ac.at/halpha2k/recent/`.

    The earliest date available is from 20-July-2000

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'kanzelhohe'

    Returns
    -------
    urls: list
    list of urls corresponding to requested time range.

    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a

    >>> results = Fido.search(a.Time('2015/12/28 00:00:00', '2015/12/28 00:03:00'), a.Instrument('kanzelhohe'))
    >>> print(results)
    >>> [<Table length=1>
        Start Time           End Time              Source        Instrument
        str19               str19                str21           str10   
    ------------------- ------------------- --------------------- ----------
    2015-12-28 00:00:00 2015-12-29 00:00:00 Global Halpha Network Kanzelhohe]
    
    >>> response = Fido.fetch(results)
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        START_DATE = datetime.datetime(2000, 7, 20, 7, 45, 46)
        if timerange.start < START_DATE:
            raise ValueError('Earliest date for which SWAP data is available is '+ str(START_DATE))
        prefix = "http://cesar.kso.ac.at/halpha2k/recent/"
        suffix = "%Y/kanz_halph_fr_%Y%m%d_%H%M%S.fts.gz"
        url_pattern = prefix + suffix
        crawler = Scraper(url_pattern)
        if not timerange:
            return []
        result = crawler.filelist(timerange)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Global Halpha Network'
        self.map_['instrument'] = 'Kanzelhohe'
##        self.map_['phyobs'] = 'irradiance'
##        self.map_['provider'] = 'esa'
        self.map_['wavelength'] = '6563 AA'

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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if (x.__class__.__name__ == 'Instrument' and x.value.lower() == 'kanzelhohe'):
                return all(chklist)
        return False
            
        
    
