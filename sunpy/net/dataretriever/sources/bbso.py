#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import urllib2
from bs4 import BeautifulSoup
from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ['BBSOClient']

class BBSOClient(GenericClient):
    """
    Returns a list of URLS to BBSO files corresponding to value of input timerange.
    URL source: `http://www.bbso.njit.edu/pub/archive/`.


    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'bbso'

    Level: Level can take only 'fl' or 'fr' as arguments.

    Returns
    -------
    urls: list
    list of urls corresponding to requested time range.

    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a

    >>> results = Fido.search(Time('2016/5/18 15:28:00','2016/5/18 16:30:00'), a.Instrument('bbso'), a.Level('fr'))
    >>> print(results)
    [<Table length=2>
         Start Time           End Time              Source        Instrument
           str19               str19                str21            str4   
    ------------------- ------------------- --------------------- ----------
    2016-05-18 00:00:00 2016-05-19 00:00:00 Global Halpha Network       bbso
    2016-05-19 00:00:00 2016-05-20 00:00:00 Global Halpha Network       bbso]
    
    >>> response = Fido.fetch(results)
    """
    
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        level = kwargs.get('level', 'fr')
        START_DATE = datetime.datetime(2000, 11, 6)
        if timerange.start < START_DATE:
            raise ValueError('Earliest date for which BBSO data is available is {:%Y-%m-%d}'.format(START_DATE))
        prefix = 'http://www.bbso.njit.edu'
        suffix = '/pub/archive/%Y/%m/%d/bbso_halph_{level}_%Y%m%d_%H%M%S.fts'
        suffix_gz = suffix + '.gz' #Downlaod compressed files as well, if available.
        total_days = (timerange.end - timerange.start).days + 1
        all_dates = timerange.split(total_days)
        crawler = Scraper(suffix, level = level)
        crawler_gz = Scraper(suffix_gz, level = level) 
        result = list()
        for day in all_dates:
            url_to_open = 'http://www.bbso.njit.edu/pub/archive/{date:%Y/%m/%d/}'.format(date = day.end)
            html_page = urllib2.urlopen(url_to_open)
            soup = BeautifulSoup(html_page)
            for link in soup.findAll('a'):
                url = str(link.get('href'))
                if crawler._URL_followsPattern(url):
                    datehref = crawler._extractDateURL(url)
                    if datehref>=timerange.start and datehref<=timerange.end:
                        result.append(prefix + url)
                if crawler_gz._URL_followsPattern(url):
                    datehref = crawler_gz._extractDateURL(url)
                    if datehref>=timerange.start and datehref<=timerange.end:
                        result.append(prefix + url)           
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Global Halpha Network'
        self.map_['instrument'] = 'bbso'
        self.map_['phyobs'] = 'irradiance'
##        self.map_['provider'] = 'esa'
        self.map_['wavelength'] = '174 AA'

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
        chkattr = ['Time', 'Instrument', 'Level']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'bbso':
                chk_var += 1
            if x.__class__.__name__ == 'Level' and type(x.value) is str and x.value.lower() in ('fl', 'fr'):
                chk_var += 1
        if (chk_var == 2):
            return True
        return False
