#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import urllib2
import re


from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ['XRTClient']

class XRTClient(GenericClient):
    """
    Returns a list of URLS to XRT files corresponding to value of input timerange.
    URL source: `http://solar.physics.montana.edu/HINODE/XRT/QL/syn_comp_fits/`.
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'xrt'
    Filter: Filter can take only 'almesh', 'alpoly', 'cpoly', 'tipoly', 'thinbe' as arguments.
            Arguments are case-insensitive and can include non alphabetic characters as well.
            Arguments like 'Al-mesh' or 'Al_mesh' or 'Al/mesh' would also work fine.
            
    Returns
    -------
    urls: list
    list of urls corresponding to requested time range.
    
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> results = Fido.search(Time('2016/5/18','2016/5/19'), a.Instrument('xrt'), a.Filter('al_mesh'))
    >>> print(results)
        [<Table length=2>
         Start Time           End Time      Source Instrument
           str19               str19         str6     str3   
    ------------------- ------------------- ------ ----------
    2016-05-18 00:00:00 2016-05-19 00:00:00 Hinode        XRT
    2016-05-19 00:00:00 2016-05-20 00:00:00 Hinode        XRT]
    
    >>> response = Fido.fetch(results)
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        START_DATE = datetime.datetime(2014, 1, 10, 18, 14, 42)
        filter_type = kwargs['filter']
        regex = re.compile('[^a-zA-Z]')
        filter_type = regex.sub('', filter_type).lower()
        result = list()
        if timerange.start < START_DATE:
            raise ValueError('Earliest date for which XRT data is available is {:%Y-%m-%d}'.format(START_DATE))
        filter_dict = {'almesh':'Al_mesh', 'alpoly':'Al_poly', 'cpoly':'C_poly', 'tipoly':'Ti_poly', 'thinbe':'thin_Be'}
        url_pattern = ['http://solar.physics.montana.edu/HINODE/XRT/QL/syn_comp_fits/XRT_{filter}_%Y%m%d_%H%M%S.{s}.fits'.
                       format(filter = filter_dict[filter_type], s=i) for i in range (0, 10)]
        arr = [Scraper(pattern, filter = filter_dict[filter_type]).filelist(timerange) for pattern in url_pattern]
        [result.extend(url) for url in arr if len(url)>0]
        #Integrate nascom.nasa.gov XRT data as well.
        url_pattern = ['http://sohowww.nascom.nasa.gov/sdb/hinode/xrt/l1q_synop/XRT_{filter}_%Y%m%d_%H%M%S.{s}.fits'.
                       format(filter = filter_dict[filter_type], s=i) for i in range (0, 10)]
        arr = [Scraper(pattern, filter = filter_dict[filter_type]).filelist(timerange) for pattern in url_pattern]
        [result.extend(url) for url in arr if len(url)>0]
        if not timerange:
            return []
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Hinode'
        self.map_['instrument'] = 'XRT'

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
        chkattr = ['Time', 'Instrument', 'Filter']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        types = ['almesh', 'alpoly', 'cpoly', 'tipoly', 'thinbe']
        regex = re.compile('[^a-zA-Z]')
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'xrt':
                chk_var += 1
            if x.__class__.__name__ == 'Filter' and type(x.value) is str and regex.sub('',x.value).lower() in types:
                chk_var += 1
        if (chk_var == 2):
            return True
        return False
