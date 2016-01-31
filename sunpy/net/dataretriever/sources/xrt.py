''' XRTClient definition. '''

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

#Note to self: Search both sources for XRT data.

import datetime
import urllib2

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper

from sunpy.time import TimeRange

__all__ = ['XRTClient']

class XRTClient(GenericClient):

    def _get_url_for_timerange(cls, timerange, filter_arg):
        """ returns a list of urls corresponding
        to given timerange. """

        XRT_STARTDATE = datetime.datetime(2014, 01, 10)
        if timerange.start < XRT_STARTDATE:
            raise ValueError('Earliest date for which XRT data available is 2014-01-10')
        result = []
        for i in range (0,10):
            url_pattern = ('http://solar.physics.montana.edu/HINODE/XRT/QL/syn_comp_fits/'
               '{instrument}'
               '_{filter_arg}_%Y%m%d_%H%M%S.'+str(i)+'.fits')

            ans = Scraper(url_pattern, instrument='XRT', filter_arg = filter_arg)

            files = ans.filelist(timerange)
            for links in files:
                result.append(links)

        return result

    def _makeimap(self):
        self.map_['source'] = 'HINODE'
        self.map_['instrument'] = 'XRT'
        self.map_['phyobs'] = 'x-rays'
##        self.map_['provider'] =

    def _can_handle_query(cls, *query):
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'XRT':
                return all(chklist)
        return False
        
    
