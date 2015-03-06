#Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding by
#Google Summer of Code 2014

import urlparse

from sunpy.net.unifieddownloader.client import GenericClient


__all__ = ['EVEClient']


class EVEClient(GenericClient):

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns list of URLS corresponding to value of input timerange.

        Parameters
        ----------
        timerange: sunpy.time.TimeRange
            time range for which data is to be downloaded.

        Returns
        -------
        urls : list
            list of URLs corresponding to the requested time range
        """
        if not timerange:
            return []
        days = timerange.get_dates()
        urls = []
        for day in days:
            urls.append(self._get_url_for_date(day, **kwargs))
        return urls

    def _get_url_for_date(self, date, **kwargs):
        """
        Return URL for corresponding date.

        Parameters
        ----------
        date : Python datetime object

        Returns
        -------
        URL : string
        """
        base_url = 'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/'
        return urlparse.urljoin(base_url, date.strftime('%Y/%Y%m%d') + '_EVE_L0CS_DIODES_1m.txt')

    def _makeimap(self):
        """
        Helper Function: used to hold information about source.
        """
        self.map_['source'] = 'SDO'
        self.map_['provider'] ='LASP'
        self.map_['instrument'] = 'eve'
        self.map_['phyobs'] = 'irradiance'

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
        chkattr =  ['Time', 'Instrument','Level']
        chklist =  [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'eve':
                chk_var +=1
            elif x.__class__.__name__ == 'Level' and x.value == 0:
                chk_var +=1

        if(chk_var == 2):
            return all(chklist)
        return False
