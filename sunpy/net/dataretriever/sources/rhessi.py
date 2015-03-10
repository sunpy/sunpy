#Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding provided by
#Google Summer of Code 2014

from sunpy.instr import rhessi

from ..client import GenericClient

__all__ = ['RHESSIClient']

class RHESSIClient(GenericClient):

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns a URL to the RHESSI data for the specified date range.
        Parameters
        ----------
        args : TimeRange, datetimes, date strings
        Date range should be specified using a TimeRange, or start
        and end dates at datetime instances or date strings.
        """
        if not timerange:
            return []

        url = rhessi.get_obssum_filename(timerange)
        return [url]


    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = ''
        self.map_['instrument'] = 'rhessi'
        self.map_['phyobs'] = 'irradiance'
        self.map_['provider'] = 'nasa'

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
        chkattr =  ['Time', 'Instrument']
        chklist =  [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'rhessi':
                return all(chklist)
        return False
