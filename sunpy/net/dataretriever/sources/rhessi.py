# Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
# This module was developed under funding provided by
# Google Summer of Code 2014

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
        return rhessi.get_obssum_filename(timerange)

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'rhessi'
        self.map_['instrument'] = 'rhessi'
        self.map_['physobs'] = 'irradiance'
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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'rhessi':
                return all(chklist)
        return False
