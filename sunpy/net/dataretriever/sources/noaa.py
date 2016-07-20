#Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding provided by
#Google Summer of Code 2014


from ..client import GenericClient
import datetime

__all__ = ['NOAAIndicesClient', 'NOAAPredictClient', 'SRSClient']

class NOAAIndicesClient(GenericClient):

    @staticmethod
    def _get_default_uri():
        """Return the url to download indices"""
        return ["ftp://ftp.swpc.noaa.gov/pub/weekly/RecentIndices.txt"]

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Helper function:
        """
        return NOAAIndicesClient._get_default_uri()

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'sdic'
        self.map_['instrument'] = 'noaa-indices'
        self.map_['phyobs'] = 'sunspot number'
        self.map_['provider'] = 'swpc'

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
            if x.__class__.__name__ == 'Instrument' and x.value == 'noaa-indices':
                return all(chklist)
        return False


class NOAAPredictClient(GenericClient):

    @staticmethod
    def _get_default_uri():
        """Return the url to download indices"""
        return ["http://services.swpc.noaa.gov/text/predicted-sunspot-radio-flux.txt"]

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Helper function:
        """
        return NOAAPredictClient._get_default_uri()

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ises'
        self.map_['instrument'] = 'noaa-predict'
        self.map_['phyobs'] = 'sunspot number'
        self.map_['provider'] = 'swpc'

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
            if x.__class__.__name__ == 'Instrument' and x.value == 'noaa-predict':
                return all(chklist)
        return False

class SRSClient(GenericClient):

    @staticmethod
    def _get_default_uri():
        today = datetime.datetime.utcnow()
        return ['ftp://ftp.swpc.noaa.gov/pub/warehouse/{date:%Y}/SRS/{date:%Y%m%d}SRS.txt'.format(date=today)]

    def _get_url_for_timerange(self, timerange, **kwargs):

        if not timerange:
            return SRSClient._get_default_uri()
        result = list()
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/warehouse/'
        total_days = (timerange.end - timerange.start).days + 1
        all_dates = timerange.split(total_days)
        for day in all_dates:
            url = base_url + '%Y/SRS/%Y%m%dSRS.txt'.format(day)
            result.append(url)
        return result

    def _makeimap(self):
        self.map_['source'] = 'swpc'
        self.map_['instrument'] = 'SOON'
        self.map_['physobs'] = 'SRS'
        self.map_['source'] = 'NOAA/USAF'

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
        chkattr = ["Time", "Instrument"]
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == "Instrument" and x.value == "SOON":
                return True
        return False
        
