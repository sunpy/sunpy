#Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#This module was developed under funding provided by 
#Google Summer of Code 2014


from sunpy.net.unifieddownloader.client import GenericClient

__all__ = ['NOAAIndicesClient', 'NOAAPredictClient']

class NOAAIndicesClient(GenericClient):
    def _get_url_for_timerange(cls, timerange, **kwargs):
        """
        Helper function:
        """
        return ["http://www.swpc.noaa.gov/ftpdir/weekly/RecentIndices.txt"]

    def _makeimap(self):
        '''Helper Function:used to hold information about source. '''
        self.map_['source'] = 'sdic'
        self.map_['instrument'] = 'noaa'
        self.map_['phyobs'] = 'sunspot number'
        self.map_['provider'] = 'swpc'

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Boolean Function:Answers whether client can service the query.
        """
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'noaa-indices':
                return all(chklist)
        return False


class NOAAPredictClient(GenericClient):

    def _get_url_for_timerange(cls, timerange, **kwargs):
        """
        Helper function:
        """
        return ["http://www.swpc.noaa.gov/ftpdir/weekly/Predict.txt"]

    def _makeimap(self):
        '''Helper Function:used to hold information about source. '''
        self.map_['source'] = 'ises'
        self.map_['instrument'] = 'noaa'
        self.map_['phyobs'] = 'sunspot number'
        self.map_['provider'] = 'swpc'

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Boolean Function:Answers whether client can service the query.
        """
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value == 'noaa-predict':
                return all(chklist)
        return False
