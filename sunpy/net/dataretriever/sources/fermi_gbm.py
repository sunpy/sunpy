from urllib.parse import urljoin

from ..client import GenericClient

__all__ = ['GBMClient']


class GBMClient(GenericClient):
    
    def _get_url_for_timerange(self, timerange, **kwargs):
        """

        Returns the url for Fermi/GBM data for the given date.

        baseurl = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'


        Parameters
        ----------
        timerange : sunpy.time.TimeRange
            time range for which to download the data

        Returns
        -------

        url for time of interest    

        """

            
        #checks if detector keyword
        #if not defaults to detector 5
        
        if 'detector' in kwargs:
            det = _check_detector(kwargs['detector'])

        else:
            det = 'n5'

        #check for datatype keyword - either CSPEC or CTIME
        #deault type is CSPEC

        if 'datatype' in kwargs:
            data_type = _check_type(kwargs['datatype'])
            
        else:
            data_type = 'cspec'

        #c = datatype
        days = timerange.get_dates()
        urls = []
        for day in days:
            urls.append(self._get_url_for_date(day, det, data_type, **kwargs))


        return urls

    def _get_url_for_date(self, date, det, data_type,  **kwargs):
        """
        Returns URL dats of interest

        Parameters
        ----------
        date : 'datetime.date'

        Returns
        -------
        url : str
            the url at the date of interest

        """
        
        filename = 'glg_'+data_type+'_' + det + date.strftime('_%y%m%d_') + 'v00.pha'
        url_path = urljoin(date.strftime('%Y/%m/%d/') + 'current/', filename)
        base_url = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'

        return urljoin(base_url, url_path)


    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'FERMI'
        self.map_['instrument'] = 'GBM'
        self.map_['physobs'] = 'flux'
        self.map_['provider'] = 'NASA'



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
        
        chkattr =  ['Time', 'Instrument', 'Detector', 'Datatype']
        chklist =  [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'gbm':
                return all(chklist)
        return False




def _check_detector(detector, **kwargs):
    """
    checks to see if detector is in right format
    """
    detector_numbers = [str(i) for i in range(12)]
    detector_list = ['n' + i for i in detector_numbers]

    if detector in detector_list:
        return detector
    elif detector in detector_numbers:
        return 'n' + detector
    else:
        raise ValueError('Detector number needs to be a string. Available detectors are n0-n11')


def _check_type(datatype, **kwargs):
    """
    checks is datatype is either CSPEC or CTIME
    """

    if not isinstance(datatype, str):
        raise ValueError('{} is not str - either cspec or ctime'.format(datatype))

    if datatype.lower() != 'cspec' and datatype.lower() != 'ctime':
        raise ValueError('{} not value datatype - either cspec or ctime'.format(datatype))

    else:
        return datatype.lower()
