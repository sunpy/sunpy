from urllib.parse import urljoin
from sunpy.util.scraper import Scraper
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

        # Checks if detector keyword
        # If not defaults to detector 5
        if 'detector' in kwargs:
            det = _check_detector(kwargs['detector'])

        else:
            det = 'n5'

        # Check for datatype keyword - either CSPEC or CTIME
        # Default type is CSPEC

        if 'datatype' in kwargs:
            data_type = _check_type(kwargs['datatype'])
        else:
            data_type = 'cspec'

        gbm_pattern = ('https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'
                       '%Y/%m/%d/current/glg_{data_type}_{det}_%y%m%d_v00.pha')
        gbm_files = Scraper(gbm_pattern, data_type=data_type, det=det)
        urls = gbm_files.filelist(timerange)

        return urls

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
        chkattr = ['Time', 'Instrument', 'Detector', 'Datatype']
        chklist = [x.__class__.__name__ in chkattr for x in query]
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
