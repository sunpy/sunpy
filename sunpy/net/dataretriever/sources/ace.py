#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import urllib2

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ['SWEPAMClient', 'EPAMClient', 'MAGClient', 'SISClient']

class SWEPAMClient(GenericClient):
    """
    Returns a list of URLS to ACE SWEPAM files corresponding to value of input timerange.
    URL source: `ftp://ftp.swpc.noaa.gov/pub/lists/ace/`.
    The earliest date available is from 29-Jul-2015
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'swepam'
    
    Returns
    -------
    urls: list
          list of urls corresponding to requested time range.
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> results = Fido.search(a.Time('2016/5/18 00:00:00', '2016/5/20 00:03:00'), a.Instrument('swepam'))
    >>> print(results)
    [<Table length=3>
         Start Time           End Time      Source Instrument
            str19               str19         str3     str6   
    ------------------- ------------------- ------ ----------
    2016-05-18 00:00:00 2016-05-19 00:00:00    ACE     swepam
    2016-05-19 00:00:00 2016-05-20 00:00:00    ACE     swepam
    2016-05-20 00:00:00 2016-05-21 00:00:00    ACE     swepam]
    
    >>> response = Fido.fetch(results)
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """ returns list of urls corresponding
        to given TimeRange. """

        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError('Earliest date for which SWEPAM data is available is {:%Y-%m-%d}'.format(START_DATE))
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/ace/'
        total_days = (timerange.end - timerange.start).days + 1
        all_days = timerange.split(total_days)
        result = [base_url + '{date:%Y%m%d}_ace_swepam_1m.txt'.format(date = day.end) for day in all_days]
        if ((datetime.datetime.now() - timerange.end).days == 0):
            url = base_url + 'ace_swepam_1m.txt'
            result.append(url)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ACE'
        self.map_['instrument'] = 'swepam'
        self.map_['phyobs'] = 'particle-flux'

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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'swepam':
                return all(chklist)
        return False


class EPAMClient(GenericClient):
    """
    Returns a list of URLS to ACE SWEPAM files corresponding to value of input timerange.
    URL source: `ftp://ftp.swpc.noaa.gov/pub/lists/ace/`.
    The earliest date available is from 29-Jul-2015
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'epam'
    
    Returns
    -------
    urls: list
          list of urls corresponding to requested time range.
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> results = Fido.search(a.Time('2016/5/18 00:00:00', '2016/5/20 00:03:00'), a.Instrument('epam'))
    >>> print(results)
    [<Table length=3>
         Start Time           End Time      Source Instrument
           str19               str19         str3     str4   
    ------------------- ------------------- ------ ----------
    2016-05-18 00:00:00 2016-05-19 00:00:00    ACE       epam
    2016-05-19 00:00:00 2016-05-20 00:00:00    ACE       epam
    2016-05-20 00:00:00 2016-05-21 00:00:00    ACE       epam]
    
    >>> response = Fido.fetch(results)
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        START_DATE = datetime.datetime(2015, 7, 29)
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/ace/'
        if timerange.start < START_DATE:
            raise ValueError("The earliest date for which EPAM data is available is {:%Y-%m-%d}".format(START_DATE))
        total_days = (timerange.end - timerange.start).days + 1
        all_days = timerange.split(total_days)
        result = [base_url + '{date:%Y%m%d}_ace_epam_5m.txt'.format(date = day.end) for day in all_days]
        if ((datetime.datetime.now() - timerange.end).days == 0):
            url = base_url + 'ace_epam_5m.txt'
            result.append(url)
        return result 

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ACE'
        self.map_['instrument'] = 'epam'
        self.map_['phyobs'] = 'particle-flux'

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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'epam':
                return all(chklist)
        return False


class MAGClient(GenericClient):
    """
    Returns a list of URLS to ACE SWEPAM files corresponding to value of input timerange.
    URL source: `ftp://ftp.swpc.noaa.gov/pub/lists/ace/`.
    The earliest date available is from 29-Jul-2015
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'mag'
    
    Returns
    -------
    urls: list
          list of urls corresponding to requested time range.
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> results = Fido.search(a.Time('2016/5/18 00:00:00', '2016/5/20 00:03:00'), a.Instrument('mag'))
    >>> print(results)
    [<Table length=3>
         Start Time           End Time      Source Instrument
            str19               str19         str3     str6   
    ------------------- ------------------- ------ ----------
    2016-05-18 00:00:00 2016-05-19 00:00:00    ACE     mag
    2016-05-19 00:00:00 2016-05-20 00:00:00    ACE     mag
    2016-05-20 00:00:00 2016-05-21 00:00:00    ACE     mag]
    
    >>> response = Fido.fetch(results)
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """ returns list of urls corresponding
        to given TimeRange. """

        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError('Earliest date for which MAG data is available is {:%Y-%m-%d}'.format(START_DATE))
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/ace/'
        total_days = (timerange.end - timerange.start).days + 1
        all_days = timerange.split(total_days)
        result = [base_url +  '{date:%Y%m%d}_ace_mag_1m.txt'.format(date = day.end) for day in all_days]
        if ((datetime.datetime.now() - timerange.end).days == 0):
            url = base_url + 'ace_mag_1m.txt'
            result.append(url)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ACE'
        self.map_['instrument'] = 'mag'
        self.map_['phyobs'] = 'magnetic-field'

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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'mag':
                return all(chklist)
        return False

class SISClient(GenericClient):
    """
    Returns a list of URLS to ACE SWEPAM files corresponding to value of input timerange.
    URL source: `ftp://ftp.swpc.noaa.gov/pub/lists/ace/`.
    The earliest date available is from 29-Jul-2015
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'sis'
    
    Returns
    -------
    urls: list
          list of urls corresponding to requested time range.
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> results = Fido.search(a.Time('2016/5/18 00:00:00', '2016/5/20 00:03:00'), a.Instrument('sis'))
    >>> print(results)
    [<Table length=3>
         Start Time           End Time      Source Instrument
            str19               str19         str3     str6   
    ------------------- ------------------- ------ ----------
    2016-05-18 00:00:00 2016-05-19 00:00:00    ACE     sis
    2016-05-19 00:00:00 2016-05-20 00:00:00    ACE     sis
    2016-05-20 00:00:00 2016-05-21 00:00:00    ACE     sis]
    
    >>> response = Fido.fetch(results)
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """ returns list of urls corresponding
        to given TimeRange. """

        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError('Earliest date for which SIS data is available is {:%Y-%m-%d}'.format(START_DATE))
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/ace/'
        total_days = (timerange.end - timerange.start).days + 1
        all_days = timerange.split(total_days)
        result = [base_url +  '{date:%Y%m%d}_ace_sis_5m.txt'.format(date = day.end) for day in all_days]
        for day in all_days:
            url = base_url + '{date:%Y%m%d}_ace_sis_5m.txt'.format(date = day.end)
            result.append(url)
        if ((datetime.datetime.now() - timerange.end).days == 0):
            url = base_url + 'ace_sis_5m.txt'
            result.append(url)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ACE'
        self.map_['instrument'] = 'sis'
        self.map_['phyobs'] = 'particle-flux'

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
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'sis':
                return all(chklist)
        return False
