"""
This module implements SWEPAM, MAG, SIS and EPAM Clients.
"""
# This module was developed with funding provided by
# the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import os
import datetime
from urllib.parse import urlsplit

import astropy.units as u
from astropy.time import Time, TimeDelta
from sunpy.time import TimeRange
from sunpy.net.dataretriever.client import GenericClient

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
    >>> results = Fido.search(a.Time('2016/5/18 00:00:00', '2016/5/20 00:03:00'),
    ...                       a.Instrument('swepam'))   #doctest: +REMOTE_DATA
    >>> print(results)  #doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the SWEPAMClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str3     str6       str3
    ------------------- ------------------- ------ ---------- ----------
    2016-05-18 00:00:00 2016-05-18 23:59:59    ACE     swepam        nan
    2016-05-19 00:00:00 2016-05-19 23:59:59    ACE     swepam        nan
    2016-05-20 00:00:00 2016-05-20 23:59:59    ACE     swepam        nan
    <BLANKLINE>
    <BLANKLINE>

    >>> response = Fido.fetch(results)  #doctest: +REMOTE_DATA
    """

    def _get_time_for_url(self, urls):
        times = []
        for uri in urls:
            uripath = urlsplit(uri).path

            # Extract the yymmdd or yyyymmdd timestamp
            datestamp = os.path.splitext(os.path.split(uripath)[1])[0][:8]

            start = Time.strptime(datestamp, "%Y%m%d")
            almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
            times.append(TimeRange(start, start + almost_day))

        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        """ returns list of urls corresponding
        to given TimeRange. """

        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError(
                'Earliest date for which SWEPAM data is available is '
                '{:%Y-%m-%d}'.format(START_DATE))
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/ace/'
        total_days = int(timerange.days / u.d + 1)
        all_days = timerange.split(total_days)
        result = [
            base_url +
            '{date}_ace_swepam_1m.txt'.format(
                date=str(day.end).split('T')[0].replace('-', '')) for day in all_days]
        time_dif = Time(datetime.datetime.now()) - timerange.end
        time_dif.format = 'datetime'
        if time_dif.value.days == 0:
            url = base_url + '_ace_swepam_1m.txt'
            result.append(url)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ACE'
        self.map_['instrument'] = 'swepam'
        self.map_['phyobs'] = 'PARTICLE_FLUX'

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
    Returns a list of URLS to ACE EPAM files corresponding to value of input timerange.
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
    >>> results = Fido.search(a.Time('2016/5/18 00:00:00', '2016/5/20 00:03:00'),
    ...                       a.Instrument('epam')) #doctest: +REMOTE_DATA
    >>> print(results)  #doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the EPAMClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str3     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2016-05-18 00:00:00 2016-05-18 23:59:59    ACE       epam        nan
    2016-05-19 00:00:00 2016-05-19 23:59:59    ACE       epam        nan
    2016-05-20 00:00:00 2016-05-20 23:59:59    ACE       epam        nan
    <BLANKLINE>
    <BLANKLINE>

    >>> response = Fido.fetch(results)  #doctest: +REMOTE_DATA
    """

    def _get_time_for_url(self, urls):
        times = []
        for uri in urls:
            uripath = urlsplit(uri).path

            # Extract the yymmdd or yyyymmdd timestamp
            datestamp = os.path.splitext(os.path.split(uripath)[1])[0][:8]

            start = Time.strptime(datestamp, "%Y%m%d")
            almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
            times.append(TimeRange(start, start + almost_day))

        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        START_DATE = datetime.datetime(2015, 7, 29)
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/ace/'
        if timerange.start < START_DATE:
            raise ValueError(
                "The earliest date for which EPAM data is available is "
                "{:%Y-%m-%d}".format(START_DATE))
        total_days = int(timerange.days / u.d + 1)
        all_days = timerange.split(total_days)
        result = [
            base_url +
            '{date}_ace_epam_5m.txt'.format(
                date=str(day.end).split('T')[0].replace('-', '')) for day in all_days]
        time_dif = Time(datetime.datetime.now()) - timerange.end
        time_dif.format = 'datetime'
        if time_dif.value.days == 0:
            url = base_url + '_ace_epam_5m.txt'
            result.append(url)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ACE'
        self.map_['instrument'] = 'epam'
        self.map_['phyobs'] = 'PARTICLE_FLUX'

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
    Returns a list of URLS to ACE MAG files corresponding to value of input timerange.
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
    >>> results = Fido.search(a.Time('2016/5/18 00:00:00', '2016/5/20 00:03:00'),
    ...                       a.Instrument('mag'))  #doctest: +REMOTE_DATA
    >>> print(results)  #doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the MAGClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str3     str3       str3
    ------------------- ------------------- ------ ---------- ----------
    2016-05-18 00:00:00 2016-05-18 23:59:59    ACE        mag        nan
    2016-05-19 00:00:00 2016-05-19 23:59:59    ACE        mag        nan
    2016-05-20 00:00:00 2016-05-20 23:59:59    ACE        mag        nan
    <BLANKLINE>
    <BLANKLINE>

    >>> response = Fido.fetch(results)  #doctest: +REMOTE_DATA
    """

    def _get_time_for_url(self, urls):
        times = []
        for uri in urls:
            uripath = urlsplit(uri).path

            # Extract the yymmdd or yyyymmdd timestamp
            datestamp = os.path.splitext(os.path.split(uripath)[1])[0][:8]

            start = Time.strptime(datestamp, "%Y%m%d")
            almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
            times.append(TimeRange(start, start + almost_day))

        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        """ returns list of urls corresponding
        to given TimeRange. """

        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError(
                'Earliest date for which MAG data is available is {:%Y-%m-%d}'.format(START_DATE))
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/ace/'
        total_days = int(timerange.days / u.d + 1)
        all_days = timerange.split(total_days)
        result = [
            base_url +
            '{date}_ace_mag_1m.txt'.format(
                date=str(day.end).split('T')[0].replace('-', '')) for day in all_days]
        time_dif = Time(datetime.datetime.now()) - timerange.end
        time_dif.format = 'datetime'
        if time_dif.value.days == 0:
            url = base_url + '_ace_mag_1m.txt'
            result.append(url)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ACE'
        self.map_['instrument'] = 'mag'
        self.map_['phyobs'] = 'MAGNETIC_FIELD'

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
    Returns a list of URLS to ACE SIS files corresponding to value of input timerange.
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
    >>> results = Fido.search(a.Time('2016/5/18 00:00:00', '2016/5/20 00:03:00'),
    ...                       a.Instrument('sis'))  #doctest: +REMOTE_DATA
    >>> print(results)  #doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the SISClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str3     str3       str3
    ------------------- ------------------- ------ ---------- ----------
    2016-05-18 00:00:00 2016-05-18 23:59:59    ACE        sis        nan
    2016-05-19 00:00:00 2016-05-19 23:59:59    ACE        sis        nan
    2016-05-20 00:00:00 2016-05-20 23:59:59    ACE        sis        nan
    <BLANKLINE>
    <BLANKLINE>

    >>> response = Fido.fetch(results)  #doctest: +REMOTE_DATA
    """

    def _get_time_for_url(self, urls):
        times = []
        for uri in urls:
            uripath = urlsplit(uri).path

            # Extract the yymmdd or yyyymmdd timestamp
            datestamp = os.path.splitext(os.path.split(uripath)[1])[0][:8]

            start = Time.strptime(datestamp, "%Y%m%d")
            almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
            times.append(TimeRange(start, start + almost_day))

        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        """ returns list of urls corresponding
        to given TimeRange. """

        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError(
                'Earliest date for which SIS data is available is {:%Y-%m-%d}'.format(START_DATE))
        base_url = 'ftp://ftp.swpc.noaa.gov/pub/lists/ace/'
        total_days = int(timerange.days / u.d + 1)
        all_days = timerange.split(total_days)
        result = []
        for day in all_days:
            url = base_url + '{date}_ace_sis_5m.txt'.format(
                date=str(day.end).split('T')[0].replace('-', ''))
            result.append(url)
        time_dif = Time(datetime.datetime.now()) - timerange.end
        time_dif.format = 'datetime'
        if time_dif.value.days == 0:
            url = base_url + 'ace_sis_5m.txt'
            result.append(url)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'ACE'
        self.map_['instrument'] = 'sis'
        self.map_['phyobs'] = 'PARTICLE_FLUX'

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
