"""
This module implements SWEPAM, MAG, SIS and EPAM Clients.
"""

import os
import datetime
from urllib.parse import urlsplit

import astropy.units as u
from astropy.time import Time, TimeDelta

from sunpy.time import TimeRange
from sunpy.util.scraper import Scraper
from sunpy.net.dataretriever.client import GenericClient

__all__ = ['SWEPAMClient', 'EPAMClient', 'MAGClient', 'SISClient']

BASEURL_PREFIX = ("ftp://ftp.swpc.noaa.gov/pub/lists/ace/%Y%m%d")


class SWEPAMClient(GenericClient):
    """
    Returns a list of URLS to ACE SWEPAM files corresponding to value of input timerange.
    URL source: `ftp://ftp.swpc.noaa.gov/pub/lists/ace/`.
    The earliest date available is from 29th July 2015.

    Parameters
    ----------
    timerange : `sunpy.time.TimeRange`
        Time range for which data is to be downloaded.

    Returns
    -------
    `list`
          List of urls corresponding to requested time range.

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
    >>> response = Fido.fetch(results)  #doctest: +SKIP
    """

    BASEURL = BASEURL_PREFIX + "_ace_swepam_1m.txt"

    def _get_time_for_url(self, urls):
        swepam = Scraper(self.BASEURL)
        times = list()
        for url in urls:
            t0 = swepam._extractDateURL(url)
            almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
            times.append(TimeRange(t0, t0 + almost_day))
        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Return URL(s) for corresponding timerange.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The time range you want the files for.

        Returns
        -------
        `list`
            The URL(s) for the corresponding timerange.
        """
        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError(
                'Earliest date for which SWEPAM data is available is '
                '{:%Y-%m-%d}'.format(START_DATE))
        swepam = Scraper(self.BASEURL)
        return sorted(swepam.filelist(timerange))

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
        query : `list`
        A list of query objects.
        Returns
        -------
        `bool` :
            Answer as to whether client can service the query

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
    The earliest date available is from 29th July 2015.

    Parameters
    ----------
    timerange : `sunpy.time.TimeRange`
        Time range for which data is to be downloaded.

    Returns
    -------
    `list`
          List of urls corresponding to requested time range.

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
    >>> response = Fido.fetch(results)  #doctest: +SKIP
    """

    BASEURL = BASEURL_PREFIX + "_ace_epam_5m.txt"

    def _get_time_for_url(self, urls):
        epam = Scraper(self.BASEURL)
        times = list()
        for url in urls:
            t0 = epam._extractDateURL(url)
            almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
            times.append(TimeRange(t0, t0 + almost_day))
        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Return URL(s) for corresponding timerange.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The time range you want the files for.

        Returns
        -------
        `list`
            The URL(s) for the corresponding timerange.
        """
        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError(
                'Earliest date for which SWEPAM data is available is '
                '{:%Y-%m-%d}'.format(START_DATE))
        epam = Scraper(self.BASEURL)
        return sorted(epam.filelist(timerange))

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
        query : `list`
        A list of query objects.
        Returns
        -------
        `bool` :
            Answer as to whether client can service the query

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
    The earliest date available is from 29th July 2015.

    Parameters
    ----------
    timerange : `sunpy.time.TimeRange`
        Time range for which data is to be downloaded.

    Returns
    -------
    `list`
          List of urls corresponding to requested time range.

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
    >>> response = Fido.fetch(results)  #doctest: +SKIP
    """

    BASEURL = BASEURL_PREFIX + "_ace_mag_1m.txt"

    def _get_time_for_url(self, urls):
        mag = Scraper(self.BASEURL)
        times = list()
        for url in urls:
            t0 = mag._extractDateURL(url)
            almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
            times.append(TimeRange(t0, t0 + almost_day))
        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Return URL(s) for corresponding timerange.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The time range you want the files for.

        Returns
        -------
        `list`
            The URL(s) for the corresponding timerange.
        """
        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError(
                'Earliest date for which SWEPAM data is available is '
                '{:%Y-%m-%d}'.format(START_DATE))
        mag = Scraper(self.BASEURL)
        return sorted(mag.filelist(timerange))

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
        query : `list`
        A list of query objects.
        Returns
        -------
        `bool` :
            Answer as to whether client can service the query

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
    The earliest date available is from 29th July 2015.

    Parameters
    ----------
    timerange : `sunpy.time.TimeRange`
        Time range for which data is to be downloaded.

    Returns
    -------
    `list`
          List of urls corresponding to requested time range.

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
    >>> response = Fido.fetch(results)  #doctest: +SKIP
    """

    BASEURL = BASEURL_PREFIX + "_ace_sis_5m.txt"

    def _get_time_for_url(self, urls):
        sis = Scraper(self.BASEURL)
        times = list()
        for url in urls:
            t0 = sis._extractDateURL(url)
            almost_day = TimeDelta(1 * u.day - 1 * u.millisecond)
            times.append(TimeRange(t0, t0 + almost_day))
        return times

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Return URL(s) for corresponding timerange.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The time range you want the files for.

        Returns
        -------
        `list`
            The URL(s) for the corresponding timerange.
        """
        START_DATE = datetime.datetime(2015, 7, 29)
        if timerange.start < START_DATE:
            raise ValueError(
                'Earliest date for which SWEPAM data is available is '
                '{:%Y-%m-%d}'.format(START_DATE))
        sis = Scraper(self.BASEURL)
        return sorted(sis.filelist(timerange))

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
        query : `list`
        A list of query objects.
        Returns
        -------
        `bool` :
            Answer as to whether client can service the query

        """
        chkattr = ['Time', 'Instrument']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'sis':
                return all(chklist)
        return False
