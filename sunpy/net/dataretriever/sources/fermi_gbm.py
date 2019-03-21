from sunpy.util.scraper import Scraper
from ..client import GenericClient

__all__ = ['GBMClient']


class GBMClient(GenericClient):
    """
    This GBMClient is for the Gamma-Ray Burst Monitor (GBM) instrument
    aboard the Fermi satellite. Although GBMs primary objective is to
    detect gamma-ray bursts, it provides high quailty high energy solar
    flare observations.

    The instrument consists of 12 Sodium Iodide (NaI) scintillation
    detectors, which are sensitive to an energy range of 4keV to 1MeV.
    At any one time, 6 of the NaI detectors are Sunward facing.
    The detectors are numbered 'n1' to 'n11'. This client supports the user
    to choose which detector to use through the `a.Detector <sunpy.net.attrs.Detector>` attribute.
    The default detector is 'n5'.

    The GBM data comes in daily version files in two formats:

        * CSPEC - counts accumulated every  4.096 seconds in 128 energy channels for each detector.
        * CTIME - counts accumulated every 0.256 seconds in 8 energy channels

    Both of which can be accessed through the attrs `a.Resolution <sunpy.net.attrs.Resolution>`.
    The default data type is CSPEC unless the user defines.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> res = Fido.search(a.Time('2015-06-21 00:00', '2015-06-23 23:59'),
    ...                   a.Instrument('gbm'), a.Detector('n3'),
    ...                   a.Resolution('ctime')) #doctest: +REMOTE_DATA
    >>> print(res) #doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the GBMClient:
         Start Time           End Time      Source Instrument Wavelength
           str19               str19         str5     str3       str3
    ------------------- ------------------- ------ ---------- ----------
    2015-06-21 00:00:00 2015-06-23 23:59:00  FERMI        GBM        nan
    2015-06-21 00:00:00 2015-06-23 23:59:00  FERMI        GBM        nan
    2015-06-21 00:00:00 2015-06-23 23:59:00  FERMI        GBM        nan
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Returns the url for Fermi/GBM data for the given date.

        Parameters
        ----------
        timerange : `sunpy.time.TimeRange`
            The time range for which to download the data.

        Returns
        -------
        `str`:
            The url(s) for time of interest.
        """
        # Checks if detector keyword
        # If not defaults to detector 5
        if 'detector' in kwargs:
            det = _check_detector(kwargs['detector'])
        else:
            det = 'n5'

        # Check for resolution keyword - either CSPEC or CTIME
        # Default type is CSPEC
        if 'resolution' in kwargs:
            data_type = _check_type(kwargs['resolution'])
        else:
            data_type = 'cspec'

        gbm_pattern = ('https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'
                       '%Y/%m/%d/current/glg_{data_type}_{det}_%y%m%d_v00.pha')
        gbm_files = Scraper(gbm_pattern, data_type=data_type, det=det)
        urls = gbm_files.filelist(timerange)

        return urls

    def _makeimap(self):
        """
        Helper function used to hold information about source.
        """
        self.map_['source'] = 'FERMI'
        self.map_['instrument'] = 'GBM'
        self.map_['physobs'] = 'flux'
        self.map_['provider'] = 'NASA'

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Answers whether a client can service the query.

        Parameters
        ----------
        query : `list`
            A list of of query objects.

        Returns
        -------
        `bool`
            `True` if this client can service the query, otherwise `False`.
        """
        chkattr = ['Time', 'Instrument', 'Detector', 'Resolution']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'gbm':
                return all(chklist)
        return False


def _check_detector(detector, **kwargs):
    """
    checks to see if detector is in right format.
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
    checks is datatype is either "CSPEC" or "CTIME".
    """
    if not isinstance(datatype, str):
        raise ValueError('{} is not str - either cspec or ctime'.format(datatype))

    if datatype.lower() != 'cspec' and datatype.lower() != 'ctime':
        raise ValueError('{} not value datatype - either cspec or ctime'.format(datatype))
    else:
        return datatype.lower()
