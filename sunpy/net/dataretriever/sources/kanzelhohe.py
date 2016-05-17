#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import urllib2


from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper

from sunpy.time import TimeRange

__all__ = ['KanzelhoheClient']

class KanzelhoheClient(GenericClient):
    """
    Returns a list of URLS to Kanzelhohe H-alpha files corresponding to value of input timerange.
    URL source: `http://cesar.kso.ac.at/`.

    The earliest data for H-alpha 2k - 20-Jul-2000
                          Ca-II k - 31-Jul-2010
                          Continuum - 7-Jan-2011

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'kanzelhohe'

    Wavelength: Fixed argument = astropy.units.quantity.Quantity
                The physical value of wavelength will belong to any of [5460, 6563, 32768]
                and units will be Angstroms.

    Returns
    -------
    urls: list
    list of urls corresponding to requested time range.

    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a

    >>> results = Fido.search(a.Time('2015/12/28 00:00:00', '2015/12/28 00:03:00'), a.Instrument('kanzelhohe'), a.Wavelength(6563*u.AA))
    >>> print(results)
    >>> [<Table length=1>
        Start Time           End Time              Source        Instrument
        str19               str19                str21           str10   
    ------------------- ------------------- --------------------- ----------
    2015-12-28 00:00:00 2015-12-29 00:00:00 Global Halpha Network Kanzelhohe]
    
    >>> response = Fido.fetch(results)
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        wave = int(kwargs['Wavelength'].min.value)
        table = {6563:'halpha2k/recent', 32768:'caiia', 5460:'phokada'}
        table1 = {6563:'halph_fr', 32768:'caiik_fi', 5460:'bband_fi'}
        datatype = table[wave]
        datatype1 = table1[wave]
        if (wave==6563):
            START_DATE = datetime.datetime(2000, 7, 20, 7, 45, 46)
        elif (wave==5460):
            START_DATE = datetime.datetime(2011, 1, 7, 10, 7, 33)
        elif (wave==32768):
            START_DATE = datetime.datetime(2010, 7, 31, 8, 10, 59)
        if timerange.start < START_DATE:
            raise ValueError('Earliest date for which SWAP data is available is {:%Y-%m-%d}'.format(START_DATE))
        prefix = "http://cesar.kso.ac.at/{datatype}/"
        if (wave==6563):
            suffix = "%Y/kanz_{datatype1}_%Y%m%d_%H%M%S.fts.gz"
        else:
            suffix = "%Y/%Y%m%d/processed/kanz_{datatype1}_%Y%m%d_%H%M%S.fts.gz"
        url_pattern = prefix + suffix
        crawler = Scraper(url_pattern, datatype = datatype, datatype1 = datatype1)
        if not timerange:
            return []
        result = crawler.filelist(timerange)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Global Halpha Network'
        self.map_['instrument'] = 'Kanzelhohe HA2'
        self.map_['phyobs'] = 'irradiance'
        self.map_['provider'] = 'Kanzelhohe'

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
        chkattr = ['Time', 'Instrument', 'Wavelength']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        chk_var = 0
        values = [6563, 5460, 32768]
        for x in query:
            if (x.__class__.__name__ == 'Instrument' and x.value.lower() == 'kanzelhohe'):
                chk_var += 1
            if (x.__class__.__name__ == 'Wavelength' and int(x.min.value) in values and (x.unit.name).lower()=='angstrom') and int(x.max.value) in values:
                chk_var += 1
        if (chk_var==2):
            return True
        return False
            
        
    
