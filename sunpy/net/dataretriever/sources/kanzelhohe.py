#This module was developed with funding provided by
#the Google Summer of Code 2016.
"""
This module implements Kanzelhohe Client.
"""
__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import numpy as np

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper

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
    >>> timerange = a.Time('2015/12/28 00:00:00','2015/12/28 00:03:00')
    >>> results = Fido.search(timerange, a.Instrument('kanzelhohe'), a.Wavelength(6563*u.AA))
    >>> print(results)
    [<Table length=1>
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
        wave_float = kwargs['wavelength'].min.value
        table = {6563:['halpha2k/recent', 'halph_fr'], 32768:['caiia', 'caiik_fi'], 5460:['phokada', 'bband_fi']}
        #Checking if value is close enough to a wavelength value.
        #Converting from one unit to other introduces precision errors.
        try:
            #The isclose function in numpy can only check two arrays if they
            #are close to each other within some precision. Standalone values such
            #as int, doubles etc can't be checked that way. In order to do that, we put the
            #two indiviual values in two seperate arrays da and db and then apply
            #isclose on both those arrays.
            da = list()
            da.append(wave_float)
            for wave_nums in table.keys():
                db = list()
                db.append(wave_nums)
                if np.isclose(da, db, 1e-10, 1e-10):
                    wave = wave_nums

            date_table = {6563: datetime.datetime(2000, 7, 20, 7, 45, 46), 5460: datetime.datetime(2011, 1, 7, 10, 7, 33),
                          32768: datetime.datetime(2010, 7, 31, 8, 10, 59)}
            START_DATE = date_table[wave]
            if timerange.start < START_DATE:
                raise ValueError('Earliest date for which Kanzelhohe data is available is {:%Y-%m-%d}'.format(START_DATE))
            prefix = "http://cesar.kso.ac.at/{datatype}/%Y/"
            suffix = ""
            if wave != 6563:
                suffix = "%Y%m%d/processed/"
            url_pattern = prefix + suffix + "kanz_{datatype1}_%Y%m%d_%H%M%S.fts.gz"
            crawler = Scraper(url_pattern, datatype=table[wave][0], datatype1=table[wave][1])
            if not timerange:
                return []
            result = crawler.filelist(timerange)
            return result
        except:
            raise ValueError("Enter wavelength with proper values and units")

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Global Halpha Network' #TODO: check with kanzelhohe
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
            if x.__class__.__name__ == 'Instrument' and type(x.value) is str and x.value.lower() == 'kanzelhohe':
                chk_var += 1
            if x.__class__.__name__ == 'Wavelength':
                chk_var += 1
        return chk_var == 2

