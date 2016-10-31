# This module was developed with funding provided by
# the Google Summer of Code 2016.
import re
import datetime

import numpy as np
from bs4 import BeautifulSoup

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper

from sunpy.extern.six.moves.urllib.request import urlopen

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

__all__ = ['VSMClient']


class VSMClient(GenericClient):
    """
    Returns a list of URLS to SOLIS VSM files corresponding to value of input
    timerange. URL source: `http://gong2.nso.edu/pubkeep/`.

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')

    Instrument: Fixed argument - 'vsm' ,case insensitive.

    Wavelength: Any of 6302, 8542 or 10830 Angstrom units or equivalent values in different units.
                Expects a Wavelength object.
                e.g. 6302*u.AA, 630.2*u.nm etc.

    Physobs: Optional, includes 'LOS_MAGNETIC_FIELD', "VECTOR_MAGNETIC_FIELD'
             'EQUIVALENT_WIDTH'.

    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> res = Fido.search(a.Time('2016/6/4', '2016/6/4 00:10:00'), a.Instrument('vsm'),
                          a.Wavelength(6302*u.AA), a.Physobs('VECTOR_MAGNETIC_FIELD'))
    >>> print (res)
    [<Table length=6>
         Start Time           End Time      Source Instrument
           str19               str19         str5     str3
    ------------------- ------------------- ------ ----------
    2015-01-03 00:00:00 2015-01-04 00:00:00  SOLIS        vsm
    2015-01-04 00:00:00 2015-01-05 00:00:00  SOLIS        vsm
    2015-01-05 00:00:00 2015-01-06 00:00:00  SOLIS        vsm
    2015-01-06 00:00:00 2015-01-07 00:00:00  SOLIS        vsm
    2015-01-07 00:00:00 2015-01-08 00:00:00  SOLIS        vsm
    2015-01-08 00:00:00 2015-01-09 00:00:00  SOLIS        vsm]
    """

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        table_physobs = [
            'EQUIVALENT_WIDTH', 'LOS_MAGNETIC_FIELD', 'VECTOR_MAGNETIC_FIELD'
        ]
        table_wave = {6302: 'v72', 8542: 'v82', 10830: 'v22'}

        physobs_in = ('physobs' in kwargs.keys() and
                      kwargs['physobs'] in table_physobs)
        url_pattern = 'http://gong2.nso.edu/pubkeep/{dtype}/%Y%m/k4{dtype}%y%m%d/k4{dtype}%y%m%dt%H%M%S{suffix}'

        wave_float = kwargs['wavelength'].min.value

        # The isclose function in numpy can only check two arrays if they are
        # close to each other within some precision. Standalone values such as
        # int, doubles etc can't be checked that way. In order to do that, we
        # put the two indiviual values in two seperate arrays da and db and
        # then apply isclose on both those arrays.
        da = list()
        da.append(wave_float)
        for wave_nums in table_wave.keys():
            db = list()
            db.append(wave_nums)
            if np.isclose(da, db, 1e-10, 1e-10):
                wave = wave_nums
                break

        if wave is None:
            raise ValueError("Enter correct wavelength values and units")

        def download_from_nso(dtype, suffix, timerange):
            tmp_pattern = 'k4{dtype}%y%m%dt%H%M%S{suffix}'
            result = list()
            base_url = 'http://gong2.nso.edu/pubkeep/{dtype}/{date:%Y%m}/k4{dtype}{date:%y%m%d}/'
            total_days = (timerange.end - timerange.start).days + 1
            all_dates = timerange.split(total_days)
            for day in all_dates:
                try:
                    html = urlopen(base_url.format(dtype=dtype, date=day.end))
                    soup = BeautifulSoup(html)
                    for link in soup.findAll("a"):
                        crawler = Scraper(
                            tmp_pattern, dtype=dtype, suffix=suffix)
                        simple_pattern = '%y%m%d_%H%M%S'
                        url = str(link.get('href'))
                        if (crawler._URL_followsPattern(url)):
                            tmp = re.sub('k4' + dtype, '', url)
                            tmp = re.sub('[a-zA-Z]', '', tmp)
                            crawler_ = Scraper(simple_pattern)
                            if (timerange.start <=
                                    crawler_._extractDateURL(tmp) <=
                                    timerange.end):
                                result.append(
                                    base_url.format(
                                        dtype=dtype, date=day.end) + url)
                except:
                    pass
            return result

        start_date = {
            6302: datetime.datetime(2003, 8, 21),
            8542: datetime.datetime(2003, 8, 26),
            10830: datetime.datetime(2004, 11, 4)
        }
        START_DATE = start_date[wave]

        if timerange.start < START_DATE:
            raise ValueError(
                'Earliest date for which SOLIS VSM data is available is {:%Y-%m-%d}'.
                format(START_DATE))

        result = list()
        if wave == 6302:
            if not physobs_in:
                result.extend(download_from_nso('v72', '.fts.gz', timerange))
                result.extend(
                    download_from_nso('v93', '_FDISK.fts.gz', timerange))
            else:
                if kwargs['physobs'] == 'VECTOR_MAGNETIC_FIELD':
                    result.extend(
                        download_from_nso('v93', '_FDISK.fts.gz', timerange))
                else:
                    result.extend(
                        download_from_nso('v72', '.fts.gz', timerange))
        else:
            result.extend(
                download_from_nso(table_wave[wave]), '.fts.gz', timerange)
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'SOLIS'
        self.map_['instrument'] = 'vsm'

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
        chk_var = 0
        for x in query:
            if x.__class__.__name__ == 'Instrument' and type(
                    x.value) is str and x.value.lower() == 'vsm':
                chk_var += 1
            if x.__class__.__name__ == 'Wavelength':
                chk_var += 1
        return chk_var == 2
