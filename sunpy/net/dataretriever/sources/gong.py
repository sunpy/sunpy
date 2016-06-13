#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import re
import urllib2

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ['GONGClient', 'FARSIDEClient']

class GONGClient(GenericClient):
    """
    Returns a list of URLS to GONG files corresponding to value of input timerange.
    URL source: `ftp://gong2.nso.edu/QR/` and `http://gong2.nso.edu/HA/haf/`.
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Arguments are any one from ['bb', 'ct', 'le', 'ml', 'td', 'ud','z']

    Physobs: Two arguments or none 'intensity' and 'los_magnetic_field'

    GONGClient expects at least one argument from Physobs or Instrument

    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> res = Fido.search(a.Time('2016/6/4', '2016/6/4 00:10:00'), a.Physobs('intensity'))
    >>> print (res)
    [<Table length=32>
     Start Time           End Time      Source     Instrument    
       str19               str19         str4        str18       
    ------------------- ------------------- ------ ------------------
    2016-06-04 00:00:00 2016-06-05 00:00:00   GONG Data not Available
    2016-06-05 00:00:00 2016-06-06 00:00:00   GONG Data not Available
    2016-06-06 00:00:00 2016-06-07 00:00:00   GONG Data not Available
    2016-06-07 00:00:00 2016-06-08 00:00:00   GONG Data not Available
    2016-06-08 00:00:00 2016-06-09 00:00:00   GONG Data not Available
    2016-06-09 00:00:00 2016-06-10 00:00:00   GONG Data not Available
    2016-06-10 00:00:00 2016-06-11 00:00:00   GONG Data not Available
    2016-06-11 00:00:00 2016-06-12 00:00:00   GONG Data not Available
    2016-06-12 00:00:00 2016-06-13 00:00:00   GONG Data not Available
    2016-06-13 00:00:00 2016-06-14 00:00:00   GONG Data not Available
                    ...                 ...    ...                ...
    2016-06-26 00:00:00 2016-06-27 00:00:00   GONG Data not Available
    2016-06-27 00:00:00 2016-06-28 00:00:00   GONG Data not Available
    2016-06-28 00:00:00 2016-06-29 00:00:00   GONG Data not Available
    2016-06-29 00:00:00 2016-06-30 00:00:00   GONG Data not Available
    2016-06-30 00:00:00 2016-07-01 00:00:00   GONG Data not Available
    2016-07-01 00:00:00 2016-07-02 00:00:00   GONG Data not Available
    2016-07-02 00:00:00 2016-07-03 00:00:00   GONG Data not Available
    2016-07-03 00:00:00 2016-07-04 00:00:00   GONG Data not Available
    2016-07-04 00:00:00 2016-07-05 00:00:00   GONG Data not Available
    2016-07-05 00:00:00 2016-07-06 00:00:00   GONG Data not Available]
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        #TO-DO: Figure out where scraper would and wouldn't work
        table_physobs = {'INTENSITY' : 'i', 'LOS_MAGNETIC_FIELD': 'b'} #legitimate physical observations.
        table_instruments = ['bb', 'ct', 'le', 'ml', 'td', 'ud','z'] #For magnetogram and intensity
        
        physobs_in = True if ('physobs' in kwargs.keys() and kwargs['physobs'] in table_physobs.keys()) else False #Is PhysObs entered
        instrument_in = True if ('instrument' in kwargs.keys() and kwargs['instrument'] in table_instruments) else False
        wavelength_in = True if 'wavelength' in kwargs.keys() else False
        #Is instrument entered.
        url_pattern_1 =  'ftp://gong2.nso.edu/QR/{id}qa/%Y%m/{ObsID}bqa%y%m%d/{ObsID}bqa%y%m%dt%H%M.fits.gz'
        url_pattern_2 = 'http://gong2.nso.edu/HA/haf/%Y%m/%Y%m%d/%Y%m%d%H%M%S{ObsID}h.fits.fz'
        
        result = list() #final list of urls
        patterns = list()
        if not physobs_in:
            patterns.append(url_pattern_1), patterns.append(url_pattern_2)
        else:
            if kwargs['physobs'] == 'LOS_MAGNETIC_FIELD':
                patterns.append(url_pattern_1)
            else:
                #Differentiate on basis of wavelength
                if not wavelength_in:
                    patterns.append(url_pattern_1), patterns.append(url_pattern_2)
                else:
                    wave = int(kwargs['wavelength'].min.value)
                    if (6562 <= wave <= 6563):
                        patterns.append(url_pattern_2)
                    elif wave == 6768:
                        patterns.append(url_pattern_1)


        #All valid patterns to be downloaded are in the patterns list.
        instruments_to = list() #The instruments from which user wants to download
        if not instrument_in:
            instruments_to.extend(table_instruments)
        else:
            instruments_to.append(kwargs['instrument'])

        
        def download_from_nso(id, ObsID, time_range):
            #Cannot scrape urls from NSO. The logic for downloading is
            #dependent on NSO's policy of uploading data. Every day staring from 00:04
            #fits files are uploaded regularly at an interval of 10 minutes.
            result = list()
            base_url = 'ftp://gong2.nso.edu/QR/{id}qa/{date:%Y%m}/{ObsID}{id}qa{date:%y%m%d}/{ObsID}{id}qa{date:%y%m%d}t{date:%H%M}.fits.gz'
            total_days = (time_range.end - time_range.start).days + 1
            all_dates = time_range.split(total_days)
            for day in all_dates:
                today = datetime.datetime(day.end.year, day.end.month, day.end.day,
                                          0, 4, 0)
                tomorrow = today
                tomorrow += datetime.timedelta(days=1)
                while (today < tomorrow):
                    if (time_range.start <= today <= time_range.end):
                        result.append(base_url.format(id=id, ObsID=ObsID, date=today))
                    today += datetime.timedelta(seconds = 600)#10 minutes, 600 seconds.
            result = list(set(result)) #Remove duplicates, for safety.
            return result
                
            
        
        for pattern_ in patterns:
            urls = list()
            if (pattern_ == url_pattern_1):
                if not physobs_in:
                    for instr in instruments_to:
                        arr = download_from_nso('i',instr,timerange)
                        urls.extend(arr)
                        arr = download_from_nso('b',instr,timerange)
                        urls.extend(arr)
                else:
                    for instr in instruments_to:
                        urls.extend(download_from_nso(table_physobs[kwargs['physobs']],instr,timerange))
                result.extend(urls)
            elif (pattern_ == url_pattern_2):
                urls = [url_pattern_2.format(ObsID=id[0].upper()) for id in instruments_to]
                arr = [Scraper(pattern__).filelist(timerange) for pattern__ in urls]
                [result.extend(url) for url in arr if len(url)>0]
        
        if not timerange:
            return []
        final_result = list()
        #Check for illegitimate urls.
        for urls in result:
            try:
                check = urllib2.urlopen(urls)
                final_result.append(urls)
            except urllib2.HTTPError,e:
                pass
            except urllib2.URLError,e:
                pass
        return final_result

    def _makeimap(self):
        """
        Helper Function: used to hold information about source.
        """   
        self.map_['source'] = 'GONG'
        
        
    
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
        chkattr = ['Time', 'Instrument', 'Physobs', 'Wavelength']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        physobs = ['INTENSITY', 'LOS_MAGNETIC_FIELD'] #from VSO
        instruments = ['bb', 'ct', 'le', 'ml', 'td', 'ud', 'z'] #for Magnetogram and intensity
        chk_instr, chk_physobs = 0, 0
        chk_wavelength = 0
        values = [6562, 6563, 6768]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() in instruments:
                chk_instr += 1
            if x.__class__.__name__ == 'Physobs' and x.value in physobs:
                chk_physobs += 1
            if (x.__class__.__name__ == 'Wavelength' and int(x.min.value) in values and int(x.max.value) in values and (x.unit.name).lower()=='angstrom'):
                chk_wavelength += 1
        if chk_instr == 1 or chk_physobs == 1:
            return True
        return False


class FARSIDEClient(GenericClient):
    """
    Returns a list of URLS to FARSIDE files corresponding to value of input timerange.
    URL source: `http://farside.nso.edu/oQR/fqo/`.
    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value - TimeRange('2015-12-30 00:00:00','2015-12-31 00:01:00')
    
    Instrument: Fixed argument = 'farside'
            
    Returns
    -------
    urls: list
    list of urls corresponding to requested time range.
    
    Examples
    --------
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a
    >>> results = Fido.search(Time('2015/4/2','2015/4/4'), a.Instrument('farside'))
    >>> print(results)
    [<Table length=4>
         Start Time           End Time      Source Instrument
           str19               str19         str4     str7   
    ------------------- ------------------- ------ ----------
    2015-04-02 00:00:00 2015-04-03 00:00:00   GONG    farside
    2015-04-03 00:00:00 2015-04-04 00:00:00   GONG    farside
    2015-04-04 00:00:00 2015-04-05 00:00:00   GONG    farside
    2015-04-05 00:00:00 2015-04-06 00:00:00   GONG    farside]
    
    >>> response = Fido.fetch(results)
    """
    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        START_DATE = datetime.datetime(2006, 5, 13, 12, 0, 0)
        if timerange.start < START_DATE:
            raise ValueError('Earliest date for which FARSIDE data is available is {:%Y-%m-%d}'.format(START_DATE))
        url_pattern = 'http://farside.nso.edu/oQR/fqo/{:%Y%m}/mrfqo{:%y%m%d}/mrfqo{:%y%m%d}t{:%H%M}.fits'
        tot_days = (timerange.end - timerange.start).days
        all_days = timerange.split(tot_days)
        result = list()
        for dt in all_days:
            times = [datetime.datetime(dt.start.year, dt.start.month, dt.start.day, 0, 0),
                     datetime.datetime(dt.start.year, dt.start.month, dt.start.day, 12, 0)]
            [result.append(url_pattern.format(dt_, dt_, dt_, dt_)) for dt_ in times]
        if not timerange:
            return []
        return result

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'GONG'
        self.map_['instrument'] = 'farside'

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
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'farside':
                return all(chklist)
        return False
        
