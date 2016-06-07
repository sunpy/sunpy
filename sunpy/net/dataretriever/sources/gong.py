#This module was developed with funding provided by
#the Google Summer of Code 2016.

__author__ = "Sudarshan Konge"
__email__ = "sudk1896@gmail.com"

import datetime
import re

from sunpy.net.dataretriever.client import GenericClient
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ['GONGClient', 'FARSIDEClient']

class GONGClient(GenericClient):

    def _get_url_for_timerange(self, timerange, **kwargs):
        
        
        table_physobs = {'intensity' : 'i', 'los_magnetic_field': 'b'} #legitimate physical observations.
        table_instruments = ['bb', 'ct', 'le', 'ml', 'td', 'ud', 'z'] #all possible instruments
        
        physobs_in = True if ('physobs' in kwargs.keys() and kwargs['physobs'] in table_physobs) else False #Is PhysObs entered
        instrument_in = True if 'instrument' in kwargs.keys() else False
        #Is instrument entered.
        #The or is for checking for H-alpha instruments. The observatories are single letter 
        url_pattern_1 =  'ftp://gong2.nso.edu/QR/{id}qa/%Y%m/{ObsID}bqa%y%m%d/{ObsID}bqa%y%m%dt%H%M.fits.gz'
        url_pattern_2 = 'http://gong2.nso.edu/HA/{id}haf/%Y%m/%Y%m%d/%Y%m%d%H%M%S{ObsID}h.fits.fz' #'id' here is a dummy argument
                                                                                                    #put id = '', an empty string

        generate = lambda url, id_, obs: url.format(id = id_, ObsID = obs) #Function generates urls based on instrument
        #Better than writing two regexes for each url_pattern 1 and 2.
        
        result = list() #final list of urls
        patterns = list()
        if not physobs_in:
            patterns.append(url_pattern_1), patterns.append(url_pattern_2)
        else:
            if kwargs['physobs'] == 'los_magnetic_field':
                patterns.append(url_pattern_1)
            else:
                #Differentiate on basis of wavelength
                wavelength_in = True if 'wavelength' in kwargs.keys() else False
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

        for pattern_ in patterns:
            urls = list()
            if (pattern_ == url_pattern_1):
                urls =  [generate(url_pattern_1, id, id[0].upper()) for id in instruments_to]
            elif (pattern_ == url_pattern_2):
                urls = [generate(url_pattern_2, '', id[0].upper()) for id in instruments_to]
            arr = [Scraper(pattern__).filelist(timerange) for pattern__ in urls]
            [result.extend(url) for url in arr if len(url)>0]
        
        if not timerange:
            return []
        return result

    def _makeimap(self):
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
        physobs = ['intensity', 'los_magnetic_field'] #from VSO
        instruments = ['bb', 'ct', 'le', 'ml', 'td', 'ud', 'z'] #for Magnetogram and intensity
        chk_instr, chk_physobs = 0, 0
        chk_wavelength = 0
        values = [6562, 6563, 6768]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() in instruments:
                chk_instr += 1
            if x.__class__.__name__ == 'Physobs' and x.value.lower() in physobs:
                chk_physobs += 1
            if (x.__class__.__name__ == 'Wavelength' and int(x.min.value) in values and int(x.max.value) in values and (x.unit.name).lower()=='angstrom'):
                chk_wavelength += 1
        if chk_instr == 1 or chk_physobs == 1:
            return True
        return False


class FARSIDEClient(GenericClient):

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
        #TODO: Figure out all cases where Scraper in sunpy.util would
        # and wouldn't work. Doesn't work for this website.
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
        
    
