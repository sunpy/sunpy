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
        instrument_in = True if 'instrument' in kwargs.keys() else False #Is instrument entered.
        
        url_pattern_1 =  'ftp://gong2.nso.edu/QR/{id}qa/%Y%m/{ObsID}bqa%y%m%d/{ObsID}bqa%y%m%dt%H%M.fits.gz'
        url_pattern_2 = 'http://gong2.nso.edu/HA/{id}haf/%Y%m/%Y%m%d/%Y%m%d%H%M%S{ObsID}h.fits.fz' #'id' here is a dummy argument

        generate = lambda url, id_, obs: url.format(id = id_, ObsID = obs) #Function generates urls based on instrument
        
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
                
                
                
        # Create a lambda for generating urls for 1. and for 2.
        
        #table = {'bb':'B', 'ct':'ct', 'le', 'ml', 'td', 'ud', 'z'} #for H-alpha


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
        chkattr = ['Time', 'Instrument', 'Physobs']
        chklist = [x.__class__.__name__ in chkattr for x in query]
        physobs = ['intensity', 'los_magnetic_field'] #from VSO
        instruments = ['bb', 'ct', 'le', 'ml', 'td', 'ud', 'z'] #for Magnetogram and intensity
        chk_instr, chk_physobs = 0, 0
        chk_source, chk_wavelength = 0, 0
        values = [6562, 6563, 6768]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() in instruments:
                chk_instr += 1
            if x.__class__.__name__ == 'Physobs' and x.value.lower() in physobs:
                chk_physobs += 1
            if x.__class__.__name__ == 'Source' and x.value.lower() == 'gong':
                chk_source += 1
            if (x.__class__.__name__ == 'Wavelength' and int(x.min.value) in values and int(x.max.value) in values and (x.unit.name).lower()=='angstrom'):
                chk_wavelength += 1
        if chk_instr == 1 or chk_physobs == 1:
            return True
        else:
            return chk_source == 1


class FARSIDEClient(GenericClient):

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        returns list of urls corresponding to given TimeRange.
        """
        START_DATE = datetime.datetime(2006, 5, 13, 12, 0, 0)
        if timerange.start < START:
            raise ValueError('Earliest date for which XRT data is available is {:%Y-%m-%d}'.format(START_DATE))
        url_pattern = 'http://farside.nso.edu/oQR/fqo/%Y%m/mrfqo%y%m%dt%H%M.fits'
        crawler = Scraper(url_pattern)
        if not timerange:
            return []
        result = Scraper.filelist(timerange)
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
        
    
