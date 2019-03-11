from urllib.parse import urljoin

from ..client import GenericClient

__all__ = ['GBMClient']

class GBMCLient(GenericClient):
	def _get_url_for_timerange(self, timerange, **kwargs):
		"""

		Returns the url for Fermi/GBM data for the given date.
       
        baseurl = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'


		Parameters
		----------
		timerange : sunpy.time.TimeRange
			time range for which to download the data

		Returns
		-------

		url for time of interest	
        
        """

        
        if 'detector' in kwargs:
        	det = _check_det(kwargs['detector'])


        	days = timerange.get_dates()
        	urls = []
        	for day in days:
        		urls.append(self._get_url_for_date(day, **kwargs))

        	
        	return urls

    def _get_url_for_date(self, date, **kwargs):
    	"""
    	Returns URL dats of interest

    	Parameters
    	----------
    	date : 'datetime.date'

    	Returns
    	-------
		url : str
			the url at the date of interest

		"""
		det = 'n5'
		filename = 'glg_cspec_' + det + date.strftime('_%y%m%d_') + 'v00.pha'
		url_path = urljoin(date.strftime('%Y/%m/%d/') + 'current/', filename)
		base_url = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'

		return urljoin(base_url, url_path)




    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['source'] = 'Fermi'
        self.map_['instrument'] = 'GBM'
        self.map_['physobs'] = 'flux'
        self.map_['provider'] = 'nasa'



    @classmethod
    def _can_handle_query(cls, *query):
        """
        Answers whether client can service the query.
        Parameters
        ----------
        query : list of query objects
        Returns
        -------
        boolean
            answer as to whether client can service the query
        """
        chkattr =  ['Time', 'Instrument']
        chklist =  [x.__class__.__name__ in chkattr for x in query]
        for x in query:
            if x.__class__.__name__ == 'Instrument' and x.value.lower() == 'gbm':
                return all(chklist)
        return False





