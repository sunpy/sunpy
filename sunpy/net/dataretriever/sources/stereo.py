"""STEREO Dataretriever sources subclass definitions"""

__author__ = "Ankit Kumar"
__email__  = "ankitkmr.iitk@gmail.com"

#This module was developed under funding provided by
#Google Summer of Code 2015

import datetime

from astropy import units as u
from ..client import GenericClient

from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange



__all__ = ['SEPTClient', 'HETClient', 'SITClient', 'PLASTICClient', 'MAGClient', 'LETClient']


class SEPTClient(GenericClient):
	def _get_url_for_timerange(cls,timerange, stereo_spacecraft = 'ahead', duration_of_average = 10*u.min,
							   specie = 'element', sensor_pointing = 'asun'):
		"""
		Returns list of URLS to STEREO SEPT data files corresponding to value of input timerange.
		URL Source : http://www2.physik.uni-kiel.de/stereo/data/sept/level2/

		The earliest data available is from 20-Jan-2007.

		Parameters
		----------
		timerange: sunpy.time.TimeRange
			time range for which data is to be downloaded.
			Example value -  TimeRange('2007-01-20','2015-01-01')	

		stereo_spacecraft: string	
			Default value - ahead
			Possible values - ahead, behind    # corresponding to spacecraft location

		duration_of_average: astropy units quantity
		Default value - 10*u.min
			Possible values - 1*u.min, 10*u.min, 1*u.h, 1*u.d		
			#corresponding to duration over which data is averaged

		specie:  string
			Default value - element
			Possible values - element, ion

		sensor_pointing: string
			Default value - asun
			Possible values - asun, sun, north, south, omni

		Returns
		-------
		urls : list
			list of URLs corresponding to the requested time range

		"""

		base_url = 'http://www2.physik.uni-kiel.de/stereo/data/sept/level2/'

		possible_spacecraft = ['ahead','behind']
		possible_specie = ['element', 'ion']
		possible_duration = [1*u.min, 10*u.min, 1*u.h, 1*u.d]
		possible_sensor = ['asun', 'sun', 'north','south','omni']

		dict_specie = {'element': 'ele', 'ion':'ion'}
		dict_duration = {1*u.min: '1min', 10*u.min:'10min', 1*u.h: '1h', 1*u.d: '1d'}	


		#Parameter Validation
		if timerange.start < datetime.datetime(2007,01,20):
			raise ValueError('Earliest date for which SEPT data is available is 2007-01-20')

		if stereo_spacecraft not in possible_spacecraft:
			raise ValueError('Possible stereo_spacecraft values: ' + ','.join(possible_spacecraft))

		if specie not in possible_specie:
			raise ValueError('Possible specie values: ' + ','.join(possible_specie))

		if duration_of_average not in possible_duration:
			raise ValueError('Possible duration_of_average values as astropy unit quantities: ' + ','.join([str(i) for i in possible_duration]))

		if sensor_pointing not in possible_sensor:
			raise ValueError('Possible sensor_pointing values: ' + ','.join(possible_sensor))



		base_pattern = base_url + '{stereo_spacecraft}/{duration_of_average}/%Y/'
		filename_pattern = 'sept_{stereo_spacecraft}_{specie}_{sensor_pointing}_%Y_%j_{duration_of_average}_l2_v03.dat'

		url_pattern = base_pattern + filename_pattern

		file_scraper = Scraper(url_pattern, stereo_spacecraft = stereo_spacecraft, duration_of_average = dict_duration[duration_of_average] , 
									specie = dict_specie[specie], sensor_pointing = sensor_pointing)
		
		return file_scraper.filelist(timerange)


	def _makeimap(self):
		"""
		Helper Function:used to hold information about source.
		"""
		self.map_['source']     = 'university of kiel'
		self.map_['instrument'] = 'stereo/sept'
		self.map_['phyobs']     = 'electron intensities'
		self.map_['provider']   = 'ieap'


	@classmethod
	def _can_handle_query(cls, *query, **kwargs):
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
			if x.__class__.__name__ == 'Instrument' and x.value == 'stereo/sept':
				return all(chklist)
		return False
	


class HETClient(GenericClient):


	def _get_url_for_timerange(cls,timerange, stereo_spacecraft = 'ahead', duration_of_average = 15*u.min):
		"""
		Returns list of URLS to STEREO HET data files corresponding to value of input timerange.
		URL Source : http://www.srl.caltech.edu/STEREO/DATA/HET/
		
		The earliest data available is from Decemeber 2006. It takes a couple of months for latest data to be updated on 
		to the site.
		
		Parameters
		----------
		timerange: sunpy.time.TimeRange
			time range for which data is to be downloaded.
			Example value -  TimeRange('2006-12-01','2015-03-01')	
	        
		stereo_spacecraft: string	
			Default value - ahead
			Possible values - ahead, behind    # corresponding to spacecraft location
		
		duration_of_average: astropy units quantity
			Default value - 15*u.min
			Possible values - 1*u.min, 15*u.min, 1*u.h, 12*u.h, 1*u.d	
			#corresponding to duration over which data is averaged
		
		Returns
		-------
		urls : list
			list of URLs corresponding to the requested time range
		"""
	
		base_url = 'http://www.srl.caltech.edu/STEREO/DATA/HET/'
		
		possible_spacecraft = ['ahead', 'behind']
		possible_duration   = [1*u.min, 15*u.min, 1*u.h, 12*u.h, 1*u.d]
		
		dict_duration     = { 1*u.min :'1minute', 15*u.min : '15minute', 1*u.h : '1hour' , 12*u.h : '12hour', 1*u.d :'1day'}
		dict_time         = {1*u.min:'1m', 1*u.h :'1h', 15*u.min:'15m', 1*u.d :'1d', 12*u.h:'12h'}
		dict_spacecraft   = {'Ahead': 'AeH', 'Behind': 'BeH'}
		
		#Parameter Validations
		if timerange.start < datetime.datetime(2006,12,01):
			raise ValueError('Earliest date for which HET data is available is 2006-12-01')
		
		if stereo_spacecraft not in possible_spacecraft:
			raise ValueError('Possible stereo_spacecraft values: ' + ','.join(possible_spacecraft))
		
		if duration_of_average not in possible_duration:
			raise ValueError('Possible duration_of_average values as astropy unit quantities: ' + ','.join([str(i) for i in possible_duration]))
	
	
		stereo_spacecraft = stereo_spacecraft.capitalize()
		
		url_pattern =  base_url + '{stereo_spacecraft}/{duration_of_average}/{dict_spacecraft}%y%b.{dict_time}'
		
		file_scraper = Scraper(url_pattern, stereo_spacecraft = stereo_spacecraft, 
									duration_of_average = dict_duration[duration_of_average],
									dict_spacecraft = dict_spacecraft[stereo_spacecraft], 
									dict_time = dict_time[duration_of_average])
		
		return file_scraper.filelist(timerange)

	def _makeimap(self):
		"""
		Helper Function:used to hold information about source.
		
		"""
		self.map_['source']     = 'srl caltech'
		self.map_['instrument'] = 'stereo/het'
		self.map_['phyobs']     = 'electron flux'
		self.map_['provider']   = 'solar terrestrial relations observatory '

	@classmethod
	def _can_handle_query(cls, *query, **kwargs):
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
			if x.__class__.__name__ == 'Instrument' and x.value == 'stereo/het':
				return all(chklist)
		return False


class SITClient(GenericClient):

	def _get_url_for_timerange(cls,timerange, stereo_spacecraft = 'ahead', 
							   atomic_specie       = '4He', 
							   duration_of_average = 10*u.min):
		"""
		Returns list of URLS to STEREO SIT data files corresponding to value of input timerange.
		URL Source : http://www.srl.caltech.edu/STEREO/DATA/SIT/

		The earliest data available is from Jan 2007.

		Parameters
		----------
		timerange: sunpy.time.TimeRange
			time range for which data is to be downloaded.
			Example value -  TimeRange('2007-01-01','2015-03-01')	

		stereo_spacecraft: string	
			Default value - ahead
			Possible values - ahead, behind    # corresponding to spacecraft location

		atomic_specie:  string
			Default value - 4He
			Possible values - 4He, Fe, H, O

		duration_of_average: string
		Default value - 15min
			Possible values - 1min, 10min, 1hr, 1day		#corresponding to duration over which data is averaged

		Returns
		-------
		urls : list
			list of URLs corresponding to the requested time range

		"""

		possible_spacecraft = ['ahead', 'behind']
		possible_duration   = [1*u.min, 10*u.min, 1*u.h, 1*u.d]
		possible_specie     = ['4He', 'Fe', 'H','O']

		dict_duration     =   { 1*u.min :'1min', 10*u.min : '10min', 1*u.h : '1hr', 1*u.d :'1day'}
		dict_time         =   { 1*u.min :'%j' , 10*u.min : '%m',1*u.h : '',1*u.d : '' }

		#Parameter Validations
		if timerange.start < datetime.datetime(2007,01,01):
			raise ValueError('Earliest date for which SIT data is available is 2007-01-01')

		if stereo_spacecraft not in possible_spacecraft:
			raise ValueError('Possible stereo_spacecraft values: ' + ','.join(possible_spacecraft))

		if atomic_specie not in possible_specie:
			raise ValueError('Possible atomic_specie values: ' + ','.join(possible_specie))

		if duration_of_average not in possible_duration:
			raise ValueError('Possible duration_of_average values as astropy unit quantities: ' + ','.join([str(i) for i in possible_duration]))



		base_url = 'http://www.srl.caltech.edu/STEREO/DATA/SIT/{stereo_spacecraft[0]}/{duration_of_average}/'

		if duration_of_average in possible_duration[:2]:
			url_pattern = base_url + '{atomic_specie}/SIT_{stereo_spacecraft[1]}_{duration_of_average}_{atomic_specie}_%Y_{dict_time}.txt'
		else: 
			# duration_of_average in [1*u.h, 1*u.d ]:
			url_pattern = base_url + 'SIT_{stereo_spacecraft[1]}_{duration_of_average}_{atomic_specie}_%Y.txt'


		file_scraper = Scraper(url_pattern, stereo_spacecraft = [stereo_spacecraft, stereo_spacecraft.capitalize()],
									duration_of_average      = dict_duration[duration_of_average] , 
									atomic_specie            = atomic_specie, 
									dict_time                = dict_time[duration_of_average])

		return file_scraper.filelist(timerange)

	def _makeimap(self):
		"""
		Helper Function:used to hold information about source.
		"""
		self.map_['source']     = 'srl caltech'
		self.map_['instrument'] = 'stereo/sit'
		self.map_['phyobs']     = 'atomic intensities'
		self.map_['provider']   = 'solar terrestrial relations observatory '

	@classmethod
	def _can_handle_query(cls, *query, **kwargs):
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
			if x.__class__.__name__ == 'Instrument' and x.value == 'stereo/sit':
				return all(chklist)
		return False



class PLASTICClient(GenericClient):

	def _get_url_for_timerange(cls,timerange, stereo_spacecraft = 'ahead', 
							   			duration_of_average = 10*u.min):
		"""
		Returns list of URLS to STEREO PLASTIC data files corresponding to value of input timerange.
		URL Source : http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/

		The earliest data available is from 2007.

		Parameters
		----------
		timerange: sunpy.time.TimeRange
			time range for which data is to be downloaded.
			Example value -  TimeRange('2007-02-14','2014-12-17')	

		stereo_spacecraft: string	
			Default value - ahead
			Possible values - ahead, behind    # corresponding to spacecraft location


		duration_of_average: string
		Default value - 10*u.min
			Possible values - 1*u.min, 10*u.min, 1*u.h		
			#corresponding to duration over which data is averaged

		Returns
		-------
		urls : list
			list of URLs corresponding to the requested time range

		"""

		possible_spacecraft = ['ahead', 'behind']
		possible_duration   = [1*u.min, 10*u.min, 1*u.h ]

		dict_spacecraft = { 'ahead': 'A', 'behind': 'B' }
		dict_duration   = { 1*u.min: '1min', 10*u.min: '10min', 1*u.h: '1hr'}


		#Parameter Validations
		if timerange.start < datetime.datetime(2007,02,14):
			raise ValueError('Earliest date for which PLASTIC data is available is 2007-02-14')

		if stereo_spacecraft not in possible_spacecraft:
			raise ValueError('Possible stereo_spacecraft values: ' + ','.join(possible_spacecraft))

		if duration_of_average not in possible_duration:
			raise ValueError('Possible duration_of_average values as astropy unit quantities: ' + ','.join([str(i) for i in possible_duration]))

		base_url = 'http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/'

		url_pattern =  base_url + '{duration_of_average}/{stereo_spacecraft}/%Y/ST{stereo_spacecraft}_L2_PLA_1DMax_{duration_of_average}_%Y%m%d_%j_V09.txt'

		file_scraper = Scraper(url_pattern, stereo_spacecraft   = dict_spacecraft[stereo_spacecraft], 
										   duration_of_average = dict_duration[duration_of_average])
		
		return file_scraper.filelist(timerange)

	def _makeimap(self):
		"""
		Helper Function:used to hold information about source.
		"""
		self.map_['source']     = 'nasa'
		self.map_['instrument'] = 'stereo/plastic'
		self.map_['phyobs']     = '1d maxwellian proton moments'
		self.map_['provider']   = 'stereo science center'

	@classmethod
	def _can_handle_query(cls, *query, **kwargs):
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
			if x.__class__.__name__ == 'Instrument' and x.value == 'stereo/plastic':
				return all(chklist)
		return False


class MAGClient(GenericClient):

	def _get_url_for_timerange(cls,timerange, stereo_spacecraft = 'ahead'):
		"""
		Returns list of URLS to STEREO MagPlasma data files corresponding to value of input timerange.
		URL Source : http://stereo.ssl.berkeley.edu/L2_data.php

		The earliest data available is from 2006.

		Parameters
		----------
		timerange: sunpy.time.TimeRange
			time range for which data is to be downloaded.

		stereo_spacecraft: string	
			Default value - ahead
			Possible values - ahead, behind    # corresponding to spacecraft location

		Returns
		-------
		urls : list
			list of URLs corresponding to the requested time range

		"""

		possible_spacecraft = ['ahead', 'behind']
		dict_spacecraft = {'ahead':'A', 'behind':'B'}

		#Parameter Validations
		if timerange.start < datetime.datetime(2006,01,01):
			raise ValueError('Earliest date for which PLASTIC data is available is 2006-01-01')

		if stereo_spacecraft not in possible_spacecraft:
			raise ValueError('Possible stereo_spacecraft values: ' + ','.join(possible_spacecraft))


		url_pattern = ('http://stereo.ssl.berkeley.edu/l2data/{stereo_spacecraft}/magplasma/ST{dict_spacecraft}_L2_MAGPLASMA_1m_%Y_V01.cdf')
		file_scraper = Scraper(url_pattern, stereo_spacecraft = stereo_spacecraft, dict_spacecraft = dict_spacecraft[stereo_spacecraft])
		
		return file_scraper.filelist(timerange)

	def _makeimap(self):
		"""
		Helper Function:used to hold information about source.
		"""
		self.map_['source']     = 'ssl berkeley'
		self.map_['instrument'] = 'stereo/mag'
		self.map_['phyobs']     = 'vector magnetic field'
		self.map_['provider']   = 'solar terrestrial relations observatory'

	@classmethod
	def _can_handle_query(cls, *query, **kwargs):
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
			if x.__class__.__name__ == 'Instrument' and x.value == 'stereo/mag':
				return all(chklist)
		return False

	
class LETClient(GenericClient):

	def _get_url_for_timerange(cls,timerange, duration_of_average, type_of_data , specie, stereo_spacecraft = 'ahead'):
		"""
		Returns list of URLS to STEREO LET data files corresponding to value of input timerange.
		URL Source : http://www.srl.caltech.edu/STEREO/Public/LET_public.html

		The earliest data available is from 13-Nov-2006.

		Parameters
		----------
		timerange: sunpy.time.TimeRange
			time range for which data is to be downloaded.
			Example value -  TimeRange('2007-01-01','2008-06-01')	

		duration_of_average: string
			Default value - 10*u.min
			Possible values - 1*u.min, 10*u.min, 1*u.h, 1*u.d, 27*u.d	
			#corresponding to duration over which data is averaged

		type_of_data:  string
			Possible values - depends on other parameters
			if duration_of_average = 27*u.d:
				#Possible Values: summed, narrow
			else:
				#Possible values: sectored, standard, summed

		specie:  string
			Possible values - depends on other parameters
			if type_of_data = 'Sectored' and duration_of_average in [ 1*u.min, 10*u.min, 1*u.h, 1*u.d]:
				#Possible values: CNO_hi,CNO_lo, Fe_hi, Fe_lo, H_lo, He3_lo, He4_hi, He4_lo, He_lo, NeMgSi_hi, NeMgSi_lo
			else:
				#Possible values: Al, Ar, C, Ca, Fe, H, He, He3, He4, Mg, N, Na, Ne, Ni, O, S, Si

		stereo_spacecraft: string	
			Default value - ahead
			Possible values - ahead, behind    # corresponding to spacecraft location

		
		Returns
		-------
		urls : list
			list of URLs corresponding to the requested time range

		"""
		possible_spacecraft = ['ahead', 'behind']
		possible_duration   = [1*u.min, 10*u.min, 1*u.h, 1*u.d, 27*u.d]

		dict_duration       = { 1*u.min : '1Minute', 10*u.min : '10Minute', 1*u.h: 'Hourly', 1*u.d: 'Daily', 27*u.d: '27day'}


		#Parameter Validations
		if timerange.start < datetime.datetime(2006,11,13):
			raise ValueError('Earliest date for which SEPT data is available is 2006-11-13')

		if stereo_spacecraft not in possible_spacecraft:
			raise ValueError('Possible stereo_spacecraft values: ' + ','.join(possible_spacecraft))

		if duration_of_average not in possible_duration:
			raise ValueError('Possible duration_of_average values: ' + ','.join([str(i) for i in possible_duration]))


		url_base_pattern = 'http://www.srl.caltech.edu/STEREO/DATA/Level1/Public/{stereo_spacecraft}/{duration_of_average}/'
		duration_of_average = dict_duration[duration_of_average]

		#Type of Data validation on the basis of duration of Average 
		if duration_of_average == '27day':
			possible_datatype = ['summed','narrow']
		else: 
			possible_datatype = ['sectored','standard', 'summed']
			url_base_pattern =  url_base_pattern + '%Y/{type_of_data}/'

		if type_of_data not in possible_datatype:
				raise ValueError('Possible type_of_data values for {duration} average: '.format(duration = str(duration_of_average)) + ','.join(possible_datatype))

		if duration_of_average != '27day':
			type_of_data = type_of_data.capitalize()

		#Specie validation on the basis of type of data and duration of average
		if type_of_data == 'Sectored' and duration_of_average in ['1Minute', '10Minute', 'Hourly', 'Daily']:
			possible_specie = ['CNO_hi', 'CNO_lo', 'Fe_hi', 'Fe_lo', 'H_lo', 'He3_lo', 'He4_hi', 'He4_lo', 'He_lo', 'NeMgSi_hi', 'NeMgSi_lo']
		else:
			possible_specie = ['Al','Ar','C','Ca','Fe', 'H','He', 'He3','He4','Mg','N','Na','Ne','Ni','O','S','Si']

		if specie not in possible_specie:
				raise ValueError('Invalid Specie Selection. Possible values are: ' + ','.join(possible_specie))


		if duration_of_average in ['1Minute', '10Minute']:
			url_base_pattern =  url_base_pattern + '{specie}/'


		duration_replacement = {'1Minute': '', '10Minute':'_10min', 'Hourly':'_1hr', 'Daily':'_1day', '27day':''}
		type_replacement = {'Sectored': '_sectored', 'Standard':'', 'Summed': '_summed', 'summed':'_summed', 'narrow':'_narrow'}
		date_replacement = {'1Minute':'_%j','10Minute':'_%m', 'Hourly':'', 'Daily':'','27day':''}


		#Adding file name to the final directory pattern
		if duration_of_average != '27day':
			url_pattern =  url_base_pattern + '{specie}{type_replacement}_{stereo_spacecraft}_%Y{date_replacement}{duration_replacement}_level1_11.txt'

			file_scraper = Scraper(url_pattern, stereo_spacecraft = stereo_spacecraft, 
										duration_of_average = duration_of_average , 
							 			type_of_data = type_of_data, 
							 			specie = specie, 
							 			duration_replacement = duration_replacement[duration_of_average],
							 			type_replacement = type_replacement[type_of_data], 
							 			date_replacement = date_replacement[duration_of_average])

			return file_scraper.filelist(timerange)

		else:
			url_pattern = (url_base_pattern + '{specie}{type_replacement}_{stereo_spacecraft}.txt').format(stereo_spacecraft = stereo_spacecraft,
																					duration_of_average = duration_of_average, specie = specie,
																					type_replacement = type_replacement[type_of_data])
			return url_pattern


	def _makeimap(self):
		"""
		Helper Function:used to hold information about source.
		"""
		self.map_['source']     = 'srl caltech'
		self.map_['instrument'] = 'stereo/let'
		self.map_['phyobs']     = 'atomic intensities'
		self.map_['provider']   = 'solar terrestrial relations observatory '

	@classmethod
	def _can_handle_query(cls, *query, **kwargs):
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
			if x.__class__.__name__ == 'Instrument' and x.value == 'stereo/let':
				return all(chklist)
		return False


