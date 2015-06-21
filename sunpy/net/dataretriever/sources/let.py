import datetime

from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange

def _get_url_for_timerange(timerange, duration_of_average, type_of_data , specie, stereo_spacecraft = 'ahead'):
	"""
	Returns list of URLS to STEREO LET data files corresponding to value of input timerange.
	URL Source : http://www.srl.caltech.edu/STEREO/Public/LET_public.html

	The earliest data available is from 13-Nov-2006.

	Parameters
	----------
	timerange: sunpy.time.TimeRange
	    time range for which data is to be downloaded.
	    Default value -  TimeRange('2007-01-01','2008-06-01')	

	duration_of_average: string
		Default value - 10minute
		Possible values - 1Minute, 10Minute, Hourly, Daily, 27day	
		#corresponding to duration over which data is averaged

	stereo_spacecraft: string	
		Default value - ahead
		Possible values - ahead, behind    # corresponding to spacecraft location

	type_of_data:  string
		Possible values - depends on other parameters
		if duration_of_average = '27day':
			#Possible Values: summed, narrow
		else:
			#Possible values: Sectored, Standard, Summed

	specie:  string
		Possible values - depends on other parameters
		if type_of_data = 'Sectored' and duration_of_average in ['1Minute', '10Minute', 'Hourly', 'Daily']:
			#Possible values: CNO_hi,CNO_lo, Fe_hi, Fe_lo, H_lo, He3_lo, He4_hi, He4_lo, He_lo, NeMgSi_hi, NeMgSi_lo
		else:
			#Possible values: Al, Ar, C, Ca, Fe, H, He, He3, He4, Mg, N, Na, Ne, Ni, O, S, Si


	Returns
	-------
	urls : list
	    list of URLs corresponding to the requested time range

	"""

	#Parameter Validations
	if timerange.start < datetime.datetime(2006,11,13):
		raise ValueError('Earliest date for which SEPT data is available is 2006-11-13')

	if stereo_spacecraft not in ['ahead','behind']:
		raise ValueError('Possible stereo_spacecraft values: \'ahead\', \'behind\'')

	if duration_of_average not in ['1Minute', '10Minute', 'Hourly','Daily','27day']:
		raise ValueError('Possible duration_of_average values: \'1Minute\', \'10Minute\', \'Hourly\', \'Daily\', \'27day\' ')


	solmon_base_pattern = ('http://www.srl.caltech.edu/STEREO/DATA/Level1/Public/{stereo_spacecraft}/{duration_of_average}/')


	#Type of Data validation on the basis of duration of Average 
	if duration_of_average == '27day':
		if type_of_data not in ['summed','narrow']:
			raise ValueError('Possible type_of_data values for 27day average : \'summed\', \'narrow\'')
			#Possible Values: summed, narrow
	else: 
		if type_of_data not in ['Sectored','Standard', 'Summed']:
			raise ValueError('Possible type_of_data values: \'Sectored\', \'Standard\', \'Summed\' ')
			#Possible values: Sectored, Standard, Summed
		solmon_base_pattern =  solmon_base_pattern + '%Y/{type_of_data}/'


	#Specie validation on the basis of type of data and duration of average
	if type_of_data == 'Sectored' and duration_of_average in ['1Minute', '10Minute', 'Hourly', 'Daily']:
		#Possible values: CNO_hi,CNO_lo, Fe_hi, Fe_lo, H_lo, He3_lo, He4_hi, He4_lo, He_lo, NeMgSi_hi, NeMgSi_lo
		if specie not in ['CNO_hi', 'CNO_lo', 'Fe_hi', 'Fe_lo', 'H_lo', 'He3_lo', 'He4_hi', 'He4_lo', 'He_lo', 'NeMgSi_hi', 'NeMgSi_lo']:
			raise ValueError('Invalid Specie Selection. Possible values: CNO_hi,CNO_lo, Fe_hi, Fe_lo, H_lo, He3_lo, He4_hi, He4_lo, He_lo, NeMgSi_hi, NeMgSi_lo')

	else:
		#Possible values: Al, Ar, C, Ca, Fe, H, He, He3, He4, Mg, N, Na, Ne, Ni, O, S, Si
		if specie not in ['Al','Ar','C','Ca','Fe', 'H','He', 'He3','He4','Mg','N','Na','Ne','Ni','O','S','Si']:
			raise ValueError('Invalid Specie Selection. Possible values: Al, Ar, C, Ca, Fe, H, He, He3, He4, Mg, N, Na, Ne, Ni, O, S, Si')



	if duration_of_average in ['1Minute', '10Minute']:
		solmon_base_pattern =  solmon_base_pattern + '{specie}/'


	duration_replacement = {'1Minute': '', '10Minute':'_10min', 'Hourly':'_1hr', 'Daily':'_1day', '27day':''}
	type_replacement = {'Sectored': '_sectored', 'Standard':'', 'Summed': '_summed', 'summed':'_summed', 'narrow':'_narrow'}
	date_replacement = {'1Minute':'_%j','10Minute':'_%m', 'Hourly':'', 'Daily':'','27day':''}


	#Adding file name to the final directory pattern
	if duration_of_average != '27day':
		solmon_pattern =  solmon_base_pattern + '{specie}{type_replacement}_{stereo_spacecraft}_%Y{date_replacement}{duration_replacement}_level1_11.txt'

		solmon = Scraper(solmon_pattern, stereo_spacecraft = stereo_spacecraft, duration_of_average = duration_of_average , 
					 type_of_data = type_of_data, specie = specie, duration_replacement = duration_replacement[duration_of_average],
					 type_replacement = type_replacement[type_of_data], date_replacement = date_replacement[duration_of_average])

		return solmon.filelist(timerange)

	else:
		solmon_pattern =  (solmon_base_pattern + '{specie}{type_replacement}_{stereo_spacecraft}.txt').format(stereo_spacecraft = stereo_spacecraft,
																					duration_of_average = duration_of_average, specie = specie,
																					type_replacement = type_replacement[type_of_data])
		return solmon_pattern



print _get_url_for_timerange(TimeRange('2007-01-01','2008-06-01'), stereo_spacecraft = 'ahead', duration_of_average = '10Minute', type_of_data = 'Standard', specie = 'Al')

