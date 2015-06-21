import datetime

from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange


def _get_url_for_timerange(timerange = TimeRange('2007-01-01','2015-03-01'), stereo_spacecraft = 'ahead', 
						   atomic_specie = '4He', duration_of_average = '10min'):
	"""
	Returns list of URLS to STEREO SIT data files corresponding to value of input timerange.
	URL Source : http://www.srl.caltech.edu/STEREO/DATA/SIT/

	The earliest data available is from Jan 2007.

	Parameters
	----------
	timerange: sunpy.time.TimeRange
	    time range for which data is to be downloaded.
	    Default value -  TimeRange('2007-01-01','2015-03-01')	

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

	#Parameter Validations
	if timerange.start < datetime.datetime(2007,01,01):
		raise ValueError('Earliest date for which SIT data is available is 2007-01-01')

	if stereo_spacecraft not in ['ahead','behind']:
		raise ValueError('Possible stereo_spacecraft values: \'ahead\', \'behind\'')

	if atomic_specie not in ['4He', 'Fe', 'H','O']:
		raise ValueError('Possible atomic_specie values: \'4He\', \'Fe\',\'H\',\'O\'')

	if duration_of_average not in ['1min', '10min', '1hr','1day']:
		raise ValueError('Possible duration_of_average values: \'1min\', \'10min\',\'1hr\',\'1day\'')


	dict_time = {'1min':'%j','10min':'%m','1hr':'','1day':''}


	base_url = 'http://www.srl.caltech.edu/STEREO/DATA/SIT/{stereo_spacecraft[0]}/{duration_of_average}/'
	if duration_of_average in ['1min','10min']:
		solmon_pattern = base_url + '{atomic_specie}/SIT_{stereo_spacecraft[1]}_{duration_of_average}_{atomic_specie}_%Y_{dict_time}.txt'
	else: 
		# duration_of_average in ['1hr','1day']:
		solmon_pattern = base_url + 'SIT_{stereo_spacecraft[1]}_{duration_of_average}_{atomic_specie}_%Y.txt'


	solmon = Scraper(solmon_pattern, stereo_spacecraft = [stereo_spacecraft,stereo_spacecraft.capitalize()],
						duration_of_average = duration_of_average , 
						atomic_specie = atomic_specie, 
						dict_time = dict_time[duration_of_average])

	return solmon.filelist(timerange)

"""
print _get_url_for_timerange(TimeRange('2007-02-04','2007-03-01')) 
print _get_url_for_timerange(TimeRange('2007-02-04','2008-03-01'), stereo_spacecraft = 'behind', atomic_specie ='O', duration_of_average = '1day') 

"""


