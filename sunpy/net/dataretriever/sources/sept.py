from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange
import datetime

def _get_url_for_timerange(timerange = TimeRange('2007-01-20','2015-01-01'), stereo_spacecraft = 'ahead', 
						   duration_of_average = '10min', specie = 'ele', sensor_pointing = 'asun'):
	"""
	Returns list of URLS to STEREO SEPT data files corresponding to value of input timerange.
	URL Source : http://www2.physik.uni-kiel.de/stereo/data/sept/level2/

	The earliest data available is from 20-Jan-2007.

	Parameters
	----------
	timerange: sunpy.time.TimeRange
	    time range for which data is to be downloaded.
	    Default value -  TimeRange('2007-01-20','2015-01-01')	

	stereo_spacecraft: string	
		Default value - ahead
		Possible values - ahead, behind    # corresponding to spacecraft location

	specie:  string
		Default value - ele
		Possible values - ele, ion

	duration_of_average: string
	Default value - 10min
		Possible values - 1min, 10min, 1h, 1d		
		#corresponding to duration over which data is averaged

	sensor_pointing: string
		Default value - asun
		Possible values - asun, sun, north, south, omni

	Returns
	-------
	urls : list
	    list of URLs corresponding to the requested time range

	"""

	#Parameter Validations
	if timerange.start < datetime.datetime(2007,01,20):
		raise ValueError('Earliest date for which SEPT data is available is 2007-01-20')

	if stereo_spacecraft not in ['ahead','behind']:
		raise ValueError('Possible stereo_spacecraft values: \'ahead\', \'behind\'')

	if specie not in ['ele', 'ion']:
		raise ValueError('Possible specie values: \'ele\', \'ion\'')

	if duration_of_average not in ['1min', '10min', '1hr','1day']:
		raise ValueError('Possible duration_of_average values: \'1min\', \'10min\',\'1hr\',\'1day\'')

	if sensor_pointing not in ['asun', 'sun', 'north','south','omni']:
		raise ValueError('Possible sensor_pointing values: \'asun\', \'sun\',\'north\',\'south\', \'omni\' ')


	solmon_pattern = 'http://www2.physik.uni-kiel.de/stereo/data/sept/level2/{stereo_spacecraft}/{duration_of_average}/%Y/sept_{stereo_spacecraft}_{specie}_{sensor_pointing}_%Y_%j_{duration_of_average}_l2_v03.dat'

	solmon = Scraper(solmon_pattern, stereo_spacecraft = stereo_spacecraft, duration_of_average = duration_of_average , 
						specie = specie, sensor_pointing = sensor_pointing)
	
	return solmon.filelist(timerange)


"""
print _get_url_for_timerange(TimeRange('2007-02-01','2007-03-01'),stereo_spacecraft = 'behind')

"""



