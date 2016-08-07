import datetime

from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange


def _get_url_for_timerange(timerange = TimeRange('2007-02-14','2014-12-17'), stereo_spacecraft = 'A', 
						   duration_of_average = '10min'):
	"""
	Returns list of URLS to STEREO PLASTIC data files corresponding to value of input timerange.
	URL Source : http://www2.physik.uni-kiel.de/stereo/data/sept/level2/

	The earliest data available is from 2007.

	Parameters
	----------
	timerange: sunpy.time.TimeRange
	    time range for which data is to be downloaded.
	    Default value -  TimeRange('2007-02-14','2014-12-17')	

	stereo_spacecraft: string	
		Default value - A
		Possible values - A, B    # corresponding to spacecraft location


	duration_of_average: string
	Default value - 10min
		Possible values - 1min, 10min, 1hr		
		#corresponding to duration over which data is averaged

	Returns
	-------
	urls : list
	    list of URLs corresponding to the requested time range

	"""

	#Parameter Validations
	if timerange.start < datetime.datetime(2007,02,14):
		raise ValueError('Earliest date for which PLASTIC data is available is 2007-02-14')

	if stereo_spacecraft not in ['A','B']:
		raise ValueError('Possible stereo_spacecraft values: \'A\', \'B\'')

	if duration_of_average not in ['1min', '10min', '1hr']:
		raise ValueError('Possible duration_of_average values: \'1min\', \'10min\',\'1hr\' ')


	solmon_pattern = ('http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/{duration_of_average}/{stereo_spacecraft}/%Y/ST{stereo_spacecraft}_L2_PLA_1DMax_{duration_of_average}_%Y%m%d_%j_V09.txt')
	solmon = Scraper(solmon_pattern, stereo_spacecraft = stereo_spacecraft, duration_of_average = duration_of_average)
	
	return solmon.filelist(timerange)

"""
print _get_url_for_timerange(TimeRange('2008-03-01','2008-06-01'), stereo_spacecraft = 'B')

"""

