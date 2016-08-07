import datetime

from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange


def _get_url_for_timerange(timerange = TimeRange('2006-12-01','2015-03-01'), stereo_spacecraft = 'Ahead', duration_of_average = '15minute'):
    """
    Returns list of URLS to STEREO HET data files corresponding to value of input timerange.
    URL Source : http://www.srl.caltech.edu/STEREO/DATA/HET/

    The earliest data available is from Decemeber 2006. It takes a couple of months for latest data to be updated on 
    to the site.

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Default value -  TimeRange('2006-12-01','2015-03-01')	
    
    stereo_spacecraft: string	
    	Default value - Ahead
    	Possible values - Ahead, Behind    # corresponding to spacecraft location

    duration_of_average: string
    Default value - 15minute
    	Possible values - 1minute, 15minute, 1hour, 12hour, 1day		#corresponding to duration over which data is averaged

    Returns
    -------
    urls : list
        list of URLs corresponding to the requested time range
    """
    #Parameter Validations
    if timerange.start < datetime.datetime(2006,12,01):
    	raise ValueError('Earliest date for which HET data is available is 2006-12-01')

    if stereo_spacecraft not in ['Ahead','Behind']:
    	raise ValueError('Possible stereo_spacecraft values: \'Ahead\', \'Behind\'')

    if duration_of_average not in ['1minute', '15minute', '1hour', '12hour', '1day']:
    	raise ValueError('Possible duration_of_average values: \'1minute\', \'15minute\',\'1hour\',\'12hour\',\'1day\'')



    dic_time = {'1minute':'1m', '1hour':'1h', '15minute':'15m', '1day':'1d', '12hour':'12h'}
    dic_spacecraft = {'Ahead': 'AeH', 'Behind': 'BeH'}

    solmon_pattern = ('http://www.srl.caltech.edu/STEREO/DATA/HET/{stereo_spacecraft}/{duration_of_average}/{dic_spacecraft}%y%b.{dic_time}')
    solmon = Scraper(solmon_pattern, stereo_spacecraft = stereo_spacecraft, duration_of_average = duration_of_average,
    					dic_spacecraft = dic_spacecraft[stereo_spacecraft], dic_time = dic_time[duration_of_average])

    return solmon.filelist(timerange)

"""
print _get_url_for_timerange(TimeRange('2007-03-01','2007-06-01'))
print _get_url_for_timerange(TimeRange('2007-03-01','2007-06-01'), stereo_spacecraft = 'Behind' , duration_of_average = '1day')
print _get_url_for_timerange(TimeRange('2007-03-01','2007-06-01'), duration_of_average = '12hour')

"""
