from datetime import timedelta,datetime
import calendar
from calendar import monthrange
from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time, TimeRange
from dateutil.relativedelta import relativedelta


def _get_url_for_date_range(*args):
    """Returns a list of URLs to the SIT data for the specified range of dates.
    Parameters
    ----------
    args : TimeRange, datetimes, date strings
        Date range should be specified using a TimeRange, or start
        and end dates at datetime instances or date strings.

    instrument_orientation : string
        SIT instrument orientation (default = "ahead")
    atomic_specie : string 
        atomic specie to consider ( default = "4He")
    type_of_average : string
        average type to retrieve (default = "1min")

    """
    # TimeRange
    if len(args) == 4 and isinstance(args[0], TimeRange):
        start = args[0].start
        end = args[0].end
        instrument_orientation = args[1]
        type_of_average = args[2]
        atomic_specie = args[3]
    elif len(args) == 5:
        #start, end dates given
        start = parse_time(args[0])
        end = parse_time(args[1])
        instrument_orientation = args[2]
        type_of_average = args[3]
        atomic_specie = args[4]

    if end < start:
        raise ValueError('start time > end time')

    # find out base url of data to query from the instrument orientation, average type, atomic specie 
    base_url = 'http://www.srl.caltech.edu/STEREO/DATA/SIT/' + instrument_orientation + '/' + type_of_average + '/' 

    if type_of_average == '1min' or type_of_average == '10min':
        base_url = base_url + atomic_specie + '/'
    
    #adding the file name to base url
    if instrument_orientation == "ahead":
        base_url = base_url + 'SIT_Ahead_' + type_of_average + '_' + atomic_specie + '_'
    elif instrument_orientation == "behind":
        base_url = base_url + 'SIT_Behind_' + type_of_average + '_' + atomic_specie + '_'

    #Date Generator to generate dates in between two dates. Inclusive of both end and start dates.         
    def daterange(start_date, end_date, delta = 'day'):
        """
        start_date: datetime.datetime 
            start date of TimeRange 

        end_date: datetime.datetime 
            end date of TimeRange 

        delta : string
            type of distribution in data source files (default = 'day')
            POSSIBLE VALUES: 'day', 'month', 'year'
        """
        def monthdelta(d1, d2):
            delta = 0
            while True:
                mdays = monthrange(d1.year, d1.month)[1]
                d1 += timedelta(days=mdays)
                if d1 <= d2:
                    delta += 1
                else:
                    break
            return delta

        def add_months(sourcedate,months):
            """ Adds the number of months specified in 'months' argument to the 
            sourcedate specified as the first argument"""

            month = sourcedate.month - 1 + months
            year = sourcedate.year + month / 12
            month = month % 12 + 1
            day = min(sourcedate.day,calendar.monthrange(year,month)[1])
            return datetime(year,month,day)

        def add_years(d, years):
            """Return a date that's `years` years after the date (or datetime)
            object `d`. Return the same calendar date (month and day) in the
            destination year, if it exists, otherwise use the following day
            (thus changing February 29 to March 1).

            """
            try:
                return d.replace(year = d.year + years)
            except ValueError:
                return d + (date(d.year + years, 1, 1) - date(d.year, 1, 1))

        if delta == 'day':
            for n in range((end_date - start_date).days +1):
                yield start_date + timedelta(n)
        elif delta == 'month':
            for n in range( monthdelta(start_date, end_date) + 1):
                yield  add_months(start_date,n )
        elif delta == 'year':  
            for n in range(relativedelta(end_date, start_date).years +1):
                yield add_years(start_date,n)      





    #creating a list of links to file for each date in TimeRange by looping and adding date of file to base url
    url = []
    url_meta = (start,end)
    
    
    if type_of_average == '1hr' or type_of_average == '1day': 
        for single_date in daterange(start, end, delta = 'year'):  
            url.append(base_url + "{date:%Y}.txt".format(date=single_date))
    elif type_of_average == '1min':
        for single_date in daterange(start, end, delta = 'day'):
            day_of_year = (str)(single_date.timetuple().tm_yday)
            for i in [0,1,2]:
                if len(day_of_year) < 3:
                    day_of_year = '0' + day_of_year
            url.append(base_url + "{date:%Y}_".format(date = single_date) + day_of_year + ".txt".format(
                date=single_date))
    elif type_of_average == '10min':
        for single_date in daterange(start, end, delta = 'month'):
            url.append(base_url + "{date:%Y}_{date:%m}.txt".format(
                date=single_date))

    return url, url_meta

"""

#Test Cases for above routine (will add separate test case in pytest for the lightcurve later when I complete it)

print _get_url_for_date_range('2009-10-02', '2009-10-04','ahead','1min','4He')
print _get_url_for_date_range('2008-06-02', '2010-02-04','behind','10min','Fe')
print _get_url_for_date_range('2009-10-02', '2009-10-04','ahead','1hr','O')
print _get_url_for_date_range('2009-10-02', '2009-10-04','behind','1day','H')

print _get_url_for_date_range(TimeRange('2012/06/01', '2012/06/05'),'ahead','1min','4He')
print _get_url_for_date_range(TimeRange('2008-06-01', '2010-06-05'),'behind','10min','Fe')
print _get_url_for_date_range(TimeRange('2012/06/01', '2012/06/05'),'ahead','1hr','O')
print _get_url_for_date_range(TimeRange('2012/06/01', '2012/06/05'),'behind','1day','H')


"""





