# -*- coding: utf-8 -*-
"""Provides programs to process and analyze STEREO SIT data."""

# This module was developed with funding from the Google Summer of Code 2015

from __future__ import absolute_import

__authors__ = ["Ankit Kumar"]
__email__ = "ankitkmr.iitk@gmail.com"

from datetime import timedelta,datetime
from calendar import monthrange
from sunpy.lightcurve import LightCurve
from dateutil.relativedelta import relativedelta
from sunpy.time import parse_time, TimeRange, is_time_in_given_format
from sunpy.util import net
from shlex import split
from astropy.io import ascii
from astropy.table import Table, Column

__all__ = ['SITLightCurve']


class SITLightCurve(LightCurve):
    """
    SIT LightCurve. Provides SIT data back to 2007-07.
    Most recent data is usually available one or two days late.
    args : type_of_average, stereo_spacecraft, atomic_specie
    POSSIBLE KEYWORD VALUES: (value has to be one of these only for each keyword)
    (Also the values need to be passed in this order only for the positional arguments)
    stereo_spacecraft : ahead, behind
    type_of_average        : 1min, 10min, 1hr, 1day
    atomic_specie          : 4He, Fe, H, O
    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> from sunpy.time import TimeRange
    >>> sit = lc.SITLightCurve.create(TimeRange('2012/06/01', '2012/06/05'), "ahead", "1min", "4He")
    
    References
    ----------
    | http://www.srl.caltech.edu/STEREO/Public/SIT_public.html
    """
    
    # def peek(self, title="GOES Xray Flux"):
    #     """Plots GOES light curve is the usual manner"""
    #     figure = plt.figure()
    #     axes = plt.gca()

    #     dates = matplotlib.dates.date2num(parse_time(self.data.index))

    #     axes.plot_date(dates, self.data['xrsa'], '-',
    #                  label='0.5--4.0 $\AA$', color='blue', lw=2)
    #     axes.plot_date(dates, self.data['xrsb'], '-',
    #                  label='1.0--8.0 $\AA$', color='red', lw=2)

    #     axes.set_yscale("log")
    #     axes.set_ylim(1e-9, 1e-2)
    #     axes.set_title(title)
    #     axes.set_ylabel('Watts m$^{-2}$')
    #     axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])

    #     ax2 = axes.twinx()
    #     ax2.set_yscale("log")
    #     ax2.set_ylim(1e-9, 1e-2)
    #     ax2.set_yticks((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
    #     ax2.set_yticklabels((' ', 'A', 'B', 'C', 'M', 'X', ' '))

    #     axes.yaxis.grid(True, 'major')
    #     axes.xaxis.grid(False, 'major')
    #     axes.legend()

    #     # @todo: display better tick labels for date range (e.g. 06/01 - 06/05)
    #     formatter = matplotlib.dates.DateFormatter('%H:%M')
    #     axes.xaxis.set_major_formatter(formatter)

    #     axes.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
    #     figure.autofmt_xdate()
    #     figure.show()

    #     return figure

    @classmethod
    def _get_default_url(cls):
        """
           Retrieve latest SIT data for 1 minute average of 4He specie from ahead spacecraft 
           ( if no other data is specified )
        """
        now = datetime.datetime.utcnow()
        time_range = TimeRange(datetime.datetime(now.year, now.month, now.day), now)
        url_does_exist = net.url_exists(cls._get_url_for_date_range(time_range))
        while not url_does_exist:
            time_range = TimeRange(time_range.start-datetime.timedelta(days=1),
                                   time_range.start)
            url_does_exist = net.url_exists(cls._get_url_for_date_range(time_range))
        return cls._get_url_for_date_range(time_range, "ahead", '1min', "4He")

    

    @staticmethod
    def _get_url_for_date_range(*args):
        """Returns a list of URLs to the SIT data for the specified range of dates along with the 
         tuple containing specified start and end dates .
        Parameters
        ----------
        args : TimeRange, datetimes, date strings
            Date range should be specified using a TimeRange, or start
            and end dates at datetime instances or date strings.
        stereo_spacecraft : string
            SIT stereo spacecraft (default = "ahead")
        atomic_specie : string 
            atomic specie to consider ( default = "4He")
        type_of_average : string
            average type to retrieve (default = "1min")
        """
        # # If time duration given as sunpy.time.TimeRange
        if len(args) == 4 and isinstance(args[0], TimeRange):
            start = args[0].start
            end = args[0].end
            stereo_spacecraft = args[1]
            type_of_average = args[2]
            atomic_specie = args[3]

        elif len(args) == 5:
            #If start and end dates for observation period is given
            start = parse_time(args[0])
            end = parse_time(args[1])
            stereo_spacecraft = args[2]
            type_of_average = args[3]
            atomic_specie = args[4]
            
        else: 
            raise ValueError('must recieve 4 or 5 arguments only')
            # 4 arguments incase of first argument being TimeRange and 5 in case of it being start and end dates
            # The other three arguments : stereo_spacecraft, type_of_average, atomic_specie remain same for either

        if end < start:
            raise ValueError('start time > end time')


        # find out base url of data to query from the stereo spacecraft, average type and atomic specie quried for
        base_url = 'http://www.srl.caltech.edu/STEREO/DATA/SIT/' + stereo_spacecraft + '/' + type_of_average + '/' 

        if type_of_average == '1min' or type_of_average == '10min':
            base_url = base_url + atomic_specie + '/'
        
        #adding the file name to base url
        base_url = base_url + 'SIT_' + stereo_spacecraft.capitalize() + '_' + type_of_average + '_' + atomic_specie + '_'

        #Date Generator to generate dates in between the start and end dates. Inclusive of both end and start dates.         
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
                """ 
                Returns the number of months, while moving forward from start date, in between start date and end date(included).
                d1 : datetime.datetime
                    start_date
                d2 : datetime.datetime
                    end_date
                >>>monthdelta(parse_time('2010-09-04'),parse_time('2012-04-05'))
                19
                >>>monthdelta(parse_time('2010-09-04'),parse_time('2011-04-05'))
                7
                >>>monthdelta(parse_time('2011-05-04'),parse_time('2011-09-05'))
                4
                """
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
                day = min(sourcedate.day,monthrange(year,month)[1])
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



    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a STEREO SIT file from
        http://www.srl.caltech.edu/STEREO/Public/SIT_public.html
        
        and returns header and astropy.Table object containing data
    
        """
        header = []
        data = ascii.read(filepath, delimiter = "\s", data_start = 27) 
    
        # To read in column names
        data_all = open(filepath)
        for i, line in enumerate(data_all):
            if i > 13 :
                 header = header + [line]
            if i > 23:
                break
        data_all.close()
    
        for i in range(10):
            header[i] = "Column " + str(i+2) + header[i][ header[i].index(":") :] 
    
        header = ['DateTime'] + header
    
        #To format column names 
        for key1 in range(11,21):
            if key1 <14:
                header = header + ["Column "+ str(key1) + ': ' +'4He total counts for '+ (header[key1-10])[11:28] +' energy range']  
            else:
                header = header + ["Column "+ str(key1) + ': ' +'4He total counts for '+ (header[key1-10])[12:29] +' energy range']  
                
    
        data_modify = []
        
        #Combining Date, Time columns to make a single datetime.datetime value column 
        for i in range(len(data)): 
            date = datetime(data['col1'][i], 1, 1) + timedelta(int(data['col2'][i]) - 1)
            data_modify = data_modify + [datetime(date.year, date.month, date.day, data['col3'][i], data['col4'][i], data['col5'][i])]
    
        #Adding one DateTime column and removing 5 columns with separated time info
        data.add_column(Column(data = data_modify, name='col'),1)
        data.remove_columns(['col1', 'col2','col3','col4','col5'])
    
    
        #To add the column names in the astropy table object
        for key2 in range(21):
            data.rename_column(data.colnames[key2], header[key2])        
    
        #Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build of astropy and 
        # none of the stable versions include it
        data = data.to_pandas()
    
        return header, data
    
    """
    >>> _parse_txt('SIT_Ahead_10min_4He_2007_01.txt') 
    #Assuming the file is in same directory
    
    """
