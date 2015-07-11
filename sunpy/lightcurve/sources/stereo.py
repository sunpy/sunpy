# -*- coding: utf-8 -*-
"""STEREO LightCurve sources subclass definitions"""
from __future__ import absolute_import

__authors__ = ["Ankit Kumar"]
__email__ = "ankitkmr.iitk@gmail.com"

# This module was developed with funding from 
# Google Summer of Code 2015


from datetime import timedelta,datetime

from astropy.io import ascii
from astropy.table import Table, Column

from sunpy.time import TimeRange
from sunpy.lightcurve import LightCurve

__all__ = ['LETLightCurve', 'SITLightCurve', 'PLASTICLightCurve', 'SEPTLightCurve', 'HETLightCurve']

class LETLightCurve(LightCurve):
    """
    STEREO LET LightCurve. Provides data from as back as 2006.

    Parameters
		----------
		timerange: sunpy.time.TimeRange
		    time range for which data is to be downloaded.
		    Default value -  TimeRange('2007-01-01','2008-06-01')	

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

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> from sunpy.time import TimeRange
    >>> let = lc.LETLightCurve.create(TimeRange('2012/06/01', '2012/06/05'), stereo_spacecraft = 'ahead', duration_of_average = 1*u.h, 
                                            average_type = 'standard', specie = 'He')
    >>> let.peek()

    References
    ----------
    | http://www.srl.caltech.edu/STEREO/Level1/LET_public.html

    """

    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a STEREO LET file from
        http://www.srl.caltech.edu/STEREO/Public/LET_public.html

        """ 
        type_of_data = ''
        header = []
        data_read = {'data_type_read': False, 'header_line_reached': False }

        #To Read in header and get index of line where data starts
        fp = open(filepath)
        for i,line in enumerate(fp):
            #To Find the type of data, 27 day, 1 day, 1 hour, 10 min or 1 min
            if data_read['data_type_read'] == False:
                if i == 2 and line[:3] == 'Bin':
                    type_of_data = '27day'
                    data_read['data_type_read'] = True
                elif i in [5,6] and line[:4] == 'Time': 
                    type_of_data = line[17:].strip('\n')
                    data_read['data_type_read'] = True
                else:
                    continue

            # To figure out and customise header in case of 27 day type of data
            if type_of_data == '27day':
                if line[:3] == 'Bin':
                    header = header + [line]
                    continue
                elif line[:7] == 'Columns':
                    data_start = i + 4
                    data_points = int(line[20])
                
                    for k in range(data_points):
                        header = header + ['Uncertainty for ' + header[k]]
                        header[k] = 'Flux for ' + header[k]
                    
                    header = ['Datetime'] + header
                    break
               
            # To figure out and customise header in case of any other type of data since their files are formatted 
            # similarly but different from 27 day type of data
            else:
                if data_read['header_line_reached'] == False and line[:6] != 'Column':
                    continue
                elif line[:6] == 'Column':
                    header = header + [line]
                    data_read['header_line_reached'] = True
                    continue
                elif line == 'BEGIN DATA\n':
                    data_start = i+1
                    break      
        fp.close()
        
        #Reading in Data only, using default i.e. 0 value for data_start keyword since all lines before data are commented
        data = ascii.read(filepath, delimiter = "\s", data_start = data_start) 
       

        #Storing data columns in recognizable variables (Common to all types of data)
        year_col = data['col1']
        day_of_year_col = data['col2']
       
        data_modify = []
        #Converting separate datetime element into a single datetime.datetime column
        if type_of_data == '27day':
            for i in range(len(data)):
                data_modify.append(datetime(year_col[i],1,1) + timedelta(days = day_of_year_col[i]) )
            
            data.remove_columns(['col{}'.format(i) for i in range(1,3)])
            data.add_column(Column(data = data_modify, name='Datetime'),0)
        
        #For all other types of data since they share the same format    
        else:
            #Storing data columns in recognizable variables
            hour_col = data['col3']
            minutes_col = data['col4']
            seconds_col = data['col5']

            for i in range(len(data)):
                data_modify.append(datetime(year_col[i],1,1,hour_col[i],minutes_col[i],seconds_col[i]) + timedelta(days = int(day_of_year_col[i]-1)))

            header = ['Datetime'] + header[5:]

            data.remove_columns(['col{}'.format(i) for i in range(1,6)])
            data.add_column(Column(data = data_modify, name='Datetime'),0)
        
        # Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build and none of the stable
        data = data.to_pandas()
        
        return header, data


class SITLightCurve(LightCurve):
    """
    SIT LightCurve. Provides SIT data back to 2007-07.
    Most recent data is usually available one or two days late.

    Parameters
    ----------
        timerange: sunpy.time.TimeRange
            time range for which data is to be downloaded.
            Example value -  TimeRange('2007-01-01','2015-03-01')   

        stereo_spacecraft: string   
            Default value - ahead
            Possible values - ahead, behind    # corresponding to spacecraft location

        specie:  string
            Default value - 4He
            Possible values - 4He, Fe, H, O

        duration_of_average: string
        Default value - 15min
            Possible values - 1min, 10min, 1hr, 1day        #corresponding to duration over which data is averaged


    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> from sunpy.time import TimeRange
    >>> sit = lc.SITLightCurve.create(TimeRange('2012/06/01', '2012/06/05'), stereo_spacecraft = "ahead",
    									duration_of_average =  "1min", specie = "4He")

    References
    ----------
    | http://www.srl.caltech.edu/STEREO/Public/SIT_public.html

    """

    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a STEREO SIT file from
        http://www.srl.caltech.edu/STEREO/Public/SIT_public.html
        
        and returns header and astropy.Table object containing data
    
        """
        header = []
        # To read in column names
        data_all = open(filepath)
        for i, line in enumerate(data_all):
            if i > 13 :
                 header = header + [line]
            if line == 'BEGIN DATA\n':
                break
        data_all.close()

        specie = header[-2][header[-2].index(':')+3:header[-2].index('t')-1]
        header = header[:-2]

        for i in range(len(header)):
            header[i] = header[i][ header[i].index(":") + 2:] 

        header = ['DateTime'] + header

        
        #To format column names 
        for i in range(len(header)-2):
            header = header + [specie + ' total counts for '+ (header[i+2])[1:20] +' energy range']

        data = ascii.read(filepath, delimiter = "\s", data_start = 27) 
    
        data_modify = []
        
        #Storing data columns in recognizable variables
        year_col = data['col1']
        day_of_year_col = data['col2']
        hour_col = data['col3']
        minutes_col = data['col4']
        seconds_col = data['col5']
        
        #Combining Date, Time columns to make a single datetime.datetime value column 
        for i in range(len(data)): 
            date = datetime(year_col[i], 1, 1) + timedelta(int(day_of_year_col[i]) - 1)
            data_modify.append(datetime(date.year, date.month, date.day, hour_col[i], minutes_col[i], seconds_col[i]))
    
        #Adding one DateTime column and removing 5 columns with separated time info
        data.add_column(Column(data = data_modify, name='col'),1)
        data.remove_columns(['col{}'.format(i) for i in range(1,6)])
    
    
        #To add the column names in the astropy table object
        for elem, head_key in enumerate(header):
            data.rename_column(data.colnames[elem], head_key)        
    
        #Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build of astropy and 
        # none of the stable versions include it
        data = data.to_pandas()
    
        return header, data


class PLASTICLightCurve(LightCurve):

    """
    STEREO PLASTIC LightCurve. 

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
    
    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> from sunpy.time import TimeRange
    >>> plastic = lc.PLASTICLightCurve.create(TimeRange('2012/06/01', '2012/06/05'), stereo_spacecraft = 'ahead', 
                                                duration_of_average = 10*u.min)
    >>> plastic.peek()

    References
    ----------
    | http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/

    """

    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a STEREO PLASTIC file from
        http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/

        """
        #Reading in Header  
        with open(filepath) as f:
            f.readline()
            header = f.readline().strip('\n').split('\t')

        types = {'sec':'1min', 'date and time': '10min', 'Np [1/cc]': '1hr'}

        #Determining the type of data
        type_of_data = types.get(header[4])
        if type_of_data is None:
            raise ValueError('{var} is not one of the possibilities in the header it should be one of {poss}'.format(var=header[4], poss= types.keys()))

        
        #Reading in Data only, using default i.e. 0 value for data_start keyword since all lines before data are commented
        data = ascii.read(filepath, delimiter = "\s", data_start = 2) 
       
        data_modify = []
        data_modify_other = []

        #Converting separate datetime element into a single datetime.datetime column
        #Adding the combined datetime column and removing the separate datetime elements columns
        if type_of_data == '1min':
            date_and_time_col = data['col7']
            millisec_col = data['col6']
            kev_q_time_col = data['col8']

            for i in range(len(data)):
                data_modify.append((datetime.strptime(str(date_and_time_col[i]) + str(millisec_col[i])[1:], "%Y-%m-%d/%H:%M:%S.%f") ))
                data_modify_other.append(datetime.strptime(str(kev_q_time_col[i]), "%Y-%m-%d/%H:%M:%S"))
            
            data.remove_columns(['col{}'.format(i) for i in range(1,9)])
            data.add_column(Column(data = data_modify, name='col_1'),0)
            data.add_column(Column(data = data_modify_other, name='col_2'),1)
            header = ['Datetime'] + header[7:]
            
        elif type_of_data == '10min':
            date_and_time_col = data['col5']

            date_and_time_col = [datetime.strptime(str(var), "%Y-%m-%d/%H:%M:%S") for var in date_and_time_col]

            data.remove_columns(['col{}'.format(i) for i in range(1,5)])
            header = ['Datetime'] + header[5:]

        elif type_of_data == '1hr':
            date_and_time_col = data['col4']

            date_and_time_col = [datetime.strptime(str(var), "%Y-%m-%d/%H:%M:%S") for var in date_and_time_col]

            data.remove_columns(['col{}'.format(i) for i in range(1,4)])
            header = ['Datetime'] + header[4:]

        else:
            raise ValueError('Unrecognized type of data')
        
        # To add the column names in the astropy table object
        for elem, head_key in enumerate(header):
            data.rename_column(data.colnames[elem], head_key)         

        # Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build and none of the stable
        data = data.to_pandas()

        return header, data


class SEPTLightCurve(LightCurve):
    
    """
    STEREO SEPT LightCurve. 

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
            Possible values -       
            #corresponding to duration over which data is averaged

        specie:  string
            Default value - element
            Possible values - element, ion

        sensor_pointing: string
            Default value - asun
            Possible values - asun, sun, north, south, omni

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> from sunpy.time import TimeRange
    >>> sept = lc.SEPTLightCurve.create(TimeRange('2012/06/01', '2012/06/05'), stereo_spacecraft = 'ahead', 
                                                duration_of_average = 10*u.min, specie = 'element', sensor_pointing = 'asun')
    >>> sept.peek()

    References
    ----------
    | http://www2.physik.uni-kiel.de/stereo/data/sept/level2/

    """

    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a STEREO SEPT file from
        http://www2.physik.uni-kiel.de/stereo/index.php?doc=data

        """
        #Reading in Data only, using default i.e. 0 value for data_start keyword since all lines before data are commented
        data = ascii.read(filepath, delimiter = "\s") 

        #Header
        header = [  'Column 1 : DateTime (datetime.datetime)',
                    'Column 2 : Bin 02 (  45.0-  55.0 keV) electron intensity (float)',
                    'Column 3 : Bin 03 (  55.0-  65.0 keV) electron intensity (float)',
                    'Column 4 : Bin 04 (  65.0-  75.0 keV) electron intensity (float)',
                    'Column 5: Bin 05 (  75.0-  85.0 keV) electron intensity (float)',
                    'Column 6: Bin 06 (  85.0- 105.0 keV) electron intensity (float)',
                    'Column 7: Bin 07 ( 105.0- 125.0 keV) electron intensity (float)',
                    'Column 8: Bin 08 ( 125.0- 145.0 keV) electron intensity (float)',
                    'Column 9: Bin 09 ( 145.0- 165.0 keV) electron intensity (float)',
                    'Column 10: Bin 10 ( 165.0- 195.0 keV) electron intensity (float)',
                    'Column 11: Bin 11 ( 195.0- 225.0 keV) electron intensity (float)',
                    'Column 12: Bin 12 ( 225.0- 255.0 keV) electron intensity (float)',
                    'Column 13: Bin 13 ( 255.0- 295.0 keV) electron intensity (float)',
                    'Column 14: Bin 14 ( 295.0- 335.0 keV) electron intensity (float)',
                    'Column 15: Bin 15 ( 335.0- 375.0 keV) electron intensity (float)',
                    'Column 16: Bin 16 ( 375.0- 425.0 keV) electron intensity (float)',
                    'Column 17: Bin 02 Uncertainty (float)',
                    'Column 18: Bin 03 Uncertainty (float)',
                    'Column 19: Bin 04 Uncertainty (float)',
                    'Column 20: Bin 05 Uncertainty (float)',
                    'Column 21: Bin 06 Uncertainty (float)',
                    'Column 22: Bin 07 Uncertainty (float)',
                    'Column 23: Bin 08 Uncertainty (float)',
                    'Column 24: Bin 09 Uncertainty (float)',
                    'Column 25: Bin 10 Uncertainty (float)',
                    'Column 26: Bin 11 Uncertainty (float)',
                    'Column 27: Bin 12 Uncertainty (float)',
                    'Column 28: Bin 13 Uncertainty (float)',
                    'Column 29: Bin 14 Uncertainty (float)',
                    'Column 30: Bin 15 Uncertainty (float)',
                    'Column 31: Bin 16 Uncertainty (float)',
                    'Column 32: Accumulation time in seconds (float)'  ]
       

        data_modify = []

        #Storing data columns in recognizable variables
        year_col    = data['col2']
        day_of_year = data['col3']
        hour_col    = data['col4']
        minutes_col = data['col5']
        seconds_col = data['col6']

        #Converting first five columns into a single datetime.datetime column
        for i in range(len(data)): 
            date = datetime(year_col[i], 1, 1) + timedelta(int(day_of_year[i]) - 1)
            data_modify.append(datetime(date.year, date.month, date.day, hour_col[i], minutes_col[i], seconds_col[i]))

        #Removing separate datetime element columns
        data.remove_columns(['col{}'.format(i) for i in range(1,7)])
        #Adding the combined datetime column created above in data_modify
        data.add_column(Column(data = data_modify, name='col1'),0)

        #To add the column names in the astropy table object
        for elem, head_key in enumerate(header):
            data.rename_column(data.colnames[elem], head_key)        

        # Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build and none of the stable
        data = data.to_pandas()

        return header, data


class HETLightCurve(LightCurve):
    """
    STEREO HET LightCurve. 

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
            Possible values - 1*u.min, 15*u.min, 1*u.h, 12*u.h, 1*u.d       #corresponding to duration over which data is averaged

    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> from sunpy.time import TimeRange
    >>> het = lc.HETLightCurve.create(TimeRange('2012/06/01', '2012/06/05'), stereo_spacecraft = 'ahead', 
                                                duration_of_average = 15*u.min)
    >>> het.peek()

    References
    ----------
    | http://www.srl.caltech.edu/STEREO/DATA/HET/

    """

    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a STEREO HET file from
        http://www.srl.caltech.edu/STEREO/Public/HET_public.html

        """
        #Reading in Data only
        data = ascii.read(filepath, delimiter = "\s", data_start = 20)

        #Determining the type of data from the number of columns in data
        number_of_columns = len(data.colnames)
        if number_of_columns == 33:
            type_of_data = '1min'
        else:
            type_of_data = 'other'
     
        #Header
        header = [  'Column 2: Electron flux, 0.7-1.4 MeV, particles/(cm2-sr-sec-MeV)', 
                    'Column 3: Uncertainty (sigma) for 0.7-1.4 MeV electron flux',
                    'Column 4: Electron flux, 1.4-2.8 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 5: Uncertainty (sigma) for 1.4-2.8 MeV electron flux',
                    'Column 6: Electron flux, 2.8-4.0 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 7: Uncertainty (sigma) for 2.8-4.0 MeV electron flux',
                    'Column 8: Proton flux, 13.6-15.1 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 9: Uncertainty (sigma) for 13.6-15.1 MeV proton flux',
                    'Column 10: Proton flux, 14.9-17.1 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 11: Uncertainty (sigma) for 14.9-17.1 MeV proton flux',
                    'Column 12: Proton flux, 17.0-19.3 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 13: Uncertainty (sigma) for 17.0-19.3 MeV proton flux',
                    'Column 14: Proton flux, 20.8-23.8 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 15: Uncertainty (sigma) for 20.8-23.8 MeV proton flux',
                    'Column 16: Proton flux, 23.8-26.4 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 17: Uncertainty (sigma) for 23.8-26.4 MeV proton flux',
                    'Column 18: Proton flux, 26.3-29.7 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 19: Uncertainty (sigma) for 26.3-29.7 MeV proton flux',
                    'Column 20: Proton flux, 29.5-33.4 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 21: Uncertainty (sigma) for 29.5-33.4 MeV proton flux' ,
                    'Column 22: Proton flux, 33.4-35.8 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 23: Uncertainty (sigma) for 33.4-35.8 MeV proton flux' ,
                    'Column 24: Proton flux, 35.5-40.5 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 25: Uncertainty (sigma) for 35.5-40.5 MeV proton flux' ,
                    'Column 26: Proton flux, 40.0-60.0 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 27: Uncertainty (sigma) for 40.0-60.0 MeV proton flux' ,
                    'Column 28: Proton flux, 60.0-100.0 MeV, particles/(cm2-sr-sec-MeV)',
                    'Column 29: Uncertainty (sigma) for 60.0-100.0 MeV proton flux' ]
       

        data_modify = []

        #Storing data columns in recognizable variables (Common to all types of data)
        start_year_col  = data['col2']
        start_month_col = data['col3']
        start_date_col  = data['col4']
        start_time_col  = data['col5']

        #Adding Time Column based on type of data
        if type_of_data == '1min':
            header = ['Verse Number', 'Column 1 : DateTime'] + header 

            for i in range(len(data)): 
                date = datetime.strptime(str(start_year_col[i])+ '-' +start_month_col[i]+ '-' +"%02d"%start_date_col[i] + '/' + ("%04d"%start_time_col[i])[:2] + ':' + ("%04d"%start_time_col[i])[2:], '%Y-%b-%d/%H:%M' )
                data_modify.append(date)

            data.remove_columns(['col{}'.format(i) for i in range(2,6)])
        else:
            header = ['Verse Number','Column 1 : TimeRange'] + header
            
            #Storing data columns in recognizable variables
            end_year_col    = data['col6']
            end_month_col   = data['col7']
            end_date_col    = data['col8']
            end_time_col    = data['col9']

            for i in range(len(data)): 
                date1 = datetime.strptime(str(start_year_col[i])+ '-' +start_month_col[i]+ '-' +"%02d"%start_date_col[i] + '/' + ("%04d"%start_time_col[i])[:2] + ':' + ("%04d"%start_time_col[i])[2:], '%Y-%b-%d/%H:%M' )
                date2 = datetime.strptime(str(end_year_col[i])+ '-' +end_month_col[i]+ '-' +"%02d"%end_date_col[i] + '/' + ("%04d"%end_time_col[i])[:2] + ':' + ("%04d"%end_time_col[i])[2:], '%Y-%b-%d/%H:%M' )
                data_modify.append(TimeRange(date1,date2))

            data.remove_columns(['col{}'.format(i) for i in range(2,10)])

        data.add_column(Column(data = data_modify, name='col'),1)

        #To add the column names in the astropy table object
        for elem, head_key in enumerate(header):
            data.rename_column(data.colnames[elem], head_key)        

        # Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build and none of the stable
        data = data.to_pandas()
        
        return data

