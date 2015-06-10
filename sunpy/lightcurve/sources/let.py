from datetime import timedelta,datetime

from astropy.io import ascii
from astropy.table import Table, Column

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

"""

_parse_txt('Ar_summed_ahead.txt')
_parse_txt('Al_summed_ahead.txt')
_parse_txt('CNO_lo_sectored_ahead_2015_1hr_level1_11.txt') 
_parse_txt('Fe_sectored_ahead_2015_001_level1_11.txt')
_parse_txt('Fe_sectored_ahead_2015_01_10min_level1_11.txt') 

"""
