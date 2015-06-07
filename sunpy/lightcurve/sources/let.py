from datetime import timedelta,datetime
from astropy.io import ascii
from astropy.table import Table, Column

def _parse_txt(filepath):
    """
    Parses a STEREO LET file from
    http://www.srl.caltech.edu/STEREO/Public/LET_public.html

    """
    #Reading in Header  

    type_of_data = ''
    header = []
    flag1 = -1
    flag2 = -1

    #To Read in header and get index of line where data starts
    fp = open(filepath)
    for i,line in enumerate(fp):
        #To Find the typr of data, 27 day, 1 day, 1 hour, 10 min or 1 min
        if flag1 == -1:
            if i == 2 and line[:3] == 'Bin':
                type_of_data = '27day'
                flag1 = 1
            elif i in [5,6] and line[:4] == 'Time': 
                type_of_data = line[17:].strip('\n')
                flag1 = 2
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
            if flag2 == -1 and line[:6] != 'Column':
                continue
            elif line[:6] == 'Column':
                header = header + [line]
                flag2 = 0
                continue
            elif line == 'BEGIN DATA\n':
                data_start = i+1
                break      
    fp.close()
    
    #Reading in Data only, using default i.e. 0 value for data_start keyword since all lines before data are commented
    data = ascii.read(filepath, delimiter = "\s", data_start = data_start) 
   
    data_modify = []
   

    #Converting separate datetime element into a single datetime.datetime column
    if type_of_data == '27day':
        for i in range(len(data)):
            data_modify = data_modify + [datetime(data['col1'][i],1,1) + timedelta(days = data['col2'][i]) ]
        
        data.remove_columns(['col1','col2'])
        data.add_column(Column(data = data_modify, name='Datetime'),0)
        
    else: 
        for i in range(len(data)):
            data_modify = data_modify + [datetime(data['col1'][i],1,1,data['col3'][i],data['col4'][i],data['col5'][i]) + timedelta(days = int(data['col2'][i]-1)) ]
        header = ['Datetime'] + header[5:]
        data.remove_columns(['col1','col2','col3','col4','col5'])
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
