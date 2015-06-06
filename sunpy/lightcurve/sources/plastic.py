from datetime import timedelta,datetime
from astropy.io import ascii
from astropy.table import Table, Column

from sunpy.time import TimeRange

def _parse_txt(filepath):
    """
    Parses a STEREO PLASTIC file from
    http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/

    """
    #Reading in Header  
    with open(filepath) as f:
        f.readline()
        header = f.readline().strip('\n').split('\t')

    if header[4] == 'sec':
        type_of_data = '1min'
    elif header[4] == 'date and time':
        type_of_data = '10min'
    elif header[4] == 'Np [1/cc]':
        type_of_data = '1hr'
    else:
        raise ValueError 

    
    #Reading in Data only, using default i.e. 0 value for data_start keyword since all lines before data are commented
    data = ascii.read(filepath, delimiter = "\s", data_start = 2) 
   
    data_modify = []
    data_modify_other = []

    #Converting separate datetime element into a single datetime.datetime column
    if type_of_data == '1min':
        for i in range(len(data)):
            data_modify = data_modify + [(datetime.strptime(str(data['col7'][i]) + str(data['col6'][i])[1:], "%Y-%m-%d/%H:%M:%S.%f") )]
            data_modify_other = data_modify_other + [datetime.strptime(str(data['col8'][i]), "%Y-%m-%d/%H:%M:%S") ]
        
        data.remove_columns(['col1','col2','col3','col4','col5','col6','col7','col8'])
        data.add_column(Column(data = data_modify, name='col_1'),0)
        data.add_column(Column(data = data_modify_other, name='col_2'),1)
        header = ['Datetime'] + header[7:]
        
    elif type_of_data == '10min':
        for i in range(len(data)):
            data['col5'][i] = datetime.strptime(str(data['col5'][i]), "%Y-%m-%d/%H:%M:%S")
        
        data.remove_columns(['col1','col2','col3','col4'])
        header = ['Datetime'] + header[5:]

    elif type_of_data == '1hr':
        for i in range(len(data)):
            data['col4'][i] = datetime.strptime(str(data['col4'][i]), "%Y-%m-%d/%H:%M:%S")
        
        data.remove_columns(['col1','col2','col3'])
        header = ['Datetime'] + header[4:]

    else:
        raise ValueError
    
    # To add the column names in the astropy table object
    for key2 in range(len(data.colnames)):
        data.rename_column(data.colnames[key2], header[key2])        

    # Converting from astropy.table.Table to pandas.Dataframe
    # to_pandas() bound method is only available in the latest development build and none of the stable
    data = data.to_pandas()

    print data
    return header, data

"""
_parse_txt('STA_L2_PLA_1DMax_1hr_20140101_001_V09.txt')
_parse_txt('STA_L2_PLA_1DMax_10min_20140101_001_V09.txt') 
_parse_txt('STA_L2_PLA_1DMax_1min_20140101_001_V09.txt') 

"""
