from datetime import timedelta,datetime

from astropy.io import ascii
from astropy.table import Table, Column


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

"""
_parse_txt('STA_L2_PLA_1DMax_1hr_20140101_001_V09.txt')
_parse_txt('STA_L2_PLA_1DMax_10min_20140101_001_V09.txt') 
_parse_txt('STA_L2_PLA_1DMax_1min_20140101_001_V09.txt') 

"""
