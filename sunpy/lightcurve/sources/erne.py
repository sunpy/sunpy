from datetime import timedelta,datetime
from astropy.io import ascii
from astropy.table import Table, Column

from sunpy.time import TimeRange

def _parse_txt(filepath):
    """
    Parses a SOHO/ERNE file from
    http://srl.utu.fi/erne_data/carrot/carrota.html
    http://srl.utu.fi/erne_data/carrot/carrotp.html

    """
    
    #Reading in Data along with header
    data = ascii.read(filepath, delimiter = "\s", data_start = 2 ) 
    header = data.colnames
    for i in range(len(header)):
        header[i] = data[0][header[i]]
    data = data[1:]

    data_modify = []
    #Converting separate datetime element into a single datetime.datetime column
    for i in range(len(data)):
        if int(data['col1'][i]) >= 96:
            year = 1900 + int(data['col1'][i])
        else:
            year = 2000 + int(data['col1'][i])
       
        data_modify = data_modify + [TimeRange(datetime(year,int(data['col2'][i]),int(data['col3'][i]),int(data['col4'][i][:2]),int(data['col4'][i][2:])), datetime(year,int(data['col2'][i]),int(data['col3'][i]),int(data['col5'][i][:2]),int(data['col5'][i][2:])))]
    
    data.remove_columns(['col1','col2','col3','col4','col5'])
    data.add_column(Column(data = data_modify, name='col_1'),0)
    
    #To modify header
    for i in range(len(header[5:])):
        header[i+5] = 'Intensities [1/(cm^2*sr*s*MeV)] in energy channel ' + header[i+5] +' [MeV] ' 
    header = ['TimeRange'] + header[5:]
    
    # To add the column names in the astropy table object
    for key2 in range(len(data.colnames)):
        data.rename_column(data.colnames[key2], header[key2])        

    # Converting from astropy.table.Table to pandas.Dataframe
    # to_pandas() bound method is only available in the latest development build and none of the stable
    data = data.to_pandas()

    return header, data

"""
_parse_txt('cr1907p.txt')
_parse_txt('cr1906p.txt') 
_parse_txt('cr1906a.txt')
_parse_txt('cr1907a.txt')

"""
