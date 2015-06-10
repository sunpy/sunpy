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

    #header is read along with data i.e the first row in data is the header
    header = data.colnames
    #So extracting the first row elements and forming header list that contains the energy bins
    header = [data[0][key] for key in data.colnames]
    #and now excluding the first line to have only the data in data variable
    data = data[1:]

    data_modify = []

    #Storing data columns in recognizable variables
    year_col = data['col1']
    month_col = data['col2']
    date_col = data['col3']
    begin_time_col = data['col4']
    end_time_col = data['col5']

    #Converting separate datetime element into a single datetime.datetime column
    for i in range(len(data)):
        if int(year_col[i]) >= 96:
            year = 1900 + int(year_col[i])
        else:
            year = 2000 + int(year_col[i])
        
        #Combining separate datetime elements into single datetime value
        #start time
        date1 = datetime(year,int(month_col[i]),int(date_col[i]),int(begin_time_col[i][:2]),int(begin_time_col[i][2:]))
        #end time
        date2 = datetime(year,int(month_col[i]),int(date_col[i]),int(end_time_col[i][:2]),int(end_time_col[i][2:]))
        #Appending the start and end time as sunpy.time TimeRange in a separate list
        data_modify.append(TimeRange(date1, date2))
    
    #Removing the columns with separate date elements like year, day, hour, minute, second
    data.remove_columns(['col{}'.format(i) for i in range(1,6)])
    #Adding a single Timerange column to the data
    data.add_column(Column(data = data_modify, name='col_1'),0)
    
    #To modify header
    header[5:] = ['Intensities [1/(cm^2*sr*s*MeV)] in energy channel {}  [MeV] '.format(val) for val in header[5:]]
    header = ['TimeRange'] + header[5:]
    
    # To add the column names in the astropy table object
    for elem, head_key in enumerate(header):
        data.rename_column(data.colnames[elem], head_key)        

    # Converting from astropy.table.Table to pandas.Dataframe
    # to_pandas() bound method is only available in the latest development build of astropy and none of the stable versions :/
    data = data.to_pandas()
    
    return header, data

"""
_parse_txt('cr1907p.txt')
_parse_txt('cr1906p.txt') 
_parse_txt('cr1906a.txt')
_parse_txt('cr1907a.txt')

"""
