from datetime import timedelta,datetime

from astropy.io import ascii
from astropy.table import Table, Column

from sunpy.time import TimeRange

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

"""
_parse_txt('sept_ahead_ele_asun_2015_001_1min_l2_v03.dat.txt')
_parse_txt('sept_ahead_ele_asun_2015_001_1h_l2_v03.dat.txt') 

"""
