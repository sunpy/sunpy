from datetime import timedelta,datetime
from astropy.io import ascii
from astropy.table import Table, Column

from sunpy.time import TimeRange

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
                'Column 9 Uncertainty (sigma) for 13.6-15.1 MeV proton flux',
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
    month_dict = {'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,'Oct':10,'Nov':11,'Dec':12}


    #Adding Time Column based on type of data
    if type_of_data == '1min':
        header = ['Verse Number', 'DateTime'] + header

        for i in range(len(data)): 
            date = datetime(data['col2'][i], month_dict[ data['col3'][i] ], data['col4'][i], int(("%04d" % (data['col5'][i],))[:2]), int(("%04d" % (data['col5'][i],))[2:]) )
            data_modify = data_modify + [date]

        data.remove_columns(['col2','col3','col4','col5'])
    else:
        header = ['Verse Number','TimeRange'] + header

        for i in range(len(data)): 
            date1 = datetime(data['col2'][i], month_dict[ data['col3'][i] ], data['col4'][i], int(("%04d" % (data['col5'][i],))[:2]), int(("%04d" % (data['col5'][i],))[2:]) )
            date2 = datetime(data['col6'][i], month_dict[ data['col7'][i] ], data['col8'][i], int(("%04d" % (data['col9'][i],))[:2]), int(("%04d" % (data['col9'][i],))[2:]) )
            data_modify = data_modify + [TimeRange(date1,date2)]

        data.remove_columns(['col2','col3','col4','col5','col6','col7','col8','col9'])

    data.add_column(Column(data = data_modify, name='col'),1)

    #To add the column names in the astropy table object
    for key2 in range(len(data.colnames)):
        data.rename_column(data.colnames[key2], header[key2])        

    # Converting from astropy.table.Table to pandas.Dataframe
    # to_pandas() bound method is only available in the latest development build and none of the stable
    data = data.to_pandas()
    
    return header, data

"""
_parse_txt('AeH07Apr.1m.txt')
_parse_txt('AeH06Dec.15m.txt') 

"""
