from datetime import timedelta,datetime

from astropy.io import ascii
from astropy.table import Column


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
            
            #I don't know why but including it outside interferes with the upper case and results in segmetation fault 11 :/
            from sunpy.time import TimeRange
            data_modify.append(TimeRange(date1,date2))

        data.remove_columns(['col{}'.format(i) for i in range(2,10)])

    data.add_column(Column(data = data_modify, name='col'),1)

    #To add the column names in the astropy table object
    for elem, head_key in enumerate(header):
        data.rename_column(data.colnames[elem], head_key)        

    # Converting from astropy.table.Table to pandas.Dataframe
    # to_pandas() bound method is only available in the latest development build and none of the stable
    data = data.to_pandas()
    print data
    return data

""""
_parse_txt('AeH07Apr.1m.txt')
_parse_txt('AeH06Dec.15m.txt') 

"""
