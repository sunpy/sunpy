# -*- coding: utf-8 -*-
"""SOHO Dataretriever sources subclass definitions"""

from __future__ import absolute_import


__authors__ = ["Ankit Kumar"]
__email__   = "ankitkmr.iitk@gmail.com"

# This module was developed with funding from 
# Google Summer of Code 2015


from datetime import timedelta,datetime,time

from astropy.io import ascii
from astropy.table import Table, Column

from sunpy.time import TimeRange
from sunpy.lightcurve import LightCurve



__all__ = ['ERNELightCurve']

class ERNELightCurve(LightCurve):
    """
    SOHO ERNE LightCurve. Provides Two-hour averaged intensities in 1.8-51 MeV data from as back as 1996-Feb-13.
    ERNE data can be download as a 2-hour resolution Carrington rotation sets. 
    
    Where: 1906 is the Carrington rotation number (a running number of full solar rotations starting 
    from November 9, 1853, using a rotation period of 27.2753 days)
    p is for protons, and
    a is for Helium (alpha particles).

    args:   
        Timerange or start and end date, 
        atomic_specie

    KEYWORD VALUES:-  
        atomic_specie: 'proton' or 'alpha'
    
    Examples
    --------
    >>> from sunpy import lightcurve as lc
    >>> from sunpy.time import TimeRange
    >>> erne = lc.ERNELightCurve.create(TimeRange('2012/06/01', '2012/06/05'),'proton')
    >>> erne.peek()

    References
    ----------
    They are available at the srl server, and have file names of type
    | http://srl.utu.fi/erne_data/carrot/1906/cr1906p.txt
    | http://srl.utu.fi/erne_data/carrot/1906/cr1906a.txt

    """

    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a SOHO/ERNE file from
        http://srl.utu.fi/erne_data/carrot/carrota.html
        http://srl.utu.fi/erne_data/carrot/carrotp.html

        and returns header as a list and ERNE data as pandas dataframe
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
            date = datetime.strptime(year_col[i] + '-' + month_col[i] + '-' + day_col[i],"%y-%m-%d")
            date1 = datetime.combine(date, time(int(begin_time_col[i][:2])))               #start time
            
            if end_time_col[i] == '0000':
                date = date + timedelta(days = 1)

            date2 = datetime.combine(date, time(int(end_time_col[i][:2])))                 #end time
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
