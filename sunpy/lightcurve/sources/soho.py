# -*- coding: utf-8 -*-
"""SOHO Dataretriever sources subclass definitions"""

from __future__ import absolute_import


__authors__ = ["Ankit Kumar"]
__email__   = "ankitkmr.iitk@gmail.com"

# This module was developed with funding from 
# Google Summer of Code 2015


from datetime import timedelta,datetime,time

import pandas as pd
from pandas import DataFrame
import numpy as np

from astropy.io import ascii
from astropy.table import Table, Column
from astropy.utils import OrderedDict
from astropy.table.column import MaskedColumn

from sunpy.time import TimeRange
from sunpy.lightcurve import LightCurve



__all__ = ['ERNELightCurve']

def _to_pandas(self):
    """
    Return a :class:`pandas.DataFrame` instance

    Returns
    -------
    dataframe : :class:`pandas.DataFrame`
        A pandas :class:`pandas.DataFrame` instance

    Raises
    ------
    ImportError
        If pandas is not installed
    ValueError
        If the Table contains mixin or multi-dimensional columns
    """
    from pandas import DataFrame

    if self.has_mixin_columns:
        raise ValueError("Cannot convert a table with mixin columns to a pandas DataFrame")

    if any(getattr(col, 'ndim', 1) > 1 for col in self.columns.values()):
        raise ValueError("Cannot convert a table with multi-dimensional columns to a pandas DataFrame")

    out = OrderedDict()

    for name, column in self.columns.items():
        if isinstance(column, MaskedColumn):
            if column.dtype.kind in ['i', 'u']:
                out[name] = column.astype(float).filled(np.nan)
            elif column.dtype.kind in ['f', 'c']:
                out[name] = column.filled(np.nan)
            else:
                out[name] = column.astype(np.object).filled(np.nan)
        else:
            out[name] = column

        if out[name].dtype.byteorder not in ('=', '|'):
            out[name] = out[name].byteswap().newbyteorder()

    return DataFrame(out)
    

class ERNELightCurve(LightCurve):
    """
    SOHO ERNE LightCurve. Provides Two-hour averaged intensities in 1.8-51 MeV data from as back as 1996-Feb-13.
    ERNE data can be download as a 2-hour resolution Carrington rotation sets. 
    
    Where: 1906 is the Carrington rotation number (a running number of full solar rotations starting 
    from November 9, 1853, using a rotation period of 27.2753 days)
    p is for protons, and
    a is for Helium (alpha particles).

    Parameters
    ----------
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

    def peek(self, title="ERNE Two-hour averaged intensities"):
        """Plots ERNE light curve in the usual manner"""

        figure = plt.figure()
        ax = plt.gca()

        timerange_start = self.data['TimeRange'].apply(lambda col: col.start)
        dates = matplotlib.dates.date2num(timerange_start.astype(datetime))

        colors = ['Green','Red','Chocolate', 'Blue','SeaGreen','Tomato',
                            'SlateBlue','Orange','Purple','Magenta','MediumVioletRed']
        figure.delaxes(ax)
        axes = figure.add_axes([0.1, 0.15, 0.55, 0.8])

        for i,line in enumerate(self.header):
            if i >= 1:
                axes.plot_date(dates, self.data[line].ffill(), '-',
                     label=line[line.index('l')+2:], color=colors[i], lw=1)
        
        axes.set_yscale("log",nonposy='mask')
        axes.set_title(title)
        axes.set_ylabel('1/(cm^2*sr*s*MeV) [per nucleon in case of protons]')
        axes.set_xlabel('UTC TimeZone')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        figure.autofmt_xdate()
        plt.show()

        return figure

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
            date = datetime.strptime(year_col[i] + '-' + month_col[i] + '-' + date_col[i],"%y-%m-%d")
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
        header = ['energy channel {} MeV'.format(val) for val in header[5:]]
        header = ['TimeRange'] + header
        
        # To add the column names in the astropy table object
        for elem, head_key in enumerate(header):
            data.rename_column(data.colnames[elem], head_key)        

        # Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build of astropy and none of the stable versions :/
        data = _to_pandas(data)
        for i,line in enumerate(header[1:]): 
            data[line] = data[line].apply(lambda col: float(col))

        return header, data

