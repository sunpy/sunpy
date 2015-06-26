# -*- coding: utf-8 -*-
"""SOHO LightCurve sources subclass definitions"""

__authors__ = ["Ankit Kumar"]
__email__   = "ankitkmr.iitk@gmail.com"

# This module was developed with funding from 
# Google Summer of Code 2015

from __future__ import absolute_import
from datetime import timedelta,datetime

from astropy.io import ascii
from astropy.table import Table, Column

from sunpy.time import TimeRange


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
        
        # def peek(self, title="GOES Xray Flux"):
        #     """Plots GOES light curve is the usual manner"""
        #     figure = plt.figure()
        #     axes = plt.gca()
        
        #     dates = matplotlib.dates.date2num(parse_time(self.data.index))
        
        #     axes.plot_date(dates, self.data['xrsa'], '-',
        #                  label='0.5--4.0 $\AA$', color='blue', lw=2)
        #     axes.plot_date(dates, self.data['xrsb'], '-',
        #                  label='1.0--8.0 $\AA$', color='red', lw=2)
        
        #     axes.set_yscale("log")
        #     axes.set_ylim(1e-9, 1e-2)
        #     axes.set_title(title)
        #     axes.set_ylabel('Watts m$^{-2}$')
        #     axes.set_xlabel(datetime.datetime.isoformat(self.data.index[0])[0:10])
        
        #     ax2 = axes.twinx()
        #     ax2.set_yscale("log")
        #     ax2.set_ylim(1e-9, 1e-2)
        #     ax2.set_yticks((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
        #     ax2.set_yticklabels((' ', 'A', 'B', 'C', 'M', 'X', ' '))
        
        #     axes.yaxis.grid(True, 'major')
        #     axes.xaxis.grid(False, 'major')
        #     axes.legend()
        
        #     # @todo: display better tick labels for date range (e.g. 06/01 - 06/05)
        #     formatter = matplotlib.dates.DateFormatter('%H:%M')
        #     axes.xaxis.set_major_formatter(formatter)
        
        #     axes.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
        #     figure.autofmt_xdate()
        #     figure.show()
        
        #     return figure


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
                    
                        #Combining separate datetime elements into single datetime value
                        #start time
                        date1 = datetime.strptime(year_col[i] + '-' + month_col[i] + '-' + date_col[i] + '/' + begin_time_col[i][:2] + ':' +begin_time_col[i][2:],"%y-%m-%d/%H:%M")
                        #end time
                        date2 = datetime.strptime(year_col[i] + '-' + month_col[i] + '-' + date_col[i] + '/' + end_time_col[i][:2] + ':' + end_time_col[i][2:], "%y-%m-%d/%H:%M")
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
