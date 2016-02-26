# -*- coding: utf-8 -*-
"""STEREO LightCurve sources subclass definitions"""

from __future__ import absolute_import

__authors__ = ["Ankit Kumar"]
__email__ = "ankitkmr.iitk@gmail.com"

# This module was developed with funding from 
# Google Summer of Code 2015


from datetime import timedelta,datetime

import pandas as pd
from pandas import DataFrame
import matplotlib.dates
from matplotlib import pyplot as plt
import numpy as np

from astropy.io import ascii
from astropy.table import Table, Column
from astropy.utils import OrderedDict
from astropy.table.column import MaskedColumn

from sunpy.time import TimeRange, parse_time
from sunpy.lightcurve import LightCurve
from sunpy.lightcurve.table_to_dataframe import _to_pandas


__all__ = ['LETLightCurve', 'SITLightCurve', 'PLASTICLightCurve', 'SEPTLightCurve', 'HETLightCurve']


class LETLightCurve(LightCurve):
    """
    STEREO LET LightCurve. Provides data from as back as 2006. Currently the LightCurve supports only Single File Load

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
       time range for which data is to be downloaded.
       Example value -  TimeRange('2007-01-01','2008-06-01')   

    duration_of_average: string
       Default value - 10 * u.min
       Possible values - 1 * u.min, 10 * u.min, 1 * u.h, 1 * u.d, 27 * u.d   
       corresponding to duration over which data is averaged

    type_of_data:  string
       Possible values - depends on other parameters
       if duration_of_average = 27 * u.d:
       Possible Values: summed, narrow
       else:
       Possible values: sectored, standard, summed

    species:  string
       Possible values - depends on other parameters
       if type_of_data = 'Sectored' and duration_of_average in [ 1 * u.min, 10 * u.min, 1 * u.h, 1 * u.d]:
       Possible values: CNO_hi,CNO_lo, Fe_hi, Fe_lo, H_lo, He3_lo, He4_hi, He4_lo, He_lo, NeMgSi_hi, NeMgSi_lo
       else:
       Possible values: Al, Ar, C, Ca, Fe, H, He, He3, He4, Mg, N, Na, Ne, Ni, O, S, Si

    stereo_spacecraft: string   
       Default value - ahead
       Possible values - ahead, behind    
       corresponding to spacecraft location

    Examples
    --------
    
    .. plot::


        import os
        import sunpy.data.test
        filepath = sunpy.data.test.rootdir
        from sunpy import lightcurve as lc
        [header,data] = lc.LETLightCurve._parse_txt(os.path.join(filepath , 'let', 'Ar_ahead_2006_318_level1_11.txt'))
        let = lc.LETLightCurve(data,header)
        let.peek()

    References
    ----------
    | http://www.srl.caltech.edu/STEREO/Level1/LET_public.html

    """

    def peek(self, title="LET electron/proton Flux"):
        """Plots LET light curve in the usual manner"""
        figure = plt.figure()
        ax = plt.gca()

        dates = matplotlib.dates.date2num(self.data['Datetime'].astype(datetime))

        colors = ['Green','Red', 'Blue','SeaGreen','Tomato','SlateBlue','Orange',
        'Purple','Magenta','Chocolate','MediumVioletRed', 'Teal','Navy','Indigo']

        figure.delaxes(ax)
        axes = figure.add_axes([0.1, 0.15, 0.55, 0.8])

        if self.header[1][:4] == 'Flux':
            #27 day data
            num_energy_bins = (len(self.header)-1)/2
            for i,line in enumerate(self.header):
                if i >= 1 and i < num_energy_bins:
                    axes.plot_date(dates, self.data['col' +str(3+i)].ffill(), '-', label= line[9:], color=colors[i], lw=0.5)
            axes.set_yscale("log",nonposy='mask')
            
        elif self.header[1][:22] == 'Column 6: LET Livetime' and self.header[2][9] == '-':
            #Sectored data
            num_energy_bins = self.header[2][self.header[2].index(':')+2:self.header[2].index('s')-1]
            for i in range(int(num_energy_bins)):
                axes.plot_date(dates, self.data['col' + str(7+i)].ffill(), '-', label= 'Sector '+str(i), color=colors[i/2-2], lw=2)

        elif self.header[2][:22] == 'Column 7: LET Livetime':
            #Non sectored standard data
            num_energy_bins = len(self.header) - 4
            for i in range(num_energy_bins):
                label = self.header[3+i][self.header[3+i].index(':')+1:self.header[3+i].index('/')+2]
                axes.plot_date(dates, self.data['col' + str(8+i)].ffill(), '-', label= label, color=colors[i], lw=2)
        else:
            #Non sectored non standard data
            num_energy_bins = len(self.header) - 4
            for i in range(num_energy_bins):
                label = self.header[2+i][self.header[2+i].index(':')+1:self.header[2+i].index('/')+4]
                axes.plot_date(dates, self.data['col' + str(7+i)].ffill(), '-', label= label, color=colors[i], lw=2)


        
        axes.set_title(title)
        axes.set_ylabel('1/(cm^2 s sr MeV/nuc)')
        axes.set_xlabel('UTC Time')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')

        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        figure.autofmt_xdate()
        plt.show()

        return figure

    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a STEREO LET file from
        http://www.srl.caltech.edu/STEREO/Public/LET_public.html

        """ 
        type_of_data = ''
        header = []
        data_read = {'data_type_read': False, 'header_line_reached': False }

        #To Read in header and get index of line where data starts
        fp = open(filepath)
        for i,line in enumerate(fp):
            #To Find the type of data, 27 day, 1 day, 1 hour, 10 min or 1 min
            if not data_read['data_type_read']:
                if i == 2 and line[:3] == 'Bin':
                    type_of_data = '27day'
                    data_read['data_type_read'] = True
                elif i in [5,6] and line[:4] == 'Time': 
                    type_of_data = line[17:].strip('\n')
                    data_read['data_type_read'] = True
                else:
                    continue

            # To figure out and customise header in case of 27 day type of data
            if type_of_data == '27day':
                if line[:3] == 'Bin':
                    header = header + [line]
                    continue
                elif line[:7] == 'Columns':
                    data_start = i + 4
                    data_points = int(line[20:line.index('f')-1])
                
                    for k in range(data_points):
                        header = header + ['Uncertainty for ' + header[k]]
                        header[k] = 'Flux for ' + header[k]
                    
                    header = ['Datetime'] + header
                    break
               
            # To figure out and customise header in case of any other type of data since their files are formatted 
            # similarly but different from 27 day type of data
            else:
                if line[:6] == 'Column':
                    header = header + [line]
                    data_read['header_line_reached'] = True
                    continue
                elif line == 'BEGIN DATA\n':
                    data_start = i+1
                    break      
        fp.close()
        
        #Reading in Data only, using default i.e. 0 value for data_start keyword since all lines before data are commented
        data = ascii.read(filepath, delimiter = "\s", data_start = data_start) 
       

        #Storing data columns in recognizable variables (Common to all types of data)
        year_col = data['col1']
        day_of_year_col = data['col2']
       
        data_modify = []
        #Converting separate datetime element into a single datetime.datetime column
        if type_of_data == '27day':
            data_modify = [datetime.strptime('{0} {1:03d}'.format(year, int(day)), '%Y %j') for year, day in zip(year_col, day_of_year_col)]
            
            data.remove_columns(['col{}'.format(i) for i in range(1,3)])
            data.add_column(Column(data = data_modify, name='Datetime'),0)
        
        #For all other types of data since they share the same format    
        else:
            #Storing data columns in recognizable variables
            hour_col = data['col3']
            minutes_col = data['col4']
            seconds_col = data['col5']

            data_modify = [datetime.strptime('{0} {1:03d} {2} {3} {4}'.format(year, int(day), hour, minute, second), '%Y %j %H %M %S') \
                        for year, day, hour, minute, second in zip(year_col, day_of_year_col, hour_col, minutes_col, seconds_col)]
            
            header = ['Datetime'] + header[5:]

            data.remove_columns(['col{}'.format(i) for i in range(1,6)])
            data.add_column(Column(data = data_modify, name='Datetime'),0)
        
        # Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build and none of the stable
        data = _to_pandas(data)
        data = data.replace([-1.0e+31, -9999.9],float('nan'))

        return OrderedDict(enumerate(header)), data


class SITLightCurve(LightCurve):
    """
    SIT LightCurve. Provides SIT data back to 2007-07.
    Most recent data is usually available one or two days late. 
    Currently the LightCurve supports only Single File Load

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value -  TimeRange('2007-01-01','2015-03-01')   

    stereo_spacecraft: string   
        Default value - ahead
        Possible values - ahead, behind    # corresponding to spacecraft location

    species:  string
        Default value - 4He
        Possible values - 4He, Fe, H, O

    duration_of_average: string
    Default value - 15min
        Possible values - 1min, 10min, 1hr, 1day        #corresponding to duration over which data is averaged


    Examples
    --------

    .. plot::

    
        import os
        import sunpy.data.test
        filepath = sunpy.data.test.rootdir
        from sunpy import lightcurve as lc
        [header,data] = lc.SITLightCurve._parse_txt(os.path.join(filepath , 'sit', 'SIT_Ahead_10min_H_2007_01.txt'))
        sit = lc.SITLightCurve(data,header)
        sit.peek()


    References
    ----------
    | http://www.srl.caltech.edu/STEREO/Public/SIT_public.html

    """

    def peek(self, title="SIT elemental Intensities"):
        """Plots SIT light curve in the usual manner"""

        figure = plt.figure()
        ax = plt.gca()

        dates = matplotlib.dates.date2num(self.data['DateTime'].astype(datetime))

        num_energy_bins = (len(self.header.values())-2)/2
        colors = ['Green','Red','Chocolate', 'Blue','SeaGreen','Tomato','SlateBlue','Orange','Purple','Magenta','MediumVioletRed']

        figure.delaxes(ax)
        axes = figure.add_axes([0.1, 0.15, 0.55, 0.8])

        for i,line in enumerate(self.header.values()):
            if i >= 2 and i <= num_energy_bins + 1:
                axes.plot_date(dates, self.data[line].ffill(), '-',
                     label=line[1:20], color=colors[i-2], lw=0.5)
        
        axes.set_yscale("log",nonposy='mask')
        axes.set_ylim(1e-3, 1e+3)
        axes.set_title('{title} : {species}'.format(title = title, species = self.header.values()[-1][:self.header.values()[-1].index(' ')]))
        axes.set_ylabel('1/(cm^2 s sr MeV/nuc)')
        axes.set_xlabel('UTC Time')


        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        figure.autofmt_xdate()
        plt.show()

        return figure


    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a STEREO SIT file from
        http://www.srl.caltech.edu/STEREO/Public/SIT_public.html
        
        and returns header and Pandas.Dataframe containing data
    
        """
        header = []
        # To read in column names
        data_all = open(filepath)
        for i, line in enumerate(data_all):
            if i > 13 :
                 header.append(line)
            if line == 'BEGIN DATA\n':
                data_start = i+1
                break
        data_all.close()

        # To extract the species (4He) from a line like
        # Columns 17 - 26:  4He total counts for same energy ranges as above
        specie = header[-2].split(": ")[1].split()[0]

        header = header[:-2]
        header = [ h.split(": ")[1] for h in header]
        header = ['DateTime'] + header

        
        #To format column names 
        header = header + ['{} total counts for {} energy range'.format(specie, h[1:20]) for h in header[2:]]

        data = ascii.read(filepath, delimiter = "\s", data_start = data_start) 
    
        data_modify = []
        
        #Storing data columns in recognizable variables
        year_col = data['col1']
        day_of_year_col = data['col2']
        hour_col = data['col3']
        minutes_col = data['col4']
        seconds_col = data['col5']
        
        #Combining Date, Time columns to make a single datetime.datetime value column 
        data_modify = [datetime.strptime('{0} {1:03d} {2} {3} {4}'.format(year, int(day), hour, minute, second), '%Y %j %H %M %S') \
                    for year, day, hour, minute, second in zip(year_col, day_of_year_col, hour_col, minutes_col, seconds_col)]

        #Adding one DateTime column and removing 5 columns with separated time info
        data.add_column(Column(data = data_modify, name='col'),1)
        data.remove_columns(['col{}'.format(i) for i in range(1,6)])
    
    
        #To add the column names in the astropy table object
        for elem, head_key in enumerate(header):
            data.rename_column(data.colnames[elem], head_key)        
    
        #Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build of astropy and 
        # none of the stable versions include it
        data = _to_pandas(data)
    
        return OrderedDict(enumerate(header)), data


class PLASTICLightCurve(LightCurve):
    """
    STEREO PLASTIC LightCurve. Currently the LightCurve supports only Single File Load

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value -  TimeRange('2007-02-14','2014-12-17')   

    stereo_spacecraft: string   
        Default value - ahead
        Possible values - ahead, behind    # corresponding to spacecraft location


    duration_of_average: astropy units ( Default value - 10 * u.min )
        Possible values - 1 * u.min, 10 * u.min, 1 * u.h      
        #corresponding to duration over which data is averaged

    Examples
    --------

    .. plot::

    
        import os
        import sunpy.data.test
        filepath = sunpy.data.test.rootdir
        from sunpy import lightcurve as lc
        [header,data] = lc.PLASTICLightCurve._parse_txt(os.path.join(filepath , 'plastic', 'STA_L2_PLA_1DMax_1min_20140101_001_V09.txt'))
        plastic = lc.PLASTICLightCurve(data,header)
        plastic.peek()

    References
    ----------
    | http://stereo-ssc.nascom.nasa.gov/data/ins_data/plastic/level2/Protons/ASCII/

    """

    def peek(self, title=" PLASTIC Data"):
        """Plots PLASTIC light curve in the usual manner"""

        figure = plt.figure()
        
        dates = matplotlib.dates.date2num(self.data['Datetime'].apply(lambda col: parse_time(str(col).replace('/',' '))).astype(datetime))
        
        colors = ['Green','Red','Chocolate', 'Blue','SeaGreen','Tomato','SlateBlue','Orange',
        'Purple','Magenta','MediumVioletRed', 'Teal','Navy','Indigo']

        ax1 = figure.add_subplot(3,1,1)
        plt.plot_date(dates, self.data['Np [1/cc]'].ffill(), '-',label= 'Np [1/cc]' , color='Red', lw=2)
        ax2 = figure.add_subplot(3,1,2)
        plt.plot_date(dates, self.data['Bulk Speed [km/s]'].ffill(), '-',label= 'Bulk Speed [km/s]' , color= 'Blue', lw=2)
        ax3 = figure.add_subplot(3,1,3)
        plt.plot_date(dates, self.data['Tkin [deg K]'].ffill(), '-',label= 'Tkin [deg K]' , color='Teal', lw=2)
        
        ax1.set_yscale("log",nonposy = "mask")
        ax1.set_ylabel('[1/cc]')
        ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2),
              ncol=1, fancybox=True, shadow=True)
        ax1.yaxis.grid(True, 'major')
        ax1.xaxis.grid(False, 'major')

        ax2.set_yscale("log",nonposy = "mask")
        ax2.set_ylabel('[km/s]')
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2),
              ncol=1, fancybox=True, shadow=True)
        ax2.yaxis.grid(True, 'major')
        ax2.xaxis.grid(False, 'major')

        ax3.set_yscale("log", nonposy = "mask")
        ax3.set_ylabel('[deg K]')
        ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2),
              ncol=1, fancybox=True, shadow=True)
        ax3.yaxis.grid(True, 'major')
        ax3.xaxis.grid(False, 'major')
        ax3.set_xlabel('UTC Time')

        figure.autofmt_xdate()
        plt.show()

        return figure


    @staticmethod
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
            raise ValueError('{var} is not one of the possibilities in the header it \
                should be one of {poss}'.format(var=header[4], poss= types.keys()))

        
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
                data_modify.append(datetime.strptime(date_and_time_col[i] + str(millisec_col[i])[1:], "%Y-%m-%d/%H:%M:%S.%f"))
                data_modify_other.append(datetime.strptime(kev_q_time_col[i], "%Y-%m-%d/%H:%M:%S"))
            
            data.remove_columns(['col{}'.format(i) for i in range(1,9)])
            data.add_column(Column(data = data_modify, name='col_1'),0)
            data.add_column(Column(data = data_modify_other, name='col_2'),1)
            header = ['Datetime'] + header[7:]
            
        elif type_of_data == '10min':
            date_and_time_col = [datetime.strptime(var, "%Y-%m-%d/%H:%M:%S") for var in data['col5']]

            data.remove_columns(['col{}'.format(i) for i in range(1,5)])
            header = ['Datetime'] + header[5:]

        elif type_of_data == '1hr':
            date_and_time_col = [datetime.strptime(var, "%Y-%m-%d/%H:%M:%S") for var in data['col4']]

            data.remove_columns(['col{}'.format(i) for i in range(1,4)])
            header = ['Datetime'] + header[4:]

        else:
            raise ValueError('Unrecognized type of data')
        
        # To add the column names in the astropy table object
        for elem, head_key in enumerate(header):
            data.rename_column(data.colnames[elem], head_key)         

        # Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build and none of the stable
        data = _to_pandas(data)
        data = data.replace([-1.0E+31, -9999.9],float('nan'))


        return OrderedDict(enumerate(header)), data


class SEPTLightCurve(LightCurve):
    """
    STEREO SEPT LightCurve. Currently the LightCurve supports only Single File Load

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value -  TimeRange('2007-01-20','2015-01-01')   

    stereo_spacecraft: string   
        Default value - ahead
        Possible values - ahead, behind    # corresponding to spacecraft location

    duration_of_average: astropy units quantity ( Default value - 10 * u.min )
        Possible values -  1 * u.min, 10 * u.min, 1 * u.h, 1 * u.d    
        #corresponding to duration over which data is averaged

    species: string
        Default value - element
        Possible values - element, ion

    sensor_pointing: string
        Default value - asun
        Possible values - asun, sun, north, south, omni

    Examples
    --------

    .. plot::

    
        import os
        import sunpy.data.test
        filepath = sunpy.data.test.rootdir
        from sunpy import lightcurve as lc
        [header,data] = lc.SEPTLightCurve._parse_txt(os.path.join(filepath , 'sept', 'sept_ahead_ele_asun_2015_001_1min_l2_v03.dat.txt'))
        sept = lc.SEPTLightCurve(data,header)
        sept.peek()

    References
    ----------
    | http://www2.physik.uni-kiel.de/stereo/data/sept/level2/

    """

    def peek(self, title="SEPT Electron Intensities"):
        """Plots SEPT light curve in the usual manner"""
        figure = plt.figure()
        ax = plt.gca()

        dates = matplotlib.dates.date2num(self.data['DateTime'].astype(datetime))
        
        colors = ['Green','Red','Chocolate', 'Blue','SeaGreen','Tomato','SlateBlue','Orange',
        'Purple','Magenta','MediumVioletRed', 'Teal','Navy','Indigo']

        figure.delaxes(ax)
        axes = figure.add_axes([0.1, 0.15, 0.55, 0.8])

        for i,line in enumerate(self.header.values()):
            if i >= 1 and i <= 15:
                axes.plot_date(dates, self.data[line].ffill(), '-',
                     label= line[line.index('(')+1:line.index('V')+1] , color=colors[i/2-2], lw=0.5)
        
        axes.set_yscale("log",nonposy='mask')
        axes.set_title(title)
        axes.set_ylabel('1/(cm^2 s sr MeV)')
        axes.set_xlabel('UTC Time')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')

        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        figure.autofmt_xdate()
        plt.show()

        return figure


    @staticmethod
    def _parse_txt(filepath):
        """
        Parses a STEREO SEPT file from
        http://www2.physik.uni-kiel.de/stereo/index.php?doc=data

        """
        #Reading in Data only, using default i.e. 0 value for data_start keyword since all lines before data are commented
        data = ascii.read(filepath, delimiter = "\s") 

        #Header
        energy_levels = [45, 55, 65, 75, 85, 105, 125, 145, 165, 195, 225, 255, 295, 335, 375, 425]
        header = ['DateTime'] + \
        ['Bin {b:02d} ({x:3.1f} - {y:3.1f} keV) electron intensity'.format(b=b, x=x, y=y) for b, x, y in zip(range(2,17), \
            energy_levels[:-1], energy_levels[1:])] + ['Bin {b:02d} Uncertainty'.format(b=b) for b in range(2,17)]
       
        header.append('Accumulation time in seconds')
        data_modify = []

        #Storing data columns in recognizable variables
        year_col    = data['col2']
        day_of_year_col = data['col3']
        hour_col    = data['col4']
        minutes_col = data['col5']
        seconds_col = data['col6']

        #Converting first five columns into a single datetime.datetime column
        data_modify = [datetime.strptime('{0} {1:03d} {2} {3} {4}'.format(year, int(day), hour, minute, second), '%Y %j %H %M %S') \
                    for year, day, hour, minute, second in zip(year_col, day_of_year_col, hour_col, minutes_col, seconds_col)]


        #Removing separate datetime element columns
        data.remove_columns(['col{}'.format(i) for i in range(1,7)])
        #Adding the combined datetime column created above in data_modify
        data.add_column(Column(data = data_modify, name='col1'),0)

        #To add the column names in the astropy table object
        for elem, head_key in enumerate(header):
            data.rename_column(data.colnames[elem], head_key)        

        # Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build and none of the stable
        data = _to_pandas(data)
        data = data.replace(-9999.9,float('nan')) 

        return OrderedDict(enumerate(header)), data


class HETLightCurve(LightCurve):
    """
    STEREO HET LightCurve. Currently the LightCurve supports only Single File Load 

    Parameters
    ----------
    timerange: sunpy.time.TimeRange
        time range for which data is to be downloaded.
        Example value -  TimeRange('2006-12-01','2015-03-01')   
    
    stereo_spacecraft: string   
        Default value - ahead
        Possible values - ahead, behind    ( corresponding to spacecraft location )

    duration_of_average: astropy units quantity
        Default value - 15 * u.min
        Possible values - 1 * u.min, 15 * u.min, 1 * u.h, 12 * u.h, 1 * u.d        (corresponding to duration over which data is averaged)

    Examples
    --------

    .. plot::

    
        import os
        import sunpy.data.test
        filepath = sunpy.data.test.rootdir
        from sunpy import lightcurve as lc
        [header,data] = lc.HETLightCurve._parse_txt(os.path.join(filepath , 'het', 'AeH06Dec.1m.txt'))
        het = lc.HETLightCurve(data,header)
        het.peek()

    References
    ----------
    | http://www.srl.caltech.edu/STEREO/DATA/HET/

    """

    def peek(self, title="HET electron/proton Flux"):
        """Plots HET light curve in the usual manner"""

        figure = plt.figure()
        ax = plt.gca()

        if self.header.values()[1] == 'DateTime':
            dates = matplotlib.dates.date2num(self.data['DateTime'].astype(datetime))
        else:
            timerange_start = self.data['TimeRange'].apply(lambda col: col.start)
            dates = matplotlib.dates.date2num(timerange_start.astype(datetime))

        colors = ['Green','Red','Chocolate', 'Blue','SeaGreen','Tomato','SlateBlue','Orange',
        'Purple','Magenta','MediumVioletRed', 'Teal','Navy','Indigo']

        figure.delaxes(ax)
        axes = figure.add_axes([0.1, 0.15, 0.55, 0.8])

        for i,line in enumerate(self.header.values()):
            if i >= 2 and i%2 == 0:
                axes.plot_date(dates, self.data[line].ffill(), '-',
                label= line[:line.index('n')+2] + line[line.index(',')+2:line.index('V')+1], color=colors[i/2-2], lw=0.5)
        
        axes.set_yscale("log",nonposy='mask')
        axes.set_title(title)
        axes.set_ylabel('particles/(cm2-sr-sec-MeV)')
        axes.set_xlabel('UTC Time')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')

        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        figure.autofmt_xdate()
        plt.show()

        return figure


    @staticmethod
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
        header = [  'Electron flux, 0.7-1.4 MeV, particles/(cm2-sr-sec-MeV)', 
                    'Uncertainty (sigma) for 0.7-1.4 MeV electron flux',
                    'Electron flux, 1.4-2.8 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 1.4-2.8 MeV electron flux',
                    'Electron flux, 2.8-4.0 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 2.8-4.0 MeV electron flux',
                    'Proton flux, 13.6-15.1 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 13.6-15.1 MeV proton flux',
                    'Proton flux, 14.9-17.1 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 14.9-17.1 MeV proton flux',
                    'Proton flux, 17.0-19.3 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 17.0-19.3 MeV proton flux',
                    'Proton flux, 20.8-23.8 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 20.8-23.8 MeV proton flux',
                    'Proton flux, 23.8-26.4 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 23.8-26.4 MeV proton flux',
                    'Proton flux, 26.3-29.7 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 26.3-29.7 MeV proton flux',
                    'Proton flux, 29.5-33.4 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 29.5-33.4 MeV proton flux' ,
                    'Proton flux, 33.4-35.8 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 33.4-35.8 MeV proton flux' ,
                    'Proton flux, 35.5-40.5 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 35.5-40.5 MeV proton flux' ,
                    'Proton flux, 40.0-60.0 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 40.0-60.0 MeV proton flux' ,
                    'Proton flux, 60.0-100.0 MeV, particles/(cm2-sr-sec-MeV)',
                    'Uncertainty (sigma) for 60.0-100.0 MeV proton flux' ]
       

        data_modify = []

        #Storing data columns in recognizable variables (Common to all types of data)
        start_year_col  = data['col2']
        start_month_col = data['col3']
        start_date_col  = data['col4']
        start_time_col  = data['col5']

        #Adding Time Column based on type of data
        if type_of_data == '1min':
            header = ['Verse Number', 'DateTime'] + header 

            for i in range(len(data)): 
                date = datetime.strptime('{0} {1} {2:02d} {3:04d}'.format(start_year_col[i], start_month_col[i], start_date_col[i], \
                    start_time_col[i]), '%Y %b %d %H%M')
                data_modify.append(date)

            data.remove_columns(['col{}'.format(i) for i in range(2,6)])
        else:
            header = ['Verse Number','TimeRange'] + header
            
            #Storing data columns in recognizable variables
            end_year_col    = data['col6']
            end_month_col   = data['col7']
            end_date_col    = data['col8']
            end_time_col    = data['col9']


            for i in range(len(data)): 
                date1 = datetime.strptime('{0} {1} {2:02d} {3:04d}'.format(start_year_col[i], start_month_col[i], start_date_col[i], \
                    start_time_col[i]), '%Y %b %d %H%M')
                date2 = datetime.strptime('{0} {1} {2:02d} {3:04d}'.format(end_year_col[i], end_month_col[i], end_date_col[i], \
                    end_time_col[i]), '%Y %b %d %H%M' )
                data_modify.append(TimeRange(date1,date2))

            data.remove_columns(['col{}'.format(i) for i in range(2,10)])

        data.add_column(Column(data = data_modify, name='col'),1)

        #To add the column names in the astropy table object
        for elem, head_key in enumerate(header):
            data.rename_column(data.colnames[elem], head_key)        

        # Converting from astropy.table.Table to pandas.Dataframe
        # to_pandas() bound method is only available in the latest development build and none of the stable
        data = _to_pandas(data)
      
        return OrderedDict(enumerate(header)), data

