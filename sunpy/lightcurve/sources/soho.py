# -*- coding: utf-8 -*-
"""SOHO LightCurve sources subclass definitions"""

from __future__ import absolute_import


__authors__ = ["Ankit Kumar"]
__email__   = "ankitkmr.iitk@gmail.com"

# This module was developed with funding from 
# Google Summer of Code 2015


from datetime import timedelta,datetime,time

import pandas as pd
from pandas import DataFrame
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates

from astropy.io import ascii
from astropy.table import Table, Column
from astropy.utils import OrderedDict
from astropy.table.column import MaskedColumn

from sunpy.time import TimeRange
from sunpy.lightcurve import LightCurve



__all__ = ['ERNELightCurve']


#This code for _to_pandas is repeated here and in stereo.py for performance reasons.
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
		timerange ( sunpy.time.TimeRange ), 
		specie ( string )

	POSSIBLE KEYWORD VALUES:-  
		specie: 'proton' or 'alpha'

	** Currently the LightCurve supports only Single File Load **
	
	Examples
	--------

	>>> import os
	>>> import sunpy.data.test
	>>> filepath = sunpy.data.test.rootdir
	>>> from sunpy import lightcurve as lc
	>>> erne = lc.ERNELightCurve._parse_txt(os.path.join(filepath , 'erne/cr1907a.txt'))
	>>> erne = lc.ERNELightCurve(erne[1],erne[0])
	>>> erne.peek()

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

		for i,line in enumerate(self.header.values()):
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

		and returns header as a list and ERNE data as pandas dataframe.
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


		n = len(header)
		return [OrderedDict(zip(range(n), header)), data]

