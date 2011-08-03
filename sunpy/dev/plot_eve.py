#-*- coding:utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
# Written: 2011/03/25
#
# <License info will go here...>
#
'''Grab the EVE GOES Proxy data and plot it in a standard GOES plot format'''

import csv
import urllib
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, date, time	
	#
# Notes:
#
# Steven (2011/03/30)
# This function downloads a file called latest_eve.txt and does not clean up after itself. 
# Should really read the file into memory rather than creating a file if possible.
#
#

def plot_eve():
	'''Grab the EVE GOES Proxy data and plot it in a standard GOES plot format'''
	
	url = 'http://lasp.colorado.edu/eve/data_access/quicklook/quicklook_data/L0CS/LATEST_EVE_L0CS_DIODES_1m.txt'
	
	urllib.urlretrieve(url, filename='latest_eve.txt')
	
	reader = csv.reader(open("latest_eve.txt", "rb"), delimiter=' ', skipinitialspace=True)
	
	i = 0
	
	t = []
	xrsb = []
	xrsa = []
	
	for row in reader:
		#skip header rows
		if row[0][0] != ';':
			#read the date line
			if i == 0:
				d = date(int(row[0]), int(row[2]), int(row[3]))
				#print date
			else:
				t.append(time(int(row[0][0:2]), int(row[0][2:4])))
				xrsb.append(float(row[1]))
				xrsa.append(float(row[2])) 
				i = i + 1
	ts = [datetime.combine(d, s) for s in t]
	dates = matplotlib.dates.date2num(ts)
	goes_plot(t, xrsa, xrsb, title="EVE GOES Proxy Xray Flux (1 minute data)")
	
#save it as a png
#matplotlib.pyplot.savefig('test.png')
#eveplot(0)

def plot_latest_goes():
	"""Grab the latest GOES data and plot it in a standard GOES plot format"""
	url = 'http://www.swpc.noaa.gov/ftpdir/lists/xray/Gp_xr_1m.txt'
	urllib.urlretrieve(url, filename='latest_goes.txt')
	
	reader = csv.reader(open("latest_goes.txt", "rb"), delimiter=' ', skipinitialspace=True)
		
	t = []
	xrsb = []
	xrsa = []
	
	for row in reader:
		#skip header rows
		if (row[0][0] != ':') and (row[0][0] != '#'):
			#read the date line
			d = date(int(row[0]), int(row[1]), int(row[2]))
			t.append(time(int(row[3][0:2]), int(row[3][2:4])))
			xrsb.append(float(row[6]))
			xrsa.append(float(row[7])) 
				
	ts = [datetime.combine(d, s) for s in t]	
	dates = matplotlib.dates.date2num(ts)
	
	goes_plot(dates, xrsa, xrsa)
	
def goes_plot(t, xrsa, xrsb, title=""):
	"""Create a standard GOES plot"""
	#now plot!
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot_date(t, xrsa, '-', label='0.5--4.0 $\AA$', color='blue', lw=2)
	ax.plot_date(t, xrsb, '-', label='1.0--8.0 $\AA$', color='red', lw=2)
	ax.set_yscale("log")
	ax.set_ylim(1e-9, 1e-2)
	ax.set_title(title)
	ax.set_ylabel('Watts m$^{-2}$')
	#ax.set_xlabel(d)
	
	ax2 = ax.twinx()
	ax2.set_yscale("log")
	ax2.set_ylim(1e-9, 1e-2)
	ax2.set_yticks((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
	ax2.set_yticklabels((' ', 'A', 'B', 'C', 'M', 'X', ' '))
	
	ax.yaxis.grid(True, 'major')
	ax.xaxis.grid(False, 'major')
	ax.legend()
	
	formatter = matplotlib.dates.DateFormatter('%H:%M')
	ax.xaxis.set_major_formatter(formatter)
	
	ax.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
	fig.autofmt_xdate()
	plt.show()

