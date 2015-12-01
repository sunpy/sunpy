
'''
Project	:	HELCATS

Name	:	cmes

Purpose	:	Read the CME catalog

Explanation:	Input the CDAW CME catalog as a text file and generate a dataframe to house the columns of interest, output as a CSV.

Use	:	$ python cmes.py

Inputs	:	filename - of CDAW text file "univ_all.txt" sourced from http://cdaw.gsfc.nasa.gov/CME_list/UNIVERSAL/text_ver/univ_all.txt

Outputs	:	CSV file of CDAW catalog "cdaw_catalog.csv"

Keywords:	

Calls	:	os, numpy, matplotlib, pandas

Written	:	Jason P Byrne, STFC/RAL Space, Dec 2015 (jason.byrne@stfc.ac.uk)

Revisions:
2015-01-01 JPB : Working in accordance with sunpy project guidelines.

'''

import os
import numpy as np
import pandas as pd
import urllib2

def cdaw():
	###
	# Loop over lines in text file, creating a list of dictionaries
	# Open the CDAW CME text file for reading
	#f=open(os.path.join(config.cdaw_path,'univ_all.txt'),'r')
	
	try:
		f = urllib2.urlopen("http://cdaw.gsfc.nasa.gov/CME_list/UNIVERSAL/text_ver/univ_all.txt")

		# Read in the header lines first
		for i in range(4):
			f.readline()
		
		# initialise the data list
		cols = ['date','time','cpa','width','lin_speed','quad_speed_init',\
			'quad_speed_final','quad_speed_20','accel','mass','kin_energy','mpa','remarks']
		source = {k:list() for k in cols}
		for line in f:
			
			line = line.replace('Halo','360').strip()
			columns = line.split()
			columns = [col.strip('*') if ("--" not in col) and ("**" not in col) else np.NaN for col in columns]
			if len(columns)<13:
				columns.append(' ')
			for i,(k,val) in enumerate(zip(cols,columns)):
				if i==12:
					source[k].append(' '.join(columns[i:]))
				else:
					source[k].append(val)
				
		f.close()
	

		#Creating a pandas dataframe
		return pd.DataFrame(source,columns=cols)

	except:
		raise


