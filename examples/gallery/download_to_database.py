# -*- coding: utf-8 -*-
"""
====================================================
Download and store files in the database (using VSO)
====================================================

A simple demonstration of downloading and storing VSO data using the database
module.
"""


##############################################################################
# Start by importing the necessary modules.

from sunpy.database import Database
from sunpy.database.tables import display_entries
from sunpy.net import vso

##############################################################################
# Instantiate the VSO and database clients.

client = vso.VSOClient()
database = Database('sqlite:///sunpydata.sqlite')

##############################################################################
# Downloading data. We can download the required files and store them so
# that they can be easily accessed for future use. Put the VSO query for which
# you want to download files as the arguments to the fetch function. An
# optional argument path can be passed to specify the download location.

database.fetch(vso.attrs.Time('2011-05-08', '2011-05-08 00:00:05'),
	vso.attrs.Instrument('AIA'))

##############################################################################
# After downloading, The number of entries that will be added depends on the
# total number of FITS headers.

print(display_entries(database,
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# Closing the database session.

database.session.close()
