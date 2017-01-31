# -*- coding: utf-8 -*-
"""
=========================================
Caching in the database
=========================================

A simple demonstration showing the database module's caching mechanism
so that the same files aren't redownloaded.
"""


##############################################################################
# Import the necessary modules.

from sunpy.database import Database
from sunpy.database.tables import display_entries
from sunpy.net import vso

##############################################################################
# Instantiate the VSO and database clients

database = Database('sqlite:///sunpydata.sqlite')


##############################################################################
# Download data for a VSO query

database.fetch(vso.attrs.Time('2011-05-08', '2011-05-08 00:00:05'),
	vso.attrs.Instrument('AIA'))

##############################################################################
# After downloading, The number of entries that will be added depends on the
# total number of FITS headers.

print(display_entries(database,
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# Clever Fetching
# The fetch() method queries the database if the given query
# has already been used once to add entries using the download()/fetch()
# methods. Otherwise, the given query is used to download and add new data via
# the download() method. This means that if you are not sure whether the
# required data already exists in the database, use fetch() to save bandwidth.

print("Number of entries in database : " + str(len(database)))

##############################################################################

entries = database.fetch(vso.attrs.Time('2011-05-08', '2011-05-08 00:00:05'),
	vso.attrs.Instrument('AIA'))

##############################################################################

print("Number of entries returned by fetch : " + str(len(entries)))
print("Number of entries in database : " + str(len(database)))

##############################################################################
# We can see that no new entries were added to the database as we had
# previously downloaded the files of the same query. Now we use the fetch()
# method to download files for a query for which we haven't downloaded files.

entries = database.fetch(vso.attrs.Time('2012-05-08', '2012-05-08 00:00:05'),
	vso.attrs.Instrument('AIA'))
print(entries is None)
print("Number of entries in database : " + str(len(database)))
print(display_entries(database,
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# As this download for this query was not done before, so no entries were
# returned. We can see that new files were downloaded and new entries were
# added to the database. If you are not sure whether the required data already
# exists in the database, use the `~sunpy.database.Database.fetch` method !!!
# Finally lets close the database session.

database.session.close()
