# -*- coding: utf-8 -*-
"""
=========================================
The Database Module
=========================================

In this example we demonstrate some of the features of the Database module.
"""

##############################################################################
# :ref:`database_guide` gives a detailed and comprehensive guide on how
# to use the Database package.

##############################################################################
# Let's start by importing the necessary modules.

from sunpy.database import Database
from sunpy.database.tables import display_entries
from sunpy.net import vso

##############################################################################
# First we need to instantiate a `~sunpy.database.Database` object. Note that
# the argument can be any URL supported by SQLAlchemy.

database = Database('sqlite:///sunpydata.sqlite')

##############################################################################
# Adding Entries from a VSO query

vso_query_result = vso.VSOClient().query(
	vso.attrs.Time('2011-05-08', '2011-05-08 00:00:05'),
	vso.attrs.Instrument('AIA'))
database.add_from_vso_query_result(vso_query_result)

##############################################################################
# Now to see what is in the database, we use the display_entries function to
# display the desired parameters in a neat way.

print(display_entries(database,
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# Downloading data. We can download the required files and store them so
# that they can be easily accessed for future use. An optional argument path
# can be passed to specify the download location.

database.download(vso.attrs.Time('2011-05-08', '2011-05-08 00:00:05'),
	vso.attrs.Instrument('AIA'))

##############################################################################
# After downloading, The number of entries that will be added depends on the
# total number of FITS headers.

print(display_entries(database,
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# Clever Fetching.
# The fetch() method queries the database if the given query
# has already been used once to add entries using the download() method.
# Otherwise, the given query is used to download and add new data via the
# download() method. This means that if you are not sure whether the required
# data already exists in the database, use fetch() to save bandwidth.

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
# exists in the database, use the fetch() method !!!

##############################################################################
# Changing an entry's atrributes. We can use the Database.edit() method to
# change any attribute of a database entry. Here we edit the
# wavemax attribute ( if it is less than 20.0, then we change it to 100.0;

for database_entry in database:
	if database_entry.wavemax < 20.0:
		database.edit(database_entry,
			wavemax=100.0)

print(display_entries(database,
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# Starring Entries. We can star(italic) database entries which are of interest
# to us. Here we star all those entries with wavemin greater than 20.

for database_entry in database:
	if database_entry.wavemin > 20:
		database.star(database_entry)

print(display_entries(
	filter(lambda entry: entry.starred, database),
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# Adding tags to entries. By tagging entries we can add more information to
# each entry in the form of keywords. The tags(bold) property of a database
# object holds all the tags saved in the database. Here we tag all entries
# whose observation_time_start was in year 2011.

for database_entry in database:
	if database_entry.observation_time_start.year == 2011:
		database.tag(database_entry, 'eleven')

print(database.tags)

sample_tag = database.tags[0]
print(display_entries(
	filter(lambda entry: sample_tag in entry.tags, database),
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# Undo and redo operations. These features are very useful as they can undo or
# redo the last n commands. Note : If we work on a Database object in an 
# interactive Python session and quit this session, the undo and redo history
# are lost. Here we undo and redo the last 2 commands.

database.undo(2)
database.redo(2)

##############################################################################
# Removing entries. We use the database.remove() method to remove a database
# entry. Here we remmove all those entries whose observation_time_start is
# in year 2012.

for database_entry in database:
	if database_entry.observation_time_start.year == 2012:
		database.remove(database_entry)
print(display_entries(database,
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# Emptying the database. The clear() method can be used to remove all the
# entries in the database.

database.clear()
print(len(database))

##############################################################################
# Committing changes to the database. While working in an interactive python
# session, the commit() method must be called to save the changes to the
# database. Otherwise all changes are lost on quitting the session.

database.commit()