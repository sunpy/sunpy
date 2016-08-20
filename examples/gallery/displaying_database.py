# -*- coding: utf-8 -*-
"""
=========================================
Displaying the contents of a database.
=========================================

A simple demonstration showing the various ways to displays the contents of a
database.
"""


##############################################################################
# Import the necessary modules.

from sunpy.database import Database
from sunpy.database.tables import display_entries
from sunpy.net import vso

##############################################################################
# Instantiate the VSO and database clients.

client = vso.VSOClient()
database = Database('sqlite:///sunpydata.sqlite')

##############################################################################
# Now to populate the database.

vso_query_result = client.query(
	vso.attrs.Time('2011-05-08', '2011-05-08 00:00:05'),
	vso.attrs.Instrument('AIA'))
database.add_from_vso_query_result(vso_query_result)

##############################################################################
# Using display_entries from sunpy.database.tables (old method). If no column
# names are specified, then by default all columns are displayed. Here a list
# containing the desired column names are passed as an argument.

print(display_entries(database,
	['id', 'observation_time_start', 'observation_time_end',
	'instrument', 'wavemin', 'wavemax']))

##############################################################################
# The contents can also be printed simply by using print statement on the
# database object. Note that the resultant table may be truncated according to
# the terminal size.

print(database)

##############################################################################
# We can also display the database entries in a nicely formatted HTML page.

database.show_in_browser(jsviewer=True)

##############################################################################
# Closing the session.

database.session.close()
