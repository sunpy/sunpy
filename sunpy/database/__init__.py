"""
Overview
--------
The database package exports the following classes, functions and exceptions:

    :classes:
        - Database
        - DatabaseEntry
        - FitsHeaderEntry
        - Tag
    :functions:
        - entries_from_query_result
        - entries_from_path
    :exceptions:
        - EntryAlreadyAddedError
        - NoSuchEntryError
        - EntryAlreadyStarredError
        - EntryAlreadyUnstarredError
"""
from __future__ import absolute_import

from sunpy.database.database import Database, EntryAlreadyAddedError,\
    EntryAlreadyStarredError, EntryAlreadyUnstarredError, NoSuchTagError
from sunpy.database.tables import DatabaseEntry, FitsHeaderEntry, Tag,\
    entries_from_query_result, entries_from_path
from sunpy.database.commands import NoSuchEntryError

__all__ = [
    'Database', 'EntryAlreadyAddedError', 'NoSuchEntryError', 'NoSuchTagError',
    'EntryAlreadyStarredError', 'EntryAlreadyUnstarredError',
    'DatabaseEntry', 'FitsHeaderEntry', 'Tag', 'entries_from_query_result',
    'entries_from_path']
