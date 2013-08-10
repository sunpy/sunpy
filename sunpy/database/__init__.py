"""
Overview
--------
The database package exports the following classes, functions and exceptions:

    :classes:
        - Database
    :functions:
        - display_entries
    :exceptions:
        - EntryAlreadyAddedError
        - NoSuchEntryError
        - EntryAlreadyStarredError
        - EntryAlreadyUnstarredError
        - EntryNotFoundError
        - TagAlreadyAssignedError
        - NoSuchTagError
"""
from __future__ import absolute_import

from sunpy.database.database import Database, EntryAlreadyAddedError,\
    EntryAlreadyStarredError, EntryAlreadyUnstarredError, NoSuchTagError,\
    EntryNotFoundError, TagAlreadyAssignedError
from sunpy.database.commands import NoSuchEntryError

__all__ = [
    'Database', 'EntryAlreadyAddedError', 'NoSuchEntryError', 'NoSuchTagError',
    'EntryAlreadyStarredError', 'EntryAlreadyUnstarredError',
    'EntryNotFoundError', 'TagAlreadyAssignedError', 'display_entries']
