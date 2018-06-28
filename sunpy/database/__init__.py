# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

"""
Overview
^^^^^^^^
The database package exports the following classes and exceptions:

    :classes:
        - Database
    :exceptions:
        - EntryAlreadyAddedError
        - NoSuchEntryError
        - EntryAlreadyStarredError
        - EntryAlreadyUnstarredError
        - EntryNotFoundError
        - TagAlreadyAssignedError
        - NoSuchTagError
        - NonRemovableTagError
    :functions:
        - disable_undo
"""

from __future__ import absolute_import

from sunpy.database.commands import NoSuchEntryError, NonRemovableTagError
from sunpy.database.database import (Database, NoSuchTagError, EntryNotFoundError,
                                     EntryAlreadyAddedError, TagAlreadyAssignedError,
                                     EntryAlreadyStarredError, EntryAlreadyUnstarredError,
                                     disable_undo, split_database)

__all__ = [
    'Database', 'EntryAlreadyAddedError', 'NoSuchEntryError', 'NoSuchTagError',
    'NonRemovableTagError', 'EntryAlreadyStarredError',
    'EntryAlreadyUnstarredError', 'EntryNotFoundError',
    'TagAlreadyAssignedError', 'disable_undo', 'split_database']
