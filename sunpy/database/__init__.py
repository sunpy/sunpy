# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).
from sunpy.database.commands import NonRemovableTagError, NoSuchEntryError
from sunpy.database.database import (
    Database,
    EntryAlreadyAddedError,
    EntryAlreadyStarredError,
    EntryAlreadyUnstarredError,
    EntryNotFoundError,
    NoSuchTagError,
    PartialFetchError,
    TagAlreadyAssignedError,
    disable_undo,
    split_database,
)
from sunpy.util.exceptions import warn_deprecated

warn_deprecated(
    "The sunpy.database module is no longer actively maintained and has a number of outstanding issues. "
    "It will be removed in sunpy 6.0. "
    "If you are using sunpy.database and would like to see a replacement, please join the discussion thread at https://community.openastronomy.org/t/deprecating-sunpy-database/495"
)

__all__ = [
    'Database', 'EntryAlreadyAddedError', 'NoSuchEntryError', 'NoSuchTagError',
    'NonRemovableTagError', 'EntryAlreadyStarredError',
    'EntryAlreadyUnstarredError', 'EntryNotFoundError',
    'TagAlreadyAssignedError', 'PartialFetchError', 'disable_undo', 'split_database']
