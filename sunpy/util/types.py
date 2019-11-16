"""
A collection of base classes to use to shortcut having to import things.

This is used to provide a base class which can be imported into two locations
to identify a type, so that either location doesn't have to import the other.
"""


class DatabaseEntryType:
    """
    A base class for `sunpy.database.tables.DatabaseEntry`.
    This class should not be used directly.
    """
    # This is currently used to prevent `sunpy.map` having to import
    # `sunpy.database`.
