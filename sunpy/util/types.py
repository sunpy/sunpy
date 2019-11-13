"""
A collection of ABCs to use to shortcut having to import things.
"""
import abc


class DatabaseEntryType(metaclass=abc.ABCMeta):
    """
    A base class for `sunpy.database.tables.DatabaseEntry`.
    """
