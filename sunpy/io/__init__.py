"""
This subpackage contains a series of file readers and only for internal use.

None of this code is intended to be used directly by a user and is not part of the sunpy public API
and is not added to the API reference pages.
"""
from sunpy.util.exceptions import warn_deprecated

warn_deprecated('The toplevel space of the io subpackage was never intended for public use and will be removed in the future.')
from sunpy.io._file_tools import read_file, write_file

__all__ = ['read_file', 'write_file']
