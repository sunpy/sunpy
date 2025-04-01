"""
This subpackage contains a series of file readers and only for internal use.

None of this code is intended to be used directly by a user and is not part of the sunpy public API
and is not added to the API reference pages.
"""
from sunpy.io._file_tools import read_file as _read_file
from sunpy.io._file_tools import write_file as _write_file
from sunpy.util.decorators import deprecated


@deprecated("6.0", message="The toplevel space of the io subpackage was never intended for public use and will be removed in the future.")
def read_file(*args, **kwargs):
    return _read_file(*args, **kwargs)

@deprecated("6.0", message="The toplevel space of the io subpackage was never intended for public use and will be removed in the future.")
def write_file(*args, **kwargs):
    return _write_file(*args, **kwargs)

__all__ = ['read_file', 'write_file']
