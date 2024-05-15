"""
This module provides a generic FileHeader object for the readers.
"""
from collections import OrderedDict

__all__ = ['FileHeader']


class FileHeader(OrderedDict):
    """
    FileHeader is designed to provide a consistent interface to all other sunpy
    classes that expect a generic file.

    Open read all file types should format their header into a
    FileHeader
    """

    def __init__(self, *args, **kwargs):
        OrderedDict.__init__(self, *args, **kwargs)
