# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains errors/exceptions and warnings of general use for
sunpy. Exceptions that are specific to a given subpackage should *not*
be here, but rather in the particular subpackage.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


class SunpyWarning(Warning):
    """
    The base warning class from which all Sunpy warnings should inherit.
    """


class SunpyUserWarning(UserWarning, SunpyWarning):
    """
    The primary warning class for Sunpy.

    Use this if you do not need a specific sub-class.
    """


class SunpyDeprecationWarning(SunpyWarning):
    """
    A warning class to indicate a deprecated feature.
    """
