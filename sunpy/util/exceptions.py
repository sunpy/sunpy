"""
This module provides errors/exceptions and warnings of general use for SunPy.

Exceptions that are specific to a given package should **not** be here,
but rather in the particular package.
"""
import sys

from astropy.utils.exceptions import AstropyWarning

__all__ = ["SunpyWarning", "SunpyUserWarning", "SunpyDeprecationWarning",
           "SunpyPendingDeprecationWarning"]


class SunpyWarning(AstropyWarning):
    """
    The base warning class from which all Sunpy warnings should inherit.

    Any warning inheriting from this class is handled by the Sunpy
    logger. This warning should not be issued in normal code. Use
    "SunpyUserWarning" instead or a specific sub-class.
    """


class SunpyUserWarning(UserWarning, SunpyWarning):
    """
    The primary warning class for Sunpy.

    Use this if you do not need a specific type of warning.
    """


class SunpyDeprecationWarning(FutureWarning, SunpyWarning):
    """
    A warning class to indicate a deprecated feature.
    """


class SunpyPendingDeprecationWarning(PendingDeprecationWarning, SunpyWarning):
    """
    A warning class to indicate a soon-to-be deprecated feature.
    """
