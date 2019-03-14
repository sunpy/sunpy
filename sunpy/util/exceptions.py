# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains errors/exceptions and warnings of general use for
sunpy. Exceptions that are specific to a given subpackage should *not*
be here, but rather in the particular subpackage.
"""
import sys

from astropy.utils.exceptions import AstropyWarning

__all__ = ["SunpyWarning", "SunpyUserWarning", "SunpyDeprecationWarning", "SunpyPendingDeprecationWarning"]



class SunpyWarning(AstropyWarning):
    """
    The base warning class from which all Sunpy warnings should inherit.

    Any warning inheriting from this class is handled by the Sunpy logger.
    Should likely not be used.
    """


class SunpyUserWarning(UserWarning, SunpyWarning):
    """
    The primary warning class for Sunpy.

    Use this if you do not need a specific sub-class.
    """


# For PEP 565 (https://www.python.org/dev/peps/pep-0565/) compliance.
DeprecationClass = DeprecationWarning if sys.version_info >= (3, 7) else FutureWarning


class SunpyDeprecationWarning(DeprecationClass, SunpyWarning):
    """
    A warning class to indicate a deprecated feature.
    """


class SunpyPendingDeprecationWarning(PendingDeprecationWarning, SunpyWarning):
    """
    A warning class to indicate a soon-to-be deprecated feature.
    """
