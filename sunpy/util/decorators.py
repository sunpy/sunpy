# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Sundry function and class decorators."""

from __future__ import print_function

from astropy.utils.decorators import deprecated as a_deprecated
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.tests.helper import catch_warnings
import warnings

from sunpy.util.exceptions import SunpyDeprecationWarning


__all__ = ['deprecated']


def deprecated(since, **kwargs):
    """
    Used to mark a function or class as deprecated.

    For accepted keyword parameters see `~astropy.utils.decorators.deprecated`.
    """

    with catch_warnings(AstropyDeprecationWarning) as w:
        a_deprecated(since, **kwargs)
        warnings.warn(w.message, SunpyDeprecationWarning)
