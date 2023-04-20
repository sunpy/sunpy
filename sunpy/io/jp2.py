"""
This module provides a JPEG 2000 file reader.
"""
from sunpy.util.exceptions import warn_deprecated
from . import _jp2

__doc__ = _jp2.__doc__
__all__ = _jp2.__all__ # NOQA: PLE0605

warn_deprecated("The `sunpy.io.jp2` module is deprecated, as it was designed "
                "for internal use.")
