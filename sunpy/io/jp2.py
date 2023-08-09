from sunpy.util.exceptions import warn_deprecated
from . import _jp2
from ._jp2 import *  # NOQA

__doc__ = _jp2.__doc__
__all__ = _jp2.__all__

warn_deprecated("The `sunpy.io.jp2` module is deprecated, as it was designed "
                "for internal use.")
