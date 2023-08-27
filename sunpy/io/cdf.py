from sunpy.util.exceptions import warn_deprecated
from . import _cdf
from ._cdf import *  # NOQA

__doc__ = _cdf.__doc__
__all__ = _cdf.__all__

warn_deprecated("The `sunpy.io.cdf` module is deprecated, as it was designed "
                "for internal use.")
