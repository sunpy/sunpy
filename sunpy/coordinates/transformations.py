from _transforms import *  # NOQA

from sunpy.util.exceptions import warn_deprecated
from . import _transformations

__doc__ = _transformations.__doc__
__all__ = _transformations.__all__

warn_deprecated("The `sunpy.coordinates.transformations` module is deprecated, as it was designed "
                "for internal use.")
