from sunpy.util.exceptions import warn_deprecated
from . import _transformations
from ._transformations import *  # NOQA

__doc__ = _transformations.__doc__
__all__ = _transformations.__all__

warn_deprecated("The `sunpy.coordinates.transformations` module is deprecated, as it was designed "
                "for internal use. The context managers `transform_with_sun_center()` and "
                "`propagate_with_solar_surface()` continue to be available under "
                "`sunpy.coordinates`.")
