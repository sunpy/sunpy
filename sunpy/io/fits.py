from sunpy.util.exceptions import warn_deprecated
from . import _fits
from ._fits import *  # NoQA

__doc__ = _fits.__doc__
__all__ = _fits.__all__

warn_deprecated("The `sunpy.io.fits` module is deprecated, as it was designed "
                "for internal use. Use the `astropy.fits.io` module instead "
                "for more generic functionality to read FITS files.")
