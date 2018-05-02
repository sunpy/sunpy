import warnings

from sunpy.util.exceptions import SunpyDeprecationWarning

deprecation_message = ("As of v0.8.0, the `sunpy.spectra` module is deprecated and will be "
                       "removed in a future version. This module is being moved to radiospectra "
                       "- http://docs.sunpy.org/projects/radiospectra/")
warnings.warn( deprecation_message, SunpyDeprecationWarning)
