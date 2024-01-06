# Check if user has installed the image extras
from sunpy.util.sysinfo import warn_missing_deps as _warn_missing_deps
_warn_missing_deps('image')
