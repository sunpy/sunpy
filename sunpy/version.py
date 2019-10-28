# This file is for compatibility with astropy_helpers
version = 'unknown.dev'
try:
    from importlib_metadata import version as _version, PackageNotFoundError
    version = _version('sunpy')
except ImportError:
    from pkg_resources import get_distribution, DistributionNotFound
    try:
        version = get_distribution("sunpy").version
    except DistributionNotFound:
        pass
except PackageNotFoundError:
    pass
