# This file is for compatibility with astropy_helpers
from pkg_resources import DistributionNotFound, get_distribution

try:
    version = get_distribution("sunpy").version
except DistributionNotFound:
    version = "unknown.dev"
