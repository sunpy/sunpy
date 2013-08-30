from __future__ import absolute_import

# for exposure to from sunpy.net.vso import *
from sunpy.net.vso.vso import VSOClient, InteractiveVSOClient

# for, e.g.,
# >>> from sunpy.net import vso
# >>> vso.search(...)
# >>> vso.get(...)
from sunpy.net.vso.vso import search, get

__all__ = ['VSOClient', 'InteractiveVSOClient']
