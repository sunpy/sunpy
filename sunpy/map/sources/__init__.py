"""
Datasource-specific classes

This is where datasource specific logic is implemented. Each mission should
have its own file with one or more classes defined. Typically, these classes
will be subclasses of the :mod`sunpy.map.Map` class.
"""

from sunpy.map.map_factory import Map
from .adapt import *
from .gong import *
from .hinode import *
from .iris import *
from .mlso import *
from .proba2 import *
from .psp import *
from .punch import *
from .rhessi import *
from .sdo import *
from .soho import *
from .solo import *
from .source_type import *
from .stereo import *
from .suit import *
from .suvi import *
from .trace import *
from .yohkoh import *
