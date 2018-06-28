import textwrap
import warnings

from sunpy.physics.solar_rotation import *

warnings.warn(textwrap.dedent("""\
                  The module 'sunpy.physics.transforms.solar_rotation.py' is deprecated and will be
                  removed in the future. Use 'sunpy.physics.solar_rotation.py' instead"""),
              stacklevel=2)
