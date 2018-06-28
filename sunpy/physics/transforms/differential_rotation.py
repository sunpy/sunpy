import textwrap
import warnings

from sunpy.physics.differential_rotation import *

warnings.warn(textwrap.dedent("""\
                  The module 'sunpy.physics.transforms.differential_rotation.py' is deprecated and
                  will be removed in the future. Use 'sunpy.physics.differential_rotation.py' instead"""),
              stacklevel=2)
