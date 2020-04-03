"""SunPy Maps"""
import sys
import textwrap

from sunpy.map.mapbase import GenericMap  # isort:skip

from sunpy.map import sources
from sunpy.map.header_helper import *
from sunpy.map.map_factory import Map
from sunpy.map.maputils import *

from .compositemap import CompositeMap
from .mapsequence import MapSequence


# If the documentation is being built using Sphinx, append output HTML to the docstrings for the
#   GenericMap and MapSequence `quicklook()` methods.
# This code lives here instead of the respective source files because the source files need to be
#   fully imported before the following code will successfully run.
# This code detects a Sphinx build via the simplistic approach of checking if Sphinx has been
#   imported, and thus there is the potential for false positives.
if 'sphinx' in sys.modules:
    import sunpy.data.sample
    html_string = textwrap.indent(Map(sunpy.data.sample.AIA_171_IMAGE)._repr_html_(), ' ' * 12)
    GenericMap.quicklook.__doc__ += f"""\

        (which will open the following content in the default web browser)

        .. raw:: html

            <div style="border:1px solid black">{html_string}</div>

        """

    seq = Map(sunpy.data.sample.HMI_LOS_IMAGE,
              sunpy.data.sample.AIA_1600_IMAGE,
              sunpy.data.sample.EIT_195_IMAGE,
              sequence=True)
    html_string = textwrap.indent(seq._repr_html_(), ' ' * 12)
    MapSequence.quicklook.__doc__ += f"""\

        (which will open the following content in the default web browser)

        .. raw:: html

            <div style="border:1px solid black">{html_string}</div>

        """
