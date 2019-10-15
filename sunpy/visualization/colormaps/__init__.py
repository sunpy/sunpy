"""
The following colormaps are provided by this module:
"""
from sunpy.visualization.colormaps.cm import *  # noqa

for cmname in cmlist.keys():  # noqa
    __doc__ += f"\n* '`~sunpy.visualisation.colormaps.{cmname}`'\n"
