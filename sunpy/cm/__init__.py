"""
The following colormaps are provided by this module:
"""
from sunpy.cm.cm import *  # noqa

for cmname in cmlist.keys():  # noqa
    __doc__ += "\n* '`~sunpy.cm.cm.{}`'\n".format(cmname)
