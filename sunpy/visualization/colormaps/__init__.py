"""
The following colormaps are provided by this module.

.. plot::

    import matplotlib.pyplot as plt
    import sunpy.visualization.colormaps as cm
    cm.show_colormaps()


"""
from sunpy.visualization.colormaps.cm import *

for cmname in cmlist.keys():  # NOQA: F405
    __doc__ += f"\n* '{cmname}'\n"
