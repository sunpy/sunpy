import pytest

import matplotlib.pyplot as plt

import sunpy.cm as cm

# Checks that colormaps are imported by MPL
def test_get_cmap():
    for cmap in cm.cmlist.keys():
        assert cm.cmlist[cmap] == plt.get_cmap(cmap)
