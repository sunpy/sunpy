import pytest
import matplotlib.pyplot as plt

import sunpy.cm as cm


def test_get_cmap():
    colormap_list = cm.cmlist
    for cmap in colormap_list.keys():
        assert colormap_list[cmap] == plt.get_cmap(cmap)


def test_show_colormaps():
    if cm.show_colormaps() is None:
        assert True
