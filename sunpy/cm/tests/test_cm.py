import sunpy.cm

import pytest
import sunpy.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def test_get_cmap():
    assert type(cm.get_cmap(name = 'sdoaia94')) == LinearSegmentedColormap


def test_show_colormaps():
    if cm.show_colormaps() is None:
        assert True
