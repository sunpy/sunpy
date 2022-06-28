import pytest
from astropy.wcs import WCS
import matplotlib.pyplot as plt

from sunpy.tests.helpers import figure_test
from sunpy.visualization import draw
from sunpy.map.tests.conftest import aia171_test_map, heliographic_test_map


@figure_test
def test_draw_equator_aia171(aia171_test_map):
    fig = plt.figure()
    axes = fig.add_subplot(projection=aia171_test_map)
    aia171_test_map.plot()
    draw.equator(axes)


@figure_test
def test_draw_prime_meridian_aia171(aia171_test_map):
    fig = plt.figure()
    axes = fig.add_subplot(projection=aia171_test_map)
    aia171_test_map.plot()
    draw.prime_meridian(axes)


@figure_test
def test_heliographic_equator_prime_meridian(heliographic_test_map):
    fig = plt.figure()
    axes = fig.add_subplot(projection=heliographic_test_map)
    heliographic_test_map.plot()
    draw.equator(axes, color="blue")
    draw.prime_meridian(axes, color="red")


def test_prime_meridian_error():
    axes = plt.subplot(projection=WCS())
    with pytest.raises(AttributeError):
        draw.prime_meridian(axes)
