"""
Test draw
"""
import pytest

from sunpy.tests.helpers import figure_test
from sunpy.visualization import draw
from sunpy.map.tests.conftest import aia171_test_map


@figure_test
def test_draw_equator_aia171(aia171_test_map):
    aia171_test_map.plot()
    draw.draw_equator(aia171_test_map._check_axes(None))


@figure_test
def test_draw_prime_meridian_aia171(aia171_test_map):
    aia171_test_map.plot()
    draw.draw_prime_meridian(aia171_test_map._check_axes(None))
