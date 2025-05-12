"""
Test Composite Map
"""
import matplotlib as mpl
import numpy as np
import pytest

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

import sunpy.data.test
import sunpy.map
from sunpy.tests.helpers import figure_test

# Ignore missing metadata warnings
pytestmark = [pytest.mark.filterwarnings('ignore:Missing metadata for observer'),
              pytest.mark.filterwarnings(r'ignore:Unable to treat `\.meta` as a FITS header')]


@pytest.fixture
def composite_test_map(adjusted_test_maps):
    aia171_test_map, hmi_test_map = adjusted_test_maps
    return sunpy.map.Map(aia171_test_map, hmi_test_map, composite=True)


@pytest.fixture
def composite_test_map_prealigned(adjusted_test_maps):
    aia171_test_map, hmi_test_map = adjusted_test_maps
    hmi_test_map = hmi_test_map.reproject_to(aia171_test_map.wcs)
    return sunpy.map.Map(aia171_test_map, hmi_test_map, composite=True)


def test_type_of_arguments_composite_map(composite_test_map):
    with pytest.raises(ValueError, match="CompositeMap expects pre-constructed map objects.") as excinfo:
        sunpy.map.CompositeMap(23, composite=True)
    assert str(excinfo.value) == 'CompositeMap expects pre-constructed map objects.'


@figure_test
def test_autoalign_needed(composite_test_map):
    # The overplot will use pcolormesh if the WCSs do not match
    images = composite_test_map.plot()
    assert isinstance(images[0], mpl.image.AxesImage)
    assert isinstance(images[1], mpl.collections.QuadMesh)


@figure_test
def test_autoalign_not_needed(composite_test_map_prealigned):
    # The overplot will use imshow if the WCSs match
    images = composite_test_map_prealigned.plot()
    assert isinstance(images[0], mpl.image.AxesImage)
    assert isinstance(images[1], mpl.image.AxesImage)


@figure_test
def test_plot_composite_map(composite_test_map):
    composite_test_map.plot()


@figure_test
def test_plot_composite_map_contours(composite_test_map):
    composite_test_map.set_levels(1, np.arange(-75, 76, 25) << u.percent)
    composite_test_map.plot()


@figure_test
def test_plot_composite_map_linewidths(composite_test_map):
    composite_test_map.set_levels(1, np.arange(-75, 76, 25) << u.percent)
    composite_test_map.plot(linewidths=0.5)


@figure_test
def test_plot_composite_map_linestyles(composite_test_map):
    composite_test_map.set_levels(1, np.arange(-75, 76, 25) << u.percent)
    composite_test_map.plot(linestyles='--')


@figure_test
def test_plot_composite_map_colors(composite_test_map):
    composite_test_map.set_levels(1, np.arange(-75, 76, 25) << u.percent)
    composite_test_map.plot(colors='red')


def test_plot_composite_map_mplkwargs(composite_test_map):
    composite_test_map.set_levels(1, np.arange(-75, 76, 25) << u.percent)
    with pytest.raises(TypeError) as e:
        composite_test_map.plot(linestyles='--', unused_a=1, unused_b=2)
    assert 'plot() got unexpected keyword arguments' in str(e.value)
    assert 'unused_a' in str(e.value)
    assert 'unused_b' in str(e.value)
    assert 'linestyles' not in str(e.value)


def test_remove_composite_map(composite_test_map):
    composite_test_map.remove_map(0)
    with pytest.raises(IndexError):
        composite_test_map.get_map(1)


def test_get_composite_map(composite_test_map, aia171_test_map, hmi_test_map):
    assert composite_test_map.get_map(0) == aia171_test_map
    assert composite_test_map.get_map(1) == hmi_test_map


def test_get_alpha_composite_map(composite_test_map, aia171_test_map, hmi_test_map):
    assert composite_test_map.get_alpha() == [aia171_test_map.alpha, hmi_test_map.alpha]


def test_get_alpha_with_index_composite_map(composite_test_map, aia171_test_map, hmi_test_map):
    assert composite_test_map.get_alpha(0) == aia171_test_map.alpha
    assert composite_test_map.get_alpha(1) == hmi_test_map.alpha


def test_get_levels_composite_map(composite_test_map, aia171_test_map, hmi_test_map):
    assert composite_test_map.get_levels() == [aia171_test_map.levels, hmi_test_map.levels]


def test_get_levels_with_index_composite_map(composite_test_map, aia171_test_map, hmi_test_map):
    assert composite_test_map.get_levels(0) == aia171_test_map.levels
    assert composite_test_map.get_levels(1) == hmi_test_map.levels


@figure_test
def test_set_alpha_composite_map(composite_test_map):
    composite_test_map.set_alpha(1, 0.5)
    composite_test_map.plot()


@pytest.mark.parametrize(('index', 'alpha'), [(0, 5.0), (1, -3.0)])
def test_set_alpha_out_of_range_composite_map(composite_test_map, index, alpha):
    with pytest.raises(Exception, match="Alpha value must be between 0 and 1") as excinfo:
        composite_test_map.set_alpha(index, alpha)
    assert str(excinfo.value) == 'Alpha value must be between 0 and 1.'


def test_set_levels_percent(composite_test_map):
    numbers = np.arange(10, 100, 10)
    composite_test_map.set_levels(0, numbers)
    np.testing.assert_allclose(composite_test_map.get_levels(0), numbers)

    implicit_percentage = np.arange(10, 100, 10)
    composite_test_map.set_levels(0, implicit_percentage, percent=True)
    assert_quantity_allclose(composite_test_map.get_levels(0), implicit_percentage << u.percent)


@figure_test
def test_peek_composite_map(composite_test_map):
    composite_test_map.peek()
