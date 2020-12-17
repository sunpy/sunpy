"""
Test Composite Map
"""

import pytest

import sunpy.data.test
import sunpy.map
from sunpy.tests.helpers import figure_test

testpath = sunpy.data.test.rootdir

# Ignore missing metadata warnings
pytestmark = [pytest.mark.filterwarnings('ignore:Missing metadata for observer'),
              pytest.mark.filterwarnings(r'ignore:Unable to treat `\.meta` as a FITS header')]


@pytest.fixture
def composite_test_map(aia171_test_map, hmi_test_map):
    return sunpy.map.Map(aia171_test_map, hmi_test_map, composite=True)


@figure_test
def test_plot_composite_map(composite_test_map):
    composite_test_map.plot()


@figure_test
def test_peek_composite_map(composite_test_map):
    composite_test_map.peek()


@figure_test
def test_plot_composite_map_linewidths(composite_test_map):
    composite_test_map.plot(linewidths=4)
