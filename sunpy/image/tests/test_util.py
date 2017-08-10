from __future__ import absolute_import, division, print_function
import pytest

import numpy as np
from numpy.testing import assert_allclose

from sunpy.image.util import to_norm, un_norm


def test_to_norm():
    array_simple = np.array([10., 20., 30., 100.])
    assert_allclose(to_norm(array_simple), np.array([0.1, 0.2, 0.3, 1.]))
    array_simple_neg = np.array([-10., 0., 10., 90.])
    assert_allclose(to_norm(array_simple_neg), np.array([0, 0.1, 0.2, 1.]))


def test_un_norm():
    array_simple = np.array([10, 20, 30, 100.])
    assert_allclose(un_norm(np.array([0.1, 0.2, 0.3, 1.]), array_simple), array_simple)
    array_simple_neg = np.array([-10, 0, 10, 90])
    assert_allclose(un_norm(np.array([0, 0.1, 0.2, 1.]), array_simple_neg), array_simple_neg)
