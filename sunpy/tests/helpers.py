# -*- coding: utf-8 -*-
# Lovingly borrowed from Astropy
# Licensed under a 3-clause BSD style license - see licences/ASTROPY.rst

import warnings

import pytest

import astropy.units as u
from astropy.utils.decorators import wraps
import numpy as np
import matplotlib.pyplot as plt

from sunpy.tests import hash

__all__ = ['skip_windows', 'skip_glymur', 'skip_ana', 'skip_wcsaxes', 'warnings_as_errors']

# SunPy's JPEG2000 capabilities rely on the glymur library.  First we check to
# make sure that glymur imports correctly before proceeding.
try:
    import glymur
except ImportError:
    SKIP_GLYMUR = True
else:
    # See if we have a C backend
    if any((glymur.lib.openjp2.OPENJP2, glymur.lib.openjpeg.OPENJPEG)):
        SKIP_GLYMUR = False
    else:
        SKIP_GLYMUR = True

# Skip ana tests if we are on Windows or we can't import the c extension.
import platform
if platform.system() == 'Windows':
    SKIP_ANA = True
else:
    SKIP_ANA = False

try:
    import sunpy.io._pyana
except ImportError:
    SKIP_ANA = True
else:
    SKIP_ANA = SKIP_ANA or False


try:
    import wcsaxes
except ImportError:
    SKIP_WCSAXES = True
else:
    SKIP_WCSAXES = False

skip_windows = pytest.mark.skipif(platform.system() == 'Windows', reason="Windows")

skip_glymur = pytest.mark.skipif(SKIP_GLYMUR, reason="Glymur can not be imported")

skip_ana = pytest.mark.skipif(SKIP_ANA, reason="ANA is not available")

skip_wcsaxes = pytest.mark.skipif(SKIP_WCSAXES, reason="wcsaxes is not available")

@pytest.fixture
def warnings_as_errors(request):
    warnings.simplefilter('error')

    request.addfinalizer(lambda *args: warnings.resetwarnings())


def assert_quantity_allclose(actual, desired, rtol=1.e-7, atol=0, err_msg='', verbose=True):
    """
    Raise an assertion if two objects are not equal up to desired tolerance.

    This is a :class:`~astropy.units.Quantity`-aware version of
    :func:`numpy.testing.assert_allclose`.
    """

    if isinstance(actual, u.Quantity) and isinstance(desired, u.Quantity):

        if atol != 0:
            if not isinstance(atol, u.Quantity):
                raise TypeError("If `actual` and `desired` are Quantities, `atol` parameter should also be a Quantity")
            else:
                atol = atol.to(actual.unit).value

        np.testing.assert_allclose(actual.value, desired.to(actual.unit).value,
                                   rtol=rtol, atol=atol, err_msg=err_msg, verbose=verbose)

    elif isinstance(actual, u.Quantity):
        raise TypeError("If `actual` is a Quantity, `desired` should also be a Quantity")

    elif isinstance(desired, u.Quantity):
        raise TypeError("If `desired` is a Quantity, `actual` should also be a Quantity")

    else:

        if isinstance(atol, u.Quantity):
            raise TypeError("If `actual` and `desired` are not Quantities, `atol` parameter should also not be a Quantity")

        np.testing.assert_allclose(actual, desired,
                                   rtol=rtol, atol=atol, err_msg=err_msg, verbose=verbose)


def figure_test(test_function):
    """
    A decorator for a test that verifies the hash of the current figure or the returned figure,
    with the name of the test function as the hash identifier in the library.

    All such decorated tests are marked with `pytest.mark.figure` for convenient filtering.

    Examples
    --------
    @figure_test
    def test_simple_plot():
        plt.plot([0,1])
    """
    @pytest.mark.figure
    @wraps(test_function)
    def wrapper(*args, **kwargs):
        plt.figure()
        name = "{0}.{1}".format(test_function.__module__, test_function.__name__)
        figure_hash = hash.hash_figure(test_function(*args, **kwargs))
        if name not in hash.hash_library:
            hash.hash_library[name] = figure_hash
            pytest.fail("Hash not present: {0}".format(name))
        else:
            assert hash.hash_library[name] == figure_hash
    return wrapper
