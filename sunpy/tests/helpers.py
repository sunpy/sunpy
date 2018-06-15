# -*- coding: utf-8 -*-
# Lovingly borrowed from Astropy
# Licensed under a 3-clause BSD style license - see licences/ASTROPY.rst

import os
import datetime
import platform
import warnings

import pytest
import matplotlib.pyplot as plt

from astropy.utils.decorators import wraps

from sunpy.tests import hash

__all__ = ['skip_windows', 'skip_glymur', 'skip_ana', 'warnings_as_errors']

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

try:
    from sunpy.io import _pyana
except ImportError:
    SKIP_ANA = True
else:
    SKIP_ANA = False

skip_windows = pytest.mark.skipif(platform.system() == 'Windows', reason="Windows")
skip_glymur = pytest.mark.skipif(SKIP_GLYMUR, reason="Glymur can not be imported")
skip_ana = pytest.mark.skipif(SKIP_ANA, reason="ANA is not available")


@pytest.fixture
def warnings_as_errors(request):
    warnings.simplefilter('error')

    request.addfinalizer(lambda *args: warnings.resetwarnings())


new_hash_library = {}
figure_base_dir = None


def figure_test(test_function):
    """
    A decorator for a test that verifies the hash of the current figure or the returned figure,
    with the name of the test function as the hash identifier in the library.
    A PNG is also created in the 'result_image' directory, which is created
    on the current path.

    All such decorated tests are marked with `pytest.mark.figure` for
    convenient filtering.

    Examples
    --------
    @figure_test
    def test_simple_plot():
        plt.plot([0,1])
    """
    @pytest.mark.figure
    @wraps(test_function)
    def wrapper(*args, **kwargs):
        if not os.path.exists(hash.HASH_LIBRARY_FILE):
            pytest.xfail('Could not find a figure hash library at {}'.format(hash.HASH_LIBRARY_FILE))
        if figure_base_dir is None:
            pytest.xfail("No directory to save figures to found")

        name = "{0}.{1}".format(test_function.__module__,
                                test_function.__name__)
        # Run the test function and get the figure
        plt.figure()
        fig = test_function(*args, **kwargs)
        if fig is None:
            fig = plt.gcf()

        # Save the image that was generated
        if not os.path.exists(figure_base_dir):
            os.makedirs(figure_base_dir)
        result_image_loc = os.path.join(figure_base_dir, '{}.png'.format(name))
        plt.savefig(result_image_loc)
        plt.close()

        # Create hash
        imgdata = open(result_image_loc, "rb")
        figure_hash = hash._hash_file(imgdata)
        imgdata.close()

        new_hash_library[name] = figure_hash
        if name not in hash.hash_library:
            pytest.fail("Hash not present: {0}".format(name))

        if hash.hash_library[name] != figure_hash:
            raise RuntimeError('Figure hash does not match expected hash.\n'
                               'New image generated and placed at {}'.format(result_image_loc))

    return wrapper
