# -*- coding: utf-8 -*-
# Lovingly borrowed from Astropy
# Licensed under a 3-clause BSD style license - see licences/ASTROPY.rst

from __future__ import absolute_import, division, print_function
import warnings
import tempfile
import platform
import os
import urllib.request

import pytest
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
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
except ImportError as e:
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
figure_test_pngfiles = {}
remote_images = 'https://raw.githubusercontent.com/dstansby/sunpy-test-images/master/images/'


def figure_test(test_function):
    """
    A decorator for a test that verifies the hash of the current figure or the returned figure,
    with the name of the test function as the hash identifier in the library.
    A PNG is also created with a temporary filename, with the lookup stored in the
    `figure_test_pngfiles` dictionary.

    All such decorated tests are marked with `pytest.mark.figure` for convenient filtering.

    Examples
    --------
    @figure_test
    def test_simple_plot():
        plt.plot([0,1])
    """
    @pytest.mark.remote_data
    @pytest.mark.figure
    @wraps(test_function)
    def wrapper(*args, **kwargs):
        if not os.path.exists(hash.HASH_LIBRARY_FILE):
            pytest.xfail('Could not find a figure hash library at {}'.format(hash.HASH_LIBRARY_FILE))

        name = "{0}.{1}".format(test_function.__module__,
                                test_function.__name__)
        plt.figure()
        pngfile = tempfile.NamedTemporaryFile(delete=False)
        figure_hash = hash.hash_figure(test_function(*args, **kwargs), out_stream=pngfile)
        figure_test_pngfiles[name] = pngfile.name
        pngfile.close()

        new_hash_library[name] = figure_hash
        if name not in hash.hash_library:
            pytest.fail("Hash not present: {0}".format(name))

        if hash.hash_library[name] != figure_hash:
            # Save the image that was actually generated
            os.makedirs('result_images', exist_ok=True)
            result_image_loc = 'result_images/{}.png'.format(name)
            plt.savefig(result_image_loc)
            # Download the expected image from the sunpy-test-images repo
            url = '{}/{}.png'.format(remote_images, name)
            expected_image_loc = 'result_images/{}_expected.png'.format(name)
            urllib.request.urlretrieve(url, expected_image_loc)
            raise RuntimeError('Figure hash does not match expected hash.\n'
                               'New image generated and placed at {}\n'
                               'Expected image downloaded to {}'.format(result_image_loc, expected_image_loc))

        plt.close()
    return wrapper
