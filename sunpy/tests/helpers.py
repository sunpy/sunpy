import os
import platform
import warnings
from pathlib import Path
from functools import wraps

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy.distutils.system_info as sysinfo
import pkg_resources
import pytest

import astropy
from astropy.wcs.wcs import FITSFixedWarning

import sunpy.map

__all__ = ['skip_windows', 'skip_glymur', 'skip_ana', 'skip_32bit',
           'warnings_as_errors', 'asdf_entry_points']

# SunPy's JPEG2000 capabilities rely on the glymur library.
# First we check to make sure that glymur imports correctly before proceeding.
try:
    import glymur
except ImportError:
    SKIP_GLYMUR = True
else:
    # See if we have a C backend
    if glymur.lib.openjp2.OPENJP2:
        SKIP_GLYMUR = False
    else:
        SKIP_GLYMUR = True

try:
    from sunpy.io import _pyana  # NOQA
except ImportError:
    SKIP_ANA = True
else:
    SKIP_ANA = False

if sysinfo.platform_bits == 64:
    SKIP_32 = False
else:
    SKIP_32 = True

skip_windows = pytest.mark.skipif(platform.system() == 'Windows', reason="Windows.")
skip_glymur = pytest.mark.skipif(SKIP_GLYMUR, reason="Glymur can not be imported.")
skip_ana = pytest.mark.skipif(SKIP_ANA, reason="ANA is not available.")
skip_32bit = pytest.mark.skipif(SKIP_32, reason="Fails on a 32 bit system.")


# Skip if the SunPy ASDF entry points are missing.
asdf_entry_points = pytest.mark.skipif(not list(pkg_resources.iter_entry_points('asdf_extensions', 'sunpy')),
                                       reason="No SunPy ASDF entry points.")


@pytest.fixture
def warnings_as_errors(request):
    warnings.simplefilter('error')

    request.addfinalizer(lambda *args: warnings.resetwarnings())


new_hash_library = {}


def get_hash_library_name():
    """
    Generate the hash library name for this env.
    """
    ft2_version = f"{mpl.ft2font.__freetype_version__.replace('.', '')}"
    mpl_version = "dev" if "+" in mpl.__version__ else mpl.__version__.replace('.', '')
    astropy_version = "dev" if "dev" in astropy.__version__ else astropy.__version__.replace('.', '')
    return f"figure_hashes_mpl_{mpl_version}_ft_{ft2_version}_astropy_{astropy_version}.json"


def figure_test(test_function):
    """
    A decorator for a test that verifies the hash of the current figure or the
    returned figure, with the name of the test function as the hash identifier
    in the library. A PNG is also created in the 'result_image' directory,
    which is created on the current path.

    All such decorated tests are marked with `pytest.mark.figure` for convenient filtering.

    Examples
    --------
    @figure_test
    def test_simple_plot():
        plt.plot([0,1])
    """
    hash_library_name = get_hash_library_name()
    hash_library_file = Path(__file__).parent / hash_library_name

    @pytest.mark.remote_data
    @pytest.mark.mpl_image_compare(hash_library=hash_library_file,
                                   savefig_kwargs={'metadata': {'Software': None}},
                                   style='default')
    @wraps(test_function)
    def test_wrapper(*args, **kwargs):
        ret = test_function(*args, **kwargs)
        if ret is None:
            ret = plt.gcf()
        return ret

    return test_wrapper


# Skip coverage on this because we test it every time the CI runs --coverage!
def _patch_coverage(testdir, sourcedir):  # pragma: no cover
    """
    This function is used by the ``setup.py test`` command to change the
    filepath of the source code from the temporary directory "setup.py"
    installs the code into to the actual directory "setup.py" was executed in.
    """
    import coverage

    coveragerc = os.path.join(os.path.dirname(__file__), "coveragerc")

    # Load the .coverage file output by pytest-cov
    covfile = os.path.join(testdir, ".coverage")
    cov = coverage.Coverage(covfile, config_file=coveragerc)
    cov.load()
    cov.get_data()

    # Change the filename for the datafile to the new directory
    if hasattr(cov, "_data_files"):
        dfs = cov._data_files
    else:
        dfs = cov.data_files

    dfs.filename = os.path.join(sourcedir, ".coverage")

    # Replace the testdir with source dir
    # Lovingly borrowed from astropy (see licences directory)
    lines = cov.data._lines
    for key in list(lines.keys()):
        new_path = os.path.relpath(
            os.path.realpath(key),
            os.path.realpath(testdir))
        new_path = os.path.abspath(
            os.path.join(sourcedir, new_path))
        lines[new_path] = lines.pop(key)

    cov.save()


def no_vso(f):
    """
    Disable the VSO client from returning results via Fido during this test.
    """
    from sunpy.net import Fido
    from sunpy.net.vso import VSOClient

    @wraps(f)
    def wrapper(*args, **kwargs):
        Fido.registry[VSOClient] = lambda *args: False
        res = f(*args, **kwargs)
        Fido.registry[VSOClient] = VSOClient._can_handle_query
        return res

    return wrapper


def fix_map_wcs(smap):
    # Helper function to fix a WCS and silence the warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=FITSFixedWarning)
        wcs = smap.wcs
        wcs.fix()
    return sunpy.map.Map(smap.data, wcs)
