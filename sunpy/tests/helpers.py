import sys
import platform
import warnings
from pathlib import Path
from platform import python_version
from functools import wraps
from importlib.metadata import entry_points

import matplotlib as mpl
import matplotlib.pyplot as plt
import pytest
from packaging.version import Version

import astropy
from astropy.wcs.wcs import FITSFixedWarning

import sunpy.map

__all__ = ['skip_windows', 'skip_glymur', 'skip_ana', 'warnings_as_errors', 'asdf_entry_points']

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

if sys.maxsize > 2**32:
    SKIP_32 = False
else:
    SKIP_32 = True

skip_windows = pytest.mark.skipif(platform.system() == "Windows", reason="Windows.")
skip_glymur = pytest.mark.skipif(SKIP_GLYMUR, reason="Glymur can not be imported.")
skip_ana = pytest.mark.skipif(SKIP_ANA, reason="ANA is not available.")
if Version(python_version()) >= Version("3.10.0"):
    asdf_entry_points = pytest.mark.skipif(
        not entry_points().select(group="asdf.resource_mappings", name="sunpy"),
        reason="No SunPy ASDF entry points.",
    )
else:
    asdf_entry_points = pytest.mark.skipif(
        not any(
            [
                enter_point.name == "sunpy"
                for enter_point in entry_points()["asdf.resource_mappings"]
            ]
        ),
        reason="No SunPy ASDF entry points.",
    )


@pytest.fixture
def warnings_as_errors():
    warnings.simplefilter('error')
    yield
    warnings.resetwarnings()


def get_hash_library_name():
    """
    Generate the hash library name for this env.
    """
    import mpl_animators
    version = mpl_animators.__version__
    animators_version = "dev" if "+" in version else version.replace('.', '')
    ft2_version = f"{mpl.ft2font.__freetype_version__.replace('.', '')}"
    mpl_version = "dev" if ("+" in mpl.__version__) or ("rc" in mpl.__version__) else mpl.__version__.replace('.', '')
    astropy_version = "dev" if "dev" in astropy.__version__ else astropy.__version__.replace('.', '')
    return f"figure_hashes_mpl_{mpl_version}_ft_{ft2_version}_astropy_{astropy_version}_animators_{animators_version}.json"


def figure_test(test_function):
    """
    A decorator for a test that verifies the hash of the current figure or the
    returned figure, with the name of the test function as the hash identifier
    in the library. A PNG is also created in the 'result_image' directory,
    which is created on the current path.

    All such decorated tests are marked with `pytest.mark.mpl_image` for convenient filtering.

    Examples
    --------
    @figure_test
    def test_simple_plot():
        plt.plot([0,1])
    """
    hash_library_name = get_hash_library_name()
    hash_library_file = Path(__file__).parent / hash_library_name

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
