import os
import urllib
import platform
import warnings
from functools import wraps

import matplotlib.pyplot as plt
import numpy.distutils.system_info as sysinfo
import pkg_resources
import pytest
from matplotlib.testing import compare

from sunpy.tests import hash

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
    @pytest.mark.figure
    @wraps(test_function)
    def wrapper(*args, **kwargs):
        if not os.path.exists(hash.HASH_LIBRARY_FILE):
            pytest.xfail(f'Could not find a figure hash library at {hash.HASH_LIBRARY_FILE}')
        # figure_base_dir is a pytest fixture defined on use.
        if figure_base_dir is None:
            pytest.xfail("No directory to save figures to found")

        name = "{}.{}".format(test_function.__module__,
                              test_function.__name__)
        # Run the test function and get the figure
        plt.figure()
        fig = test_function(*args, **kwargs)
        if fig is None:
            fig = plt.gcf()

        # Save the image that was generated
        figure_base_dir.mkdir(exist_ok=True)
        result_image_loc = figure_base_dir / f'{name}.png'
        plt.savefig(str(result_image_loc))
        plt.close()

        # Create hash
        imgdata = open(result_image_loc, "rb")
        figure_hash = hash._hash_file(imgdata)
        imgdata.close()

        new_hash_library[name] = figure_hash
        if name not in hash.hash_library:
            pytest.fail(f"Hash not present: {name}")

        if hash.hash_library[name] != figure_hash:
            raise RuntimeError('Figure hash does not match expected hash.\n'
                               'New image generated and placed at {}'.format(result_image_loc))

    return wrapper


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


html_intro = '''
<!DOCTYPE html>
<html>
<head>
<style>
table, th, td {
    border: 1px solid black;
}
</style>
</head>
<body>

<h2>Image test comparison</h2>

<table>
  <tr>
    <th>Test Name</th>
    <th>Baseline image</th>
    <th>Diff</th>
    <th>New image</th>
  </tr>
'''


def _generate_fig_html(fname):
    generated_image = figure_base_dir / (fname + '.png')

    envname = os.environ.get("TOXENV", "figure_py36")
    # Download baseline image
    baseline_url = f'https://raw.githubusercontent.com/sunpy/sunpy-figure-tests/sunpy-master/figures/{envname}/'
    baseline_image_url = baseline_url + generated_image.name
    baseline_image = figure_base_dir / "reference_images" / generated_image.name
    baseline_image_exists = baseline_image.exists()
    if not baseline_image_exists:
        baseline_image.parent.mkdir(parents=True, exist_ok=True)
        try:
            urllib.request.urlretrieve(baseline_image_url, baseline_image)
            baseline_image_exists = True
        except urllib.error.HTTPError:
            pass

    # Create diff between baseline and generated image
    diff_image = figure_base_dir / "difference_images" / generated_image.name
    diff_image.parent.mkdir(parents=True, exist_ok=True)
    if baseline_image_exists:
        compare.save_diff_image(str(baseline_image), str(generated_image), str(diff_image))

    html_block = ('<tr>'
                  '<td>{}\n'.format(generated_image.stem) +
                  f'<td><img src="{baseline_image.relative_to(figure_base_dir)}"></td>\n' +
                  f'<td><img src="{diff_image.relative_to(figure_base_dir)}"></td>\n' +
                  f'<td><img src="{generated_image.relative_to(figure_base_dir)}"></td>\n' +
                  '</tr>\n\n')
    return html_block


def generate_figure_webpage(hash_library):
    html_file = figure_base_dir / 'fig_comparison.html'
    with open(html_file, 'w') as f:
        f.write(html_intro)
        for fname in hash_library:
            f.write(_generate_fig_html(fname))
        f.write('</table>')
        f.write('</body>')
        f.write('</html>')


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
