"""
isort:skip_file.
"""
import platform
from distutils.version import LooseVersion

import numpy as np
import pytest

import astropy.units as u
asdf = pytest.importorskip('asdf', '2.0')
from asdf.tests.helpers import assert_roundtrip_tree  # noqa

import sunpy.map  # noqa
from sunpy.data.test import get_test_filepath  # noqa
from sunpy.tests.helpers import asdf_entry_points
from sunpy.io.special.asdf.extension import SunpyExtension  # noqa


@pytest.fixture
def aia171_test_map():
    aia_path = get_test_filepath("aia_171_level1.fits")
    return sunpy.map.Map(aia_path)


# Skip these two tests on windows due to a weird interaction with atomicfile
# and tmpdir
skip_windows_asdf = pytest.mark.skipif(
    (LooseVersion(asdf.__version__) < LooseVersion("2.3.1")
     and platform.system() == 'Windows'),
    reason="See https://github.com/spacetelescope/asdf/pull/632")


@skip_windows_asdf
@asdf_entry_points
def test_genericmap_basic(aia171_test_map, tmpdir):

    tree = {'smap': aia171_test_map}

    assert_roundtrip_tree(tree, tmpdir, extensions=SunpyExtension())


@skip_windows_asdf
@asdf_entry_points
def test_genericmap_mask(aia171_test_map, tmpdir):

    mask = np.zeros_like(aia171_test_map.data)
    mask[10, 10] = 1

    aia171_test_map.mask = mask
    aia171_test_map._unit = u.m

    tree = {'smap': aia171_test_map}

    assert_roundtrip_tree(tree, tmpdir, extensions=SunpyExtension())
