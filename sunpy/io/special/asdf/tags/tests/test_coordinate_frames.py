import platform
from distutils.version import LooseVersion

import numpy as np
import pytest

import astropy.units as u
from asdf.tests.helpers import assert_roundtrip_tree

from sunpy.tests.helpers import asdf_entry_points
from sunpy.coordinates.frames import (Heliocentric, HeliographicCarrington,
                                      HeliographicStonyhurst, Helioprojective)

asdf = pytest.importorskip('asdf', '2.0.2')


@pytest.fixture(
    params=(Heliocentric, HeliographicCarrington, HeliographicStonyhurst,
            Helioprojective))
@asdf_entry_points
def coordframe_scalar(request):
    frame = request.param

    narg = 2
    if frame is Heliocentric:
        narg = 3

    data = np.random.random((narg, 1)) * u.arcsec
    return frame(*data, obstime='2018-01-01T00:00:00')


@pytest.fixture(
    params=(Heliocentric, HeliographicCarrington, HeliographicStonyhurst,
            Helioprojective))
@asdf_entry_points
def coordframe_array(request):
    frame = request.param

    if frame is Heliocentric:
        data = np.random.random((3, 10)) * u.km
    else:
        data = np.random.random((2, 10)) * u.arcsec

    return frame(*data, obstime='2018-01-01T00:00:00')


# Skip these two tests on windows due to a weird interaction with atomicfile
# and tmpdir
skip_windows_asdf = pytest.mark.skipif(
    (LooseVersion(asdf.__version__) < LooseVersion("2.3.1")
     and platform.system() == 'Windows'),
    reason="See https://github.com/spacetelescope/asdf/pull/632")


@skip_windows_asdf
@asdf_entry_points
def test_saveframe(coordframe_scalar, tmpdir):
    tree = {'frame': coordframe_scalar}
    assert_roundtrip_tree(tree, tmpdir)


@skip_windows_asdf
@asdf_entry_points
def test_saveframe_arr(coordframe_array, tmpdir):
    tree = {'frame': coordframe_array}
    assert_roundtrip_tree(tree, tmpdir)
