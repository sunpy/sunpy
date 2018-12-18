import numpy as np
import pytest

import astropy.units as u
pytest.importorskip('asdf', '2.0.2')
from asdf.tests.helpers import assert_roundtrip_tree

from sunpy.tests.helpers import skip_windows
from sunpy.coordinates.frames import (Heliocentric, Helioprojective,
                                      HeliographicCarrington, HeliographicStonyhurst)


@pytest.fixture(params=(Heliocentric, HeliographicCarrington,
                        HeliographicStonyhurst, Helioprojective))
def coordframe_scalar(request):
    frame = request.param

    narg = 2
    if frame is Heliocentric:
        narg = 3

    data = np.random.random((narg, 1)) * u.arcsec
    return frame(*data, obstime='2018-01-01T00:00:00')


@pytest.fixture(params=(Heliocentric, HeliographicCarrington, HeliographicStonyhurst,
                        Helioprojective))
def coordframe_array(request):
    frame = request.param

    if frame is Heliocentric:
        data = np.random.random((3, 10)) * u.km
    else:
        data = np.random.random((2, 10)) * u.arcsec

    return frame(*data, obstime='2018-01-01T00:00:00')


# Skip these two tests on windows due to a weird interaction with atomicfile
# and tmpdir


@skip_windows
def test_saveframe(coordframe_scalar, tmpdir):
    tree = {'frame': coordframe_scalar}
    assert_roundtrip_tree(tree, tmpdir)


@skip_windows
def test_saveframe_arr(coordframe_array, tmpdir):
    tree = {'frame': coordframe_array}
    assert_roundtrip_tree(tree, tmpdir)
