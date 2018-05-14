import pytest
import numpy as np

# These only work with astropy 3.1 or newer, so skip them on everything but the
# master version of 3.1 and up.
pytest.importorskip('astropy', '3.1.0.dev0')

import astropy.units as u

from sunpy.coordinates.frames import (Heliocentric, HeliographicCarrington,
                                      HeliographicStonyhurst, Helioprojective)

from asdf.tests.helpers import assert_roundtrip_tree


@pytest.fixture(params=(Heliocentric, HeliographicCarrington,
                        HeliographicStonyhurst, Helioprojective))
def coordframe(request):
    frame = request.param

    narg = 2
    if frame is Heliocentric:
        narg = 3

    data = np.random.random((narg, 10)) * u.arcsec
    return frame(*data, obstime='2018-01-01T00:00:00')


def test_saveframe(coordframe, tmpdir):
    tree = {'frame': coordframe}
    assert_roundtrip_tree(tree, tmpdir)
