# -*- coding: utf-8 -*-

import pytest
asdf = pytest.importorskip('asdf', '2.0')

from asdf.tests.helpers import assert_roundtrip_tree

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.io.special.asdf.extension import SunpyExtension


aia_path = get_test_filepath("aia_171_level1.fits")


@pytest.fixture
def aia171_test_map():
    return sunpy.map.Map(aia_path)


def test_genericmap_basic(aia171_test_map, tmpdir):

    tree = {'smap': aia171_test_map}

    assert_roundtrip_tree(tree, tmpdir, extensions=SunpyExtension())
