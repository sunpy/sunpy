# -*- coding: utf-8 -*-

import pytest
asdf = pytest.importorskip('asdf')

import sunpy.map
from sunpy.data.test import get_test_filepath
from sunpy.io.special.asdf.extension import SunpyExtension

from asdf.tests.helpers import assert_roundtrip_tree


aia_path = get_test_filepath("aia_171_level1.fits")


def test_genericmap_basic(tmpdir):

    tree = {'smap': sunpy.map.Map(aia_path)}

    assert_roundtrip_tree(tree, tmpdir, extensions=SunpyExtension())
