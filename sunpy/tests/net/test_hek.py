# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

#pylint: disable=W0613

from __future__ import absolute_import

import pytest

from sunpy.net import hek


def test_eventtype_collide():
    with pytest.raises(TypeError):
        hek.attrs.AR & hek.attrs.CE
    with pytest.raises(TypeError):
        (hek.attrs.AR & hek.attrs.Time((2011, 1, 1), (2011, 1, 2))) & hek.attrs.CE
        with pytest.raises(TypeError):
            (hek.attrs.AR | hek.attrs.Time((2011, 1, 1), (2011, 1, 2))) & hek.attrs.CE


def test_eventtype_or():
    assert (hek.attrs.AR | hek.attrs.CE).item == "ar,ce"
