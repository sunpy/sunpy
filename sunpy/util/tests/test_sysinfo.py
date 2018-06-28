from __future__ import absolute_import

import sunpy

# -*- coding: utf-8 -*-


def test_sysinfo():

    output = sunpy.util.get_sys_dict()

    assert isinstance(output, dict)
