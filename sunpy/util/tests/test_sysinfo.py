# -*- coding: utf-8 -*-
import sunpy

def test_sysinfo():

    output = sunpy.util.get_sys_dict()

    assert isinstance(output, dict)
