# -*- coding: utf-8 -*-
import sunpy

def test_sysinfo():

    output = sunpy.until.sys_prop_print()    
    
    assert isinstance(output, dict)