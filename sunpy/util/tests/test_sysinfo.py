# -*- coding: utf-8 -*-

def test_sysinfo():
    import sunpy
    output = sunpy.until.sys_prop_print()    
    
    assert isinstance(output, dict)