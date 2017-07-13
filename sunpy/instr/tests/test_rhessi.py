"""
Because rhessi.py is code as well.
"""
import os

import sunpy.map
import sunpy.data.test
import sunpy.instr.rhessi as rhessi


testpath = sunpy.data.test.rootdir


def test_backprojection():
    amap = rhessi.backprojection(os.path.join(testpath, 'hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits'))
    assert isinstance(amap, sunpy.map.GenericMap)
