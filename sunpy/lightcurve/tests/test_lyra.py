"""
Lyra Tests
"""
from __future__ import absolute_import

import pytest

#pylint: disable=C0103,R0904,W0201,W0232,E1103
import sunpy
import matplotlib as mpl
from matplotlib.testing.decorators import cleanup

@pytest.mark.online
def test_lyra():
    lyra = sunpy.lightcurve.LYRALightCurve.create(
    "http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev2_std.fits")
    assert isinstance(lyra, sunpy.lightcurve.LYRALightCurve)

#@cleanup    
#def test_peek():
#    lyra = sunpy.lightcurve.LYRALightCurve.create(
#    "http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev2_std.fits")
#    assert isinstance(lyra.peek(),mpl.figure.Figure)
