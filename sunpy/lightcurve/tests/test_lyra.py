"""
Lyra Tests
"""
from __future__ import absolute_import

import pytest

#pylint: disable=C0103,R0904,W0201,W0232,E1103
import sunpy
import matplotlib as mpl
from matplotlib.testing.decorators import cleanup
from sunpy.time import parse_time

def test_lyra_level():
    lyra1 = sunpy.lightcurve.LYRALightCurve.create('2012/06/03',level=3)
    lyra2 = sunpy.lightcurve.LYRALightCurve.create('2012/06/03',level=2)
    assert isinstance(lyra1,sunpy.lightcurve.LYRALightCurve)
    assert isinstance(lyra2,sunpy.lightcurve.LYRALightCurve)

    assert lyra1.time_range().start() == parse_time('2012-06-03 00:00:00.047000')
    assert lyra1.time_range().end() == parse_time('2012-06-03 00:23:59.047000')


@pytest.mark.parametrize(('date','start','end'),
[('2011/02/06','2011-02-06 00:00:00.008500','2011-02-06 23:59:59.008500'),
 ('2012/08/04','2012-08-04 00:00:00.043000','2012-08-04 23:59:59.043000'),
 ])
def test_lyra_date(date,start,end):
   
   lyra=sunpy.lightcurve.LYRALightCurve.create(date)
   assert isinstance(lyra,sunpy.lightcurve.LYRALightCurve)
   assert lyra.time_range().start()==parse_time(start)
   assert lyra.time_range().end() == parse_time(end)

@pytest.mark.parametrize(('url','start','end'),
[('http://proba2.oma.be/lyra/data/bsd/2012/02/04/lyra_20120204-000000_lev3_std.fits','2012-02-04 00:00:00.004000','2012-02-04 00:23:59.004000'),
('http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev3_std.fits','2011-08-10 00:00:00.020000','2011-08-10 00:23:59.020000')])

def test_online(url,start,end):
   
   lyra=sunpy.lightcurve.LYRALightCurve.create(url)
   assert isinstance(lyra,sunpy.lightcurve.LYRALightCurve)
   assert lyra.time_range().start()==parse_time(start)
   assert lyra.time_range().end() == parse_time(end)

 

#@cleanup    
#def test_peek():
#    lyra = sunpy.lightcurve.LYRALightCurve.create(
#    "http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev2_std.fits")
#    assert isinstance(lyra.peek(),mpl.figure.Figure)
