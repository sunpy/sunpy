"""
Lyra Tests
"""
from __future__ import absolute_import

import pytest

#pylint: disable=C0103,R0904,W0201,W0232,E1103
import sunpy
from sunpy.time import parse_time

@pytest.mark.online
@pytest.mark.parametrize("date,level,start,end",
[('2012/06/03',3,'2012-06-03 00:00:00.047000','2012-06-03 23:59:00.047000'),
 ('2012/06/03',2,'2012-06-03 00:00:00.047000','2012-06-03 23:59:59.047000')
])
def test_lyra_level(date, level,start, end):
    lyra = sunpy.lightcurve.LYRALightCurve.create(date, level=level)
    assert isinstance(lyra, sunpy.lightcurve.LYRALightCurve)
    assert lyra.time_range().start() == parse_time(start)
    assert lyra.time_range().end() == parse_time(end)



@pytest.mark.online
@pytest.mark.parametrize("date,start,end",
[('2011/02/06', '2011-02-06 00:00:00.008500', '2011-02-06 23:59:59.008500'),
 ('2012/08/04', '2012-08-04 00:00:00.043000', '2012-08-04 23:59:59.043000'),
 ])

def test_lyra_date(date, start, end):
   
    lyra = sunpy.lightcurve.LYRALightCurve.create(date)
    assert isinstance(lyra, sunpy.lightcurve.LYRALightCurve)
    assert lyra.time_range().start() == parse_time(start)
    assert lyra.time_range().end() == parse_time(end)

@pytest.mark.online
@pytest.mark.parametrize(("url, start, end"),
[('http://proba2.oma.be/lyra/data/bsd/2012/02/04/lyra_20120204-000000_lev3_std.fits', '2012-02-04 00:00:00.004000', '2012-02-04 23:59:00.004000'),
('http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev3_std.fits', '2011-08-10 00:00:00.020000', '2011-08-10 23:59:00.020000')])

def test_online(url, start, end):
   
    lyra = sunpy.lightcurve.LYRALightCurve.create(url)
    assert isinstance(lyra, sunpy.lightcurve.LYRALightCurve)
    assert lyra.time_range().start() == parse_time(start)
    assert lyra.time_range().end() == parse_time(end)

 

