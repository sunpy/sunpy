from __future__ import absolute_import

from datetime import datetime

from sunpy.util import util

LANDING = datetime(1966, 2, 3)

def test_anytim_tuple():
    assert util.anytim((1966, 2, 3)) == LANDING

def test_anytim_ISO():
    assert util.anytim('1966-02-03') == LANDING
    assert (
        util.anytim('1966-02-03T20:17:40') == datetime(1966, 2, 3, 20, 17, 40)
    )
    assert (
        util.anytim('19660203T201740') == datetime(1966, 2, 3, 20, 17, 40)
    )
    
    lst = [
        ('2007-05-04T21:08:12.999999', datetime(2007, 5, 4, 21, 8, 12, 999999)),
        ('20070504T210812.999999', datetime(2007, 5, 4, 21, 8, 12, 999999)),
        ('2007/05/04 21:08:12.999999', datetime(2007, 5, 4, 21, 8, 12, 999999)),
        ('2007-05-04 21:08:12.999999', datetime(2007, 5, 4, 21, 8, 12, 999999)),
        ('2007/05/04 21:08:12', datetime(2007, 5, 4, 21, 8, 12)),
        ('2007-05-04 21:08:12', datetime(2007, 5, 4, 21, 8, 12)),
        ('2007-05-04T21:08:12', datetime(2007, 5, 4, 21, 8, 12)),
        ('20070504T210812', datetime(2007, 5, 4, 21, 8, 12)),
        ('2007-May-04 21:08:12', datetime(2007, 5, 4, 21, 8, 12)),
        ('2007-May-04', datetime(2007, 5, 4)),
        ('2007-05-04', datetime(2007, 5, 4)),
        ('2007/05/04', datetime(2007, 5, 4)),
    ]
    
    for k, v in lst:
        assert util.anytim(k) == v