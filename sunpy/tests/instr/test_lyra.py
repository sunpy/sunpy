from __future__ import absolute_import

from urllib2 import urlopen
from sunpy.time import TimeRange
from sunpy.instr import lyra

def test_lytaf_url_exists():
    url='http://proba2.oma.be/lyra/data/lytaf/annotation_ppt.db'
    assert urlopen(url).code == 200

def test_lytaf_query():
    """Have a small sample of the LYTAF database available in data/tests for testing"""
    dir='../../data/test/'
    lar=lyra.get_lytaf_events(TimeRange('2010-06-13 02:00','2010-06-13 06:00'),lytaf_dir=dir)
    assert type(lar) == list
    assert type(lar[0]) == dict
    
