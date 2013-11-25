from __future__ import absolute_import

from urllib2 import urlopen
from sunpy.time import TimeRange
from sunpy.instr import lyra
import tempfile
import os
import pytest

@pytest.mark.online
def test_lytaf_utils():
    '''test the downloading of the LYTAF file'''
    tmp_dir=tempfile.mkdtemp() + '/'
    lyra.download_lytaf_database(lytaf_dir=tmp_dir)
    assert os.path.exists(tmp_dir + 'annotation_ppt.db') == True

    #try doing a query on the temporary database
    lar=lyra.get_lytaf_events(TimeRange('2010-06-13 02:00','2010-06-13 06:00'),lytaf_dir=tmp_dir)
    assert type(lar) == list
    assert type(lar[0]) == dict
    assert type(lar[0]['start_time']) == datetime.datetime
    assert type(lar[0]['end_time']) == datetime.datetime
    assert type(lar[0]['roi_description'] == str
    assert type(lar[0]['event_type_description']) == str
    
    
    
