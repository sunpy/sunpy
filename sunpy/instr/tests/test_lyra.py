from __future__ import absolute_import

from sunpy.time import TimeRange, parse_time
from sunpy.instr import lyra
import tempfile
import os
import pytest
import datetime
import numpy as np

@pytest.mark.online
def test_lytaf_utils():
    '''test the downloading of the LYTAF file and subsequent queries'''
    tmp_dir=tempfile.mkdtemp()
    lyra.download_lytaf_database(lytaf_dir=tmp_dir)
    assert os.path.exists(os.path.join(tmp_dir,'annotation_ppt.db'))

    #try doing a query on the temporary database
    lar=lyra.get_lytaf_events(TimeRange('2010-06-13 02:00','2010-06-13 06:00'),lytaf_dir=tmp_dir)
    assert type(lar) == list
    assert type(lar[0]) == dict
    assert type(lar[0]['start_time']) == datetime.datetime
    assert type(lar[0]['end_time']) == datetime.datetime
    assert type(lar[0]['roi_description']) == str
    assert type(lar[0]['event_type_description']) == str
    assert lar[0]['start_time'] == parse_time('2010-06-13 02:07:04')
    assert lar[0]['end_time'] == parse_time('2010-06-13 02:10:04')
    assert lar[0]['event_type_description'] == 'LAR'

    #test split_series_using_lytaf
    #construct a dummy signal for testing purposes
    basetime=parse_time('2010-06-13 02:00')
    seconds=3600
    dummy_time = [basetime + datetime.timedelta(0, s) for s in range(seconds)]
    dummy_data=np.random.random(seconds)

    split=lyra.split_series_using_lytaf(dummy_time, dummy_data, lar)
    assert type(split) == list
    assert len(split) == 4
    assert split[0]['subtimes'][0] == datetime.datetime(2010, 6, 13, 2, 0)
    assert split[0]['subtimes'][-1] == datetime.datetime(2010, 6, 13, 2, 7, 2)
    assert split[3]['subtimes'][0] == datetime.datetime(2010, 6, 13, 2, 59, 41)
    assert split[3]['subtimes'][-1] == datetime.datetime(2010, 6, 13, 2, 59, 58)
    
