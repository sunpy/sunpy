from sunpy.time import parse_time
from sunpy.net.dataretriever.client import QueryResponse


def test_reprs():
    map_ = {}
    map_['Time_start'] = parse_time("2012/1/1")
    map_['Time_end'] = parse_time("2012/1/2")
    resp = QueryResponse.create(map_, [''])
    assert isinstance(resp, QueryResponse)
    strs = ["2012-01-01 00:00:00", "2012-01-02 00:00:00"]
    assert all(s in str(resp) for s in strs)
    assert all(s in repr(resp) for s in strs)
