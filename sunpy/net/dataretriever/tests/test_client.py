from sunpy.net.dataretriever.client import QueryResponse
from sunpy.time import parse_time


def test_reprs():
    rowdict = {}
    rowdict['Start Time'] = parse_time("2012/1/1")
    rowdict['End Time'] = parse_time("2012/1/2")
    resp = QueryResponse([rowdict])
    assert isinstance(resp, QueryResponse)
    strs = ["2012-01-01T00:00:00.000", "2012-01-02T00:00:00.000"]
    assert all(s in str(resp) for s in strs)
    assert all(s in repr(resp) for s in strs)
