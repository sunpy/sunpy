import sunpy.net.solarnet as solarnet
import sunpy.net.solarnet.attrs as attrs
from sunpy.net.attr import AttrOr


def test_complex_query():
    sol = solarnet.solarnet.and_(attrs.Dataset("test_data") |attrs.Limit(23) & attrs.Dataset("test_2_data") )
    assert isinstance(sol,AttrOr)
    assert len(sol.attrs) == 2
