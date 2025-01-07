import sunpy.net.solarnet as solarnet
import sunpy.net.solarnet.attrs as attrs
from sunpy.net.attr import AttrOr


def test_complex_query():
    a1 = attrs.Dataset("test_data")
    a2 = attrs.Limit(23)
    sol = solarnet.solarnet.and_(a1 | a2)
    assert isinstance(sol,AttrOr)
    assert len(sol.attrs) == 2
