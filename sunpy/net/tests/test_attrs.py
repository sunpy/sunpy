from sunpy.net import attr
from sunpy.net import attrs as a
from sunpy.time import parse_time


def test_time_subclass():
    class NewTime(a.Time):
        pass

    assert isinstance(NewTime("2020/01/01", "2020/01/02") &
                      a.Time("2020/02/02", "2020/02/03"), attr.AttrAnd)


def test_attrs_time():
    times = a.Time("2020/10/01T00:00", "2020/10/01T00:00")
    times.start == parse_time("2020/10/01T00:00")
    times.end == parse_time("2020/10/01T00:00")
