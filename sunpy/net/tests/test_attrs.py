from sunpy.net import attr
from sunpy.net import attrs as a


def test_time_subclass():
    class NewTime(a.Time):
        pass

    assert isinstance(NewTime("2020/01/01", "2020/01/02") &
                      a.Time("2020/02/02", "2020/02/03"), attr.AttrAnd)
