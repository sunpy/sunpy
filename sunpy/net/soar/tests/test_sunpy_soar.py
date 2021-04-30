import sunpy.map
from sunpy.net import Fido, attrs as a

from sunpy.net.soar import Identifier
from sunpy.net.soar.client import SOARClient


def test_search():
    id = a.Instrument('EUI')
    time = a.Time('2021-02-01', '2021-02-02')
    level = a.Level(1)
    identifier = Identifier('EUI-FSI174-IMAGE')

    res = Fido.search(id, time, level, identifier)
    assert len(res) == 1
    assert len(res[0]) == 43

    files = Fido.fetch(res[0, 0])
    assert len(files) == 1

    eui_map = sunpy.map.Map(files[0])


def test_insitu_search():
    id = a.Instrument('MAG')
    time = a.Time('2020-04-16', '2020-04-17')
    identifier = Identifier('MAG-RTN-NORMAL-1-MINUTE')

    res = Fido.search(id, time, identifier)
    assert len(res) == 1
    assert len(res[0]) == 2

    files = Fido.fetch(res[0, 0])
    assert len(files) == 1


def test_no_results():
    id = a.Instrument('EUI')
    time = a.Time('2019-02-01', '2019-02-02')
    query = id & time

    res = SOARClient().search(query)
    assert len(res) == 0
