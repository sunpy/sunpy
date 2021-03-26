import sunpy.map
from sunpy.net import Fido, attrs as a

from sunpy.net.soar import Identifier


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
