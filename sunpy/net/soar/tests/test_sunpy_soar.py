from sunpy.net import Fido, attrs as a

from sunpy.net.soar import Identifier


def test_search():
    id = a.Instrument('EUI')
    time = a.Time('2020-02-01', '2021-02-20')
    level = a.Level(2)
    identifier = Identifier('EUI-FSI174-IMAGE')

    res = Fido.search(id, time, level, identifier)
    assert len(res) == 1
