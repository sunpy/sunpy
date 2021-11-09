import os

import astropy.units as u
import pytest

import sunpy.map
from sunpy.net import Fido, attrs as a
from sunpy.util.exceptions import SunpyDeprecationWarning

from sunpy.net.soar.client import SOARClient


def test_search():
    id = a.Instrument('EUI')
    time = a.Time('2021-02-01', '2021-02-02')
    level = a.Level(1)
    product = a.soar.Product('EUI-FSI174-IMAGE')

    res = Fido.search(id, time, level, product)
    assert len(res) == 1
    assert len(res[0]) == 43
    assert u.allclose(res[0, 0]['Filesize'], 4.740*u.Mbyte)

    files = Fido.fetch(res[0, 0])
    assert len(files) == 1
    fname = files[0]
    assert u.allclose(os.path.getsize(fname) * u.byte,
                      res[0, 0]['Filesize'], atol=1e-3*u.Mbyte)
    # Smoke test that we can read this into a map
    eui_map = sunpy.map.Map(fname)


def test_deprecated_identifier():
    id = a.Instrument('EUI')
    time = a.Time('2021-02-01', '2021-02-02')
    level = a.Level(1)
    with pytest.warns(SunpyDeprecationWarning):
        identifier = a.soar.Identifier('EUI-FSI174-IMAGE')
    product = a.soar.Product('EUI-FSI174-IMAGE')
    res1 = Fido.search(id, time, level, identifier)
    res2 = Fido.search(id, time, level, product)

    assert res1.__str__() == res2.__str__()


def test_insitu_search():
    id = a.Instrument('MAG')
    time = a.Time('2020-04-16', '2020-04-17')
    identifier = a.soar.Product('MAG-RTN-NORMAL-1-MINUTE')

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


def test_no_instrument():
    # Check that a time only search returns results
    time = a.Time('2020-04-16', '2020-04-17')
    res = SOARClient().search(time)
    assert len(res) == 50
