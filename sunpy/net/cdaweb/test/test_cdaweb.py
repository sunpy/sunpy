import pytest

import sunpy.net.attrs as a
from sunpy.net import Fido
from sunpy.net.cdaweb.attrs import Dataset


@pytest.mark.remote_data
def test_query():
    res = Fido.search(a.Time('2018-11-01', '2018-11-01 01:00:00'),
                      Dataset('WI_H1_SWE') | Dataset('WI_H5_SWE'))
    assert len(res) == 2
    assert len(res[0]) == 1
    assert len(res[1]) == 2

    files = Fido.fetch(res)
    assert len(files) == 3
