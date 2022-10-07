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


@pytest.mark.remote_data
def test_download():
    # Check that clipping the max splits at 3 reliably works for CDAWeb
    # This query has 4 small files
    trange = a.Time('2021/07/01', '2021/07/05')
    dataset = a.cdaweb.Dataset('SOLO_L2_MAG-RTN-NORMAL-1-MINUTE')
    result = Fido.search(trange, dataset)
    assert len(result[0]) == 4
    files = Fido.fetch(result)
    assert len(files.errors) == 0


@pytest.mark.remote_data
def test_no_results():
    result = Fido.search(
        a.Time('2000/03/01', '2000/03/02'),
        a.cdaweb.Dataset('SOLO_L2_MAG-RTN-NORMAL-1-MINUTE')
    )
    assert len(result[0]) == 0
