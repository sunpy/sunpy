import os
import copy
import pathlib
from unittest import mock

import hypothesis.strategies as st
import pytest
from drms import DrmsQueryError
from hypothesis import assume, given, settings
from parfive import Results
from parfive.utils import FailedDownload

import astropy.units as u

from sunpy import config
from sunpy.net import Fido, attr
from sunpy.net import attrs as a
from sunpy.net.base_client import BaseClient
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.sources.goes import XRSClient
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.tests.strategies import goes_time, offline_instruments, online_instruments, time_attr
from sunpy.net.vso import QueryResponse as vsoQueryResponse
from sunpy.net.vso.vso import DownloadFailed
from sunpy.time import TimeRange, parse_time
from sunpy.util.datatype_factory_base import MultipleMatchError

TIMEFORMAT = config.get("general", "time_format")


@st.composite
def offline_query(draw, instrument=offline_instruments()):
    """
    Strategy for any valid offline query
    """
    query = draw(instrument)
    # If we have AttrAnd then we don't have GOES
    if isinstance(query, a.Instrument) and query.value == 'goes':
        query &= draw(goes_time())
    else:
        query = attr.and_(query, draw(time_attr()))
    return query


@st.composite
def online_query(draw, instrument=online_instruments()):
    query = draw(instrument)

    if isinstance(query, a.Instrument) and query.value == 'eve':
        query &= a.Level(0)
    if isinstance(query, a.Instrument) and query.value == 'norh':
        query &= a.Wavelength(17*u.GHz)

    return query


@settings(deadline=50000)
@given(offline_query())
def test_offline_fido(query):
    unifiedresp = Fido.search(query)
    check_response(query, unifiedresp)


@pytest.mark.remote_data
@pytest.mark.flaky(reruns=5)
# Until we get more mocked, we can't really do this to online clients.
# TODO: Hypothesis this again
@pytest.mark.parametrize("query", [
    (a.Instrument('eve') & a.Time('2014/7/7', '2014/7/14') & a.Level(0)),
    (a.Instrument('rhessi') & a.Time('2014/7/7', '2014/7/14')),
    (a.Instrument('norh') & a.Time('2014/7/7', '2014/7/14') & a.Wavelength(17*u.GHz)),
])
def test_online_fido(query):
    unifiedresp = Fido.search(query)
    check_response(query, unifiedresp)


def check_response(query, unifiedresp):
    """
    Common test for online or offline query
    """
    query_tr = None
    query_instr = None
    for at in query.attrs:
        if isinstance(at, a.Time):
            query_tr = TimeRange(at.start, at.end)
        elif isinstance(at, a.Instrument):
            query_instr = at.value
    if not query_tr:
        raise ValueError("No Time Specified")

    for block in unifiedresp.responses:
        res_tr = block.time_range()
        for res in block:
            assert res.time.start in res_tr
            assert query_instr.lower() == res.instrument.lower()


@pytest.mark.remote_data
def test_save_path(tmpdir):
    qr = Fido.search(a.Instrument('EVE'), a.Time("2016/10/01", "2016/10/02"), a.Level(0))

    # Test when path is str
    files = Fido.fetch(qr, path=str(tmpdir / "{instrument}" / "{level}"))
    for f in files:
        assert str(tmpdir) in f
        assert "eve{}0".format(os.path.sep) in f


@pytest.mark.remote_data
@pytest.mark.flaky(reruns=5)
def test_save_path_pathlib(tmpdir):
    qr = Fido.search(a.Instrument('EVE'), a.Time("2016/10/01", "2016/10/02"), a.Level(0))

    # Test when path is pathlib.Path
    target_dir = tmpdir.mkdir("down")
    path = pathlib.Path(target_dir, "{instrument}", "{level}")
    files = Fido.fetch(qr, path=path)
    for f in files:
        assert target_dir.strpath in f
        assert "eve{}0".format(os.path.sep) in f


"""
Factory Tests
"""


@pytest.mark.remote_data
def test_unified_response():
    start = parse_time("2012/1/1")
    end = parse_time("2012/1/2")
    qr = Fido.search(a.Instrument('EVE'), a.Level(0), a.Time(start, end))
    assert qr.file_num == 2
    strings = ['eve', 'SDO', start.strftime(TIMEFORMAT), end.strftime(TIMEFORMAT)]
    assert all(s in qr._repr_html_() for s in strings)


def test_no_time_error():
    # TODO: work out why this is failing sometimes.
    from sunpy.net import attrs as a
    query = (a.Instrument('EVE'), a.Level(0))
    with pytest.raises(ValueError) as excinfo:
        Fido.search(*query)
    assert all(str(a) in str(excinfo.value) for a in query)

    query1 = (a.Instrument('EVE') & a.Level(0))
    query2 = (a.Time("2012/1/1", "2012/1/2") & a.Instrument("AIA"))
    with pytest.raises(ValueError) as excinfo:
        Fido.search(query1 | query2)
    assert all(str(a) in str(excinfo.value) for a in query1.attrs)
    assert all(str(a) not in str(excinfo.value) for a in query2.attrs)


@pytest.mark.remote_data
def test_no_match():
    with pytest.raises(DrmsQueryError):
        Fido.search(a.Time("2016/10/01", "2016/10/02"), a.jsoc.Series("bob"),
                    a.vso.Sample(10*u.s))


def test_call_error():
    with pytest.raises(TypeError) as excinfo:
        Fido()
    # Explicitly test all this error message as it's a copy of the one in
    # Python core.
    assert "'UnifiedDownloaderFactory' object is not callable" in str(excinfo.value)


def test_multiple_match():
    """
    Using the builtin clients a multiple match is not possible so we create a
    dummy class.
    """
    new_registry = copy.deepcopy(Fido.registry)
    Fido.registry = new_registry

    class DummyClient():
        @classmethod
        def _can_handle_query(cls, *query):
            return True
    Fido.registry.update({DummyClient: DummyClient._can_handle_query})

    with pytest.raises(MultipleMatchError):
        Fido.search(a.Time("2016/10/1", "2016/10/2"), a.Instrument('lyra'))

    Fido.registry = BaseClient._registry


@pytest.mark.remote_data
def test_fetch():
        qr = Fido.search(a.Instrument('EVE'),
                         a.Time("2016/10/01", "2016/10/02"),
                         a.Level(0))
        res = Fido.fetch(qr)
        assert isinstance(res, Results)


"""
UnifiedResponse Tests

Use LYRA here because it does not use the internet to return results.
"""


def test_unifiedresponse_slicing():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument("lyra"))
    assert isinstance(results[0:2], UnifiedResponse)
    assert isinstance(results[0], UnifiedResponse)


def test_unifiedresponse_slicing_reverse():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument("lyra"))
    assert isinstance(results[::-1], UnifiedResponse)
    assert len(results[::-1]) == len(results)
    assert isinstance(results[0, ::-1], UnifiedResponse)
    assert results[0, ::-1]._list[0] == results._list[0][::-1]


def test_vso_unifiedresponse():
    vrep = vsoQueryResponse([])
    vrep.client = True
    uresp = UnifiedResponse(vrep)
    assert isinstance(uresp, UnifiedResponse)


def test_responses():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument("lyra"))

    for i, resp in enumerate(results.responses):
        assert isinstance(resp, QueryResponse)

    assert i + 1 == len(results)


def test_repr():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument("lyra"))

    rep = repr(results)
    rep = rep.split('\n')
    # 6 header lines, the results table and two blank lines at the end
    assert len(rep) == 7 + len(list(results.responses)[0]) + 2


def filter_queries(queries):
    return attr.and_(queries) not in queries


@settings(deadline=50000)
@given(st.tuples(offline_query(), offline_query()).filter(filter_queries))
def test_fido_indexing(queries):
    query1, query2 = queries

    # This is a work around for an aberration where the filter was not catching
    # this.
    assume(query1.attrs[1].start != query2.attrs[1].start)

    res = Fido.search(query1 | query2)

    assert len(res) == 2
    assert len(res[0]) == 1
    assert len(res[1]) == 1

    aa = res[0, 0]
    assert isinstance(aa, UnifiedResponse)
    assert len(aa) == 1
    assert len(aa.get_response(0)) == 1

    aa = res[:, 0]
    assert isinstance(aa, UnifiedResponse)
    assert len(aa) == 2
    assert len(aa.get_response(0)) == 1

    aa = res[0, :]
    assert isinstance(aa, UnifiedResponse)
    assert len(aa) == 1

    with pytest.raises(IndexError):
        res[0, 0, 0]

    with pytest.raises(IndexError):
        res["saldkal"]

    with pytest.raises(IndexError):
        res[1.0132]


@settings(deadline=50000)
@given(st.tuples(offline_query(), offline_query()).filter(filter_queries))
def test_fido_iter(queries):
    query1, query2 = queries

    # This is a work around for an aberration where the filter was not catching
    # this.
    assume(query1.attrs[1].start != query2.attrs[1].start)

    res = Fido.search(query1 | query2)

    for resp in res:
        assert isinstance(resp, QueryResponse)


@settings(deadline=50000)
@given(offline_query())
def test_repr2(query):
    res = Fido.search(query)

    for rep_meth in (res.__repr__, res.__str__, res._repr_html_):
        if len(res) == 1:
            assert "Provider" in rep_meth()
            assert "Providers" not in rep_meth()

        else:
            assert "Provider" not in rep_meth()
            assert "Providers" in rep_meth()


@mock.patch("parfive.Downloader.download", return_value=Results(["/tmp/test"]))
def test_retry(mock_retry):
    """
    Test that you can use Fido.fetch to retry failed downloads.
    """
    res = Results()
    res.data.append("/this/worked.fits")

    err1 = FailedDownload("This is not a filename", "http://not.url/test", None)
    err2 = FailedDownload("This is not a filename2", "http://not.url/test2", None)
    res.errors.append(err1)
    res.errors.append(err2)

    mock_retry.return_value._errors += [err2]

    res2 = Fido.fetch(res, Results(["/this/also/worked.fits"]))

    assert res2 is not res

    # Assert that the result of retry ends up in the returned Results() object
    assert res2.data == ["/this/worked.fits", "/tmp/test", "/this/also/worked.fits", "/tmp/test"]
    assert res2.errors == [err2, err2]


def results_generator(dl):
    http = list(dl.http_queue._queue)
    ftp = list(dl.ftp_queue._queue)

    outputs = []
    for url in http+ftp:
        outputs.append(pathlib.Path(url.keywords['url'].split("/")[-1]))

    return Results(outputs)


@pytest.mark.remote_data
@mock.patch("sunpy.net.vso.VSOClient.download_all",
            return_value=Results([], errors=[DownloadFailed(None)]))
@mock.patch("parfive.Downloader.download", new=results_generator)
def test_vso_errors_with_second_client(mock_download_all):
    query = a.Time("2011/01/01", "2011/01/02") & (a.Instrument("goes") | a.Instrument("EIT"))

    qr = Fido.search(query)

    res = Fido.fetch(qr)
    assert len(res.errors) == 1
    assert len(res) != qr.file_num

    # Assert that all the XRSClient records are in the output.
    for resp in qr.responses:
        if isinstance(resp, XRSClient):
            assert len(resp) == len(res)


def test_downloader_type_error():
    with pytest.raises(TypeError):
        Fido.fetch([], downloader=Results())


def test_mixed_retry_error():
    with pytest.raises(TypeError):
        Fido.fetch([], Results())


@mock.patch("sunpy.net.dataretriever.sources.goes.XRSClient.fetch",
            return_value=["hello"])
def test_client_fetch_wrong_type(mock_fetch):
    query = a.Time("2011/01/01", "2011/01/02") & a.Instrument("goes")

    qr = Fido.search(query)

    with pytest.raises(TypeError):
        Fido.fetch(qr)
