import os
import pathlib
from stat import S_IREAD, S_IRGRP, S_IROTH
from unittest import mock
from collections import OrderedDict

import hypothesis.strategies as st
import pytest
from drms import DrmsQueryError
from hypothesis import assume, given, settings
from parfive import Results
from parfive.utils import FailedDownload

import astropy.units as u
from astropy.table import Table

from sunpy import config
from sunpy.net import Fido, attr
from sunpy.net import attrs as a
from sunpy.net import jsoc
from sunpy.net.base_client import BaseQueryResponse, BaseQueryResponseTable
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.sources.goes import XRSClient
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.tests.strategies import goes_time, offline_instruments, online_instruments, srs_time, time_attr
from sunpy.net.vso import QueryResponse as vsoQueryResponse
from sunpy.net.vso.vso import DownloadFailed
from sunpy.tests.helpers import no_vso, skip_windows
from sunpy.time import TimeRange, parse_time
from sunpy.util.exceptions import SunpyUserWarning

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
    elif isinstance(query, a.Instrument) and query.value == 'soon':
        query &= draw(srs_time())
    else:
        query = attr.and_(query, draw(time_attr()))
    return query


@st.composite
def online_query(draw, instrument=online_instruments()):
    query = draw(instrument)

    if isinstance(query, a.Instrument) and query.value == 'eve':
        query &= a.Level.zero
    if isinstance(query, a.Instrument) and query.value == 'norh':
        query &= a.Wavelength(17*u.GHz)

    return query


@no_vso
@settings(deadline=50000)
@given(offline_query())
def test_offline_fido(query):
    unifiedresp = Fido.search(query)
    check_response(query, unifiedresp)


@pytest.mark.remote_data
# Until we get more mocked, we can't really do this to online clients.
# TODO: Hypothesis this again
@pytest.mark.parametrize("query", [
    (a.Instrument.eve & a.Time('2014/7/7', '2014/7/14') & a.Level.zero),
    (a.Instrument.rhessi & a.Time('2014/7/7', '2014/7/14')),
    (a.Instrument.norh & a.Time('2014/7/7', '2014/7/14') & a.Wavelength(17*u.GHz)),
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
            if isinstance(res, OrderedDict) or isinstance(res, dict):
                assert res['Time'].start in res_tr
                assert query_instr.lower() == res['Instrument'].lower()
            else:
                assert res.time.start in res_tr
                assert query_instr.lower() == res.instrument.lower()


@pytest.mark.remote_data
def test_save_path(tmpdir):
    qr = Fido.search(a.Instrument.eve, a.Time("2016/10/01", "2016/10/02"), a.Level.zero)

    # Test when path is str
    files = Fido.fetch(qr, path=str(tmpdir / "{Instrument}" / "{Level}"))
    for f in files:
        assert str(tmpdir) in f
        assert f"EVE{os.path.sep}0" in f


@pytest.mark.remote_data
def test_save_path_pathlib(tmpdir):
    qr = Fido.search(a.Instrument.eve, a.Time("2016/10/01", "2016/10/02"), a.Level.zero)

    # Test when path is pathlib.Path
    target_dir = tmpdir.mkdir("down")
    path = pathlib.Path(target_dir, "{Instrument}", "{Level}")
    files = Fido.fetch(qr, path=path)
    for f in files:
        assert target_dir.strpath in f
        assert f"EVE{os.path.sep}0" in f


@pytest.mark.remote_data
def test_save_path_cwd(tmpdir):
    qr = Fido.search(a.Instrument.eve, a.Time("2016/10/01", "2016/10/02"), a.Level.zero)

    # Test when path is ./ for current working directory
    os.chdir(tmpdir)  # move into temp directory
    files = Fido.fetch(qr, path="./")
    for f in files:
        assert pathlib.Path.cwd().joinpath(f).exists()


"""
Factory Tests
"""


@pytest.mark.remote_data
def test_unified_response():
    start = parse_time("2012/1/1")
    end = parse_time("2012/1/2")
    qr = Fido.search(a.Instrument.eve, a.Level.zero, a.Time(start, end))
    assert qr.file_num == 2
    strings = ['eve', 'SDO', start.strftime(TIMEFORMAT), end.strftime(TIMEFORMAT)]
    assert all(s in qr._repr_html_() for s in strings)


@pytest.mark.remote_data
def test_no_match():
    with pytest.raises(DrmsQueryError):
        Fido.search(a.Time("2016/10/01", "2016/10/02"), a.jsoc.Series("bob"),
                    a.Sample(10*u.s))


def test_call_error():
    with pytest.raises(TypeError) as excinfo:
        Fido()
    # Explicitly test all this error message as it's a copy of the one in
    # Python core.
    assert "'UnifiedDownloaderFactory' object is not callable" in str(excinfo.value)


@pytest.mark.remote_data
def test_fetch():
    qr = Fido.search(a.Instrument.eve,
                     a.Time("2016/10/01", "2016/10/02"),
                     a.Level.zero)
    res = Fido.fetch(qr)
    assert isinstance(res, Results)


"""
UnifiedResponse Tests

Use LYRA here because it does not use the internet to return results.
"""


@pytest.mark.remote_data
def test_unifiedresponse_slicing():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument.lyra)
    assert isinstance(results[0:2], BaseQueryResponse)
    assert isinstance(results[0], BaseQueryResponse)


@pytest.mark.remote_data
def test_unifiedresponse_slicing_reverse():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument.lyra)
    assert isinstance(results[::-1], BaseQueryResponse)
    assert len(results[::-1]) == len(results[::1])
    assert isinstance(results[0, ::-1], BaseQueryResponse)
    assert all(results[0][::-1].build_table() == results[0, ::-1].build_table())


@pytest.mark.remote_data
def test_tables_single_response():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument.lyra, a.Level.two)
    tables = results.tables

    assert isinstance(tables, list)
    assert isinstance(tables[0], Table)
    assert len(tables) == 1

    columns = ['Start Time', 'End Time', 'Instrument', 'Physobs',
               'Source', 'Provider', 'Level']
    assert columns == tables[0].colnames
    assert len(tables[0]) == 5


@pytest.mark.remote_data
def test_tables_multiple_response():
    results = Fido.search(a.Time('2012/3/4', '2012/3/6'),
                          a.Instrument.lyra | (a.Instrument.rhessi & a.Physobs.summary_lightcurve))
    tables = results.tables
    assert isinstance(tables, list)
    assert all(isinstance(t, Table) for t in tables)
    assert len(tables) == 2

    columns0 = ['Start Time', 'End Time', 'Instrument', 'Physobs',
                'Source', 'Provider', 'Level']
    columns1 = ['Start Time', 'End Time', 'Instrument', 'Physobs',
                'Source', 'Provider']
    assert columns0 == tables[0].colnames and columns1 == tables[1].colnames

    assert all(entry == 'LYRA' for entry in tables[0]['Instrument'])
    assert all(entry == 'RHESSI' for entry in tables[1]['Instrument'])


@pytest.mark.remote_data
def test_tables_all_types():
    # Data retriver response objects
    drclient = Fido.search(a.Time('2012/3/4', '2012/3/6'),
                           a.Instrument.lyra | (a.Instrument.rhessi & a.Physobs.summary_lightcurve))
    drtables = drclient.tables
    assert isinstance(drtables, list)
    assert isinstance(drtables[0], Table)

    # VSO response objects
    vsoclient = Fido.search(a.Time('2011-06-07 06:33', '2011-06-07 06:33:08'),
                            a.Instrument.aia, a.Wavelength(171 * u.AA))
    vsotables = vsoclient.tables
    assert isinstance(vsotables, list)
    assert isinstance(vsotables[0], Table)

    # JSOC response objects
    jsocclient = Fido.search(a.Time('2014-01-01T00:00:00', '2014-01-01T01:00:00'),
                             a.jsoc.Series('hmi.v_45s'), a.jsoc.Notify('sunpy@sunpy.org'))
    jsoctables = jsocclient.tables

    assert isinstance(jsoctables, list)
    assert isinstance(jsoctables[0], Table)


@mock.patch("sunpy.net.vso.vso.build_client", return_value=True)
def test_vso_unifiedresponse(mock_build_client):
    vrep = vsoQueryResponse([])
    vrep.client = True
    uresp = UnifiedResponse(vrep)
    assert isinstance(uresp, UnifiedResponse)


@pytest.mark.remote_data
def test_responses():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument.lyra)

    for i, resp in enumerate(results.responses):
        assert isinstance(resp, QueryResponse)

    assert i + 1 == len(results)


@pytest.mark.remote_data
def test_repr():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument.lyra)

    rep = repr(results)
    rep = rep.split('\n')
    # 6 header lines, the results table and two blank lines at the end
    assert len(rep) == 6 + len(list(results.responses)[0]) + 2


def filter_queries(queries):
    return attr.and_(queries) not in queries


@pytest.mark.remote_data
def test_path():
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument.lyra)

    Fido.fetch(results, path="notapath/{file}")


@pytest.mark.remote_data
@skip_windows
def test_path_read_only(tmp_path):
    results = Fido.search(
        a.Time("2012/1/1", "2012/1/5"), a.Instrument.lyra)

    # chmod dosen't seem to work correctly on the windows CI
    os.chmod(tmp_path, S_IREAD | S_IRGRP | S_IROTH)
    # Check to see if it's actually read only before running the test
    if not os.access(tmp_path, os.W_OK):
        with pytest.raises(PermissionError):
            Fido.fetch(results, path=tmp_path / "{file}")


@no_vso
@settings(deadline=50000)
@given(st.tuples(offline_query(), offline_query()).filter(filter_queries))
def test_fido_indexing(queries):
    query1, query2 = queries

    # This is a work around for an aberration where the filter was not catching
    # this.
    assume(query1.attrs[1].start != query2.attrs[1].start)

    res = Fido.search(query1 | query2)

    assert len(res) == 2
    assert len(res[0][0]) == 1
    assert len(res[1][0]) == 1

    aa = res[0, 0]
    assert isinstance(aa, BaseQueryResponse)
    assert len(aa) == 1
    assert len(aa[0]) == 1

    aa = res[:, 0]
    assert isinstance(aa, UnifiedResponse)
    assert len(aa) == 2
    assert len(aa.get_response(0)) == 1

    aa = res[0, :]
    assert isinstance(aa, BaseQueryResponse)
    assert len(aa[0]) == 1

    with pytest.raises(IndexError):
        res[0, 0, 0]

    with pytest.raises(IndexError):
        res["saldkal"]

    with pytest.raises(IndexError):
        res[1.0132]

    if isinstance(res, UnifiedResponse):
        assert len(res) != 1


@no_vso
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


@no_vso
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
    http = dl.http_queue
    ftp = dl.ftp_queue
    # Handle compatibility with parfive 1.0
    if not isinstance(dl.http_queue, list):
        http = list(dl.http_queue._queue)
        ftp = list(dl.ftp_queue._queue)

    outputs = []
    for url in http + ftp:
        outputs.append(pathlib.Path(url.keywords['url'].split("/")[-1]))

    return Results(outputs)


@pytest.mark.remote_data
@mock.patch("sunpy.net.vso.VSOClient.download_all",
            return_value=Results([], errors=[DownloadFailed(None)]))
@mock.patch("parfive.Downloader.download", new=results_generator)
def test_vso_errors_with_second_client(mock_download_all):
    query = a.Time("2011/01/01", "2011/01/02") & (a.Instrument.goes | a.Instrument.eit)
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


@pytest.mark.remote_data
@mock.patch("sunpy.net.dataretriever.sources.goes.XRSClient.fetch",
            return_value=["hello"])
def test_client_fetch_wrong_type(mock_fetch):
    query = a.Time("2011/01/01", "2011/01/02") & a.Instrument.goes
    qr = Fido.search(query)
    with pytest.raises(TypeError):
        Fido.fetch(qr)


@pytest.mark.remote_data
def test_vso_fetch_hmi(tmpdir):
    start_time = "2017-01-25"
    end_time = "2017-01-25T23:59:59"
    results = Fido.search(a.Time(start_time, end_time),
                          a.Instrument.hmi & a.Physobs.los_magnetic_field,
                          a.Sample(1 * u.minute))

    files = Fido.fetch(results[0, 0], path=tmpdir)
    assert len(files) == 1


@pytest.mark.remote_data
def test_unclosedSocket_warning():
    with pytest.warns(None):
        attrs_time = a.Time('2005/01/01 00:10', '2005/01/01 00:15')
        result = Fido.search(attrs_time, a.Instrument.eit)
        Fido.fetch(result)


def test_fido_no_time(mocker):
    jsoc_mock = mocker.patch("sunpy.net.jsoc.JSOCClient.search")
    jsoc_mock.return_value = jsoc.JSOCResponse()

    Fido.search(a.jsoc.Series("test"))

    jsoc_mock.assert_called_once()


@pytest.mark.remote_data
def test_slice_jsoc():
    tstart = '2011/06/07 06:32:45'
    tend = '2011/06/07 06:33:15'
    res = Fido.search(a.Time(tstart, tend), a.jsoc.Series('hmi.M_45s'),
                      a.jsoc.Notify('jsoc@cadair.com'))

    with pytest.warns(SunpyUserWarning):
        Fido.fetch(res[0, 0])

    with pytest.warns(SunpyUserWarning):
        Fido.fetch(res[0, 0:1])


def test_fido_repr():
    output = repr(Fido)
    assert output[:50] == '<sunpy.net.fido_factory.UnifiedDownloaderFactory o'


@pytest.mark.remote_data
def test_fido_metadata_queries():
    results = Fido.search(a.Time('2010/8/1 03:40', '2010/8/1 3:40:10'),
                          a.hek.FI | a.hek.FL & (a.hek.FL.PeakFlux > 1000) |
                          a.jsoc.Series('hmi.m_45s') & a.jsoc.Notify("jsoc@cadair.com"))

    assert len(results['hek']) == 2
    assert isinstance(results['hek'], UnifiedResponse)
    assert isinstance(results['hek'][0], BaseQueryResponse)
    assert len(results['hek'][1]) == 2
    assert results[::-1][0] == results['jsoc']
    assert isinstance(results['jsoc'], BaseQueryResponseTable)

    files = Fido.fetch(results)
    assert len(files) == len(results['jsoc'])
