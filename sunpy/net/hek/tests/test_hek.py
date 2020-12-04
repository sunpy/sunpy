import pytest

from sunpy.net import attr, attrs, hek


@pytest.fixture
def foostrwrap(request):
    return hek.attrs._StringParamAttrWrapper("foo")


@pytest.fixture
@pytest.mark.remote_data
def hek_client_creator():
    startTime = '2011/08/09 07:23:56'
    endTime = '2011/08/09 12:40:29'
    eventType = 'FL'

    hekTime = attrs.Time(startTime, endTime)
    hekEvent = attrs.hek.EventType(eventType)

    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    return hek_query


def test_eventtype_collide():
    with pytest.raises(TypeError):
        attrs.hek.AR & attrs.hek.CE
    with pytest.raises(TypeError):
        (attrs.hek.AR & attrs.Time((2011, 1, 1),
                                   (2011, 1, 2))) & attrs.hek.CE
        with pytest.raises(TypeError):
            (attrs.hek.AR | attrs.Time((2011, 1, 1),
                                       (2011, 1, 2))) & attrs.hek.CE


def test_eventtype_or():
    assert (attrs.hek.AR | attrs.hek.CE).item == "ar,ce"


def test_paramattr():
    res = hek.attrs.walker.create(hek.attrs._ParamAttr("foo", "=", "bar"), {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '=', 'param0': 'foo'}


def test_stringwrapper_eq(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap == "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '=', 'param0': 'foo'}


def test_stringwrapper_lt(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap < "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '<', 'param0': 'foo'}


def test_stringwrapper_gt(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap > "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '>', 'param0': 'foo'}


def test_stringwrapper_le(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap <= "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '<=', 'param0': 'foo'}


def test_stringwrapper_ge(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap >= "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '>=', 'param0': 'foo'}


def test_stringwrapper_ne(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap != "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': '!=', 'param0': 'foo'}


def test_stringwrapper_like(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap.like("bar"), {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'op0': 'like', 'param0': 'foo'}


def test_err_dummyattr_create():
    with pytest.raises(TypeError):
        hek.attrs.walker.create(attr.DummyAttr(), {})


def test_err_dummyattr_apply():
    with pytest.raises(TypeError):
        hek.attrs.walker.apply(attr.DummyAttr(), {})


@pytest.mark.remote_data
def test_hek_client():
    startTime = '2011/08/09 07:23:56'
    endTime = '2011/08/09 12:40:29'
    eventType = 'FL'

    hekTime = attrs.Time(startTime, endTime)
    hekEvent = attrs.hek.EventType(eventType)

    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    assert type(hek_query) == hek.hek.HEKResponse


@pytest.mark.remote_data
def test_hek_empty_search_result():
    startTime = '1985-05-04 00:00:00'
    endTime = '1985-05-04 00:00:00'
    eventType = 'FL'

    hekTime = attrs.Time(startTime, endTime)
    hekEvent = attrs.hek.EventType(eventType)

    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    assert type(hek_query) == hek.hek.HEKResponse
    assert len(hek_query) == 0


@pytest.mark.remote_data
def test_getitem(hek_client_creator):
    hc = hek_client_creator
    assert hc.__getitem__(0) == hc[0]


@pytest.mark.remote_data
def test_get_voevent(hek_client_creator):
    hc = hek_client_creator
    ve = hc[0].get_voevent()
    assert len(ve['voe:VOEvent']) == 7


@pytest.mark.remote_data
def test_vso_time(hek_client_creator):
    hc = hek_client_creator
    ve = hc[0].vso_time
    assert type(ve) == attrs.Time


@pytest.mark.remote_data
def test_vso_instrument(hek_client_creator):
    hc = hek_client_creator
    try:
        vc = hc[1].vso_instrument
        assert type(vc) == attrs.Instrument
    except ValueError:
        assert 1


@pytest.mark.remote_data
def test_HEKRow_get(hek_client_creator):
    hc = hek_client_creator
    assert hc[0]['event_peaktime'] == hc[0].get('event_peaktime')
    assert hc[0].get('') is None


@pytest.mark.remote_data
def test_mixed_results_get():
    # To check that the following bug is fixed:
    # https://github.com/sunpy/sunpy/issues/3238
    client = hek.HEKClient()
    result = client.search(attrs.Time('2013/02/01 00:00:00', '2013/02/01 23:30:00'),
                           attrs.hek.FRM.Name == 'SPoCA')
    assert isinstance(result, hek.hek.HEKResponse)
    assert len(result) == 89
    assert result[0]["SOL_standard"] == 'SOL2013-01-31T20:13:31L219C160'


@pytest.mark.remote_data
def test_mixed_results_get_2():
    # To check that the following bug is fixed:
    # # https://github.com/sunpy/sunpy/issues/3898
    client = hek.HEKClient()
    result = client.search(attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'),
                           attrs.hek.EventType("FL"))
    assert isinstance(result, hek.hek.HEKResponse)
    assert len(result) == 19
    assert result[0]["SOL_standard"] == 'SOL2011-08-08T01:30:04L247C075'


@pytest.mark.remote_data
def test_mixed_results_get_angstrom():
    # To check that the following bug is fixed:
    # https://github.com/sunpy/sunpy/issues/4087
    client = hek.HEKClient()
    tstart = '2014/10/24 20:50'
    tend = '2014/10/25 00:14'
    event_type = 'FL'
    result = client.search(attrs.Time(tstart, tend), attrs.hek.EventType(event_type))
    assert len(result) == 13
    assert result[0]["SOL_standard"] == 'SOL2014-10-24T20:53:46L247C106'
