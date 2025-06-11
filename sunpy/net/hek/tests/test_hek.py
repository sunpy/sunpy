
import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import MaskedColumn, Table
from astropy.time import Time

from sunpy.coordinates import Helioprojective, get_earth
from sunpy.net import Fido, attr, attrs, hek
from sunpy.net.hek.utils import _get_coord_attributes, _get_unit_attributes


@pytest.fixture
def foostrwrap(request):
    return hek.attrs._StringParamAttrWrapper("foo")


@pytest.fixture(scope="session")
def hek_result():
    startTime = '2011/08/09 07:23:56'
    endTime = '2011/08/09 12:40:29'
    eventType = 'FL'
    hekTime = attrs.Time(startTime, endTime)
    hekEvent = attrs.hek.EventType(eventType)
    h = hek.HEKClient()
    return h.search(hekTime, hekEvent)


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


def test_HEKAttr():
    res = hek.attrs.walker.create(hek.attrs.HEKAttr("foo", "=", "bar"), {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'operator0': '=', 'param0': 'foo'}


def test_stringwrapper_eq(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap == "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'operator0': '=', 'param0': 'foo'}


def test_stringwrapper_lt(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap < "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'operator0': '<', 'param0': 'foo'}


def test_stringwrapper_gt(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap > "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'operator0': '>', 'param0': 'foo'}


def test_stringwrapper_le(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap <= "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'operator0': '<=', 'param0': 'foo'}


def test_stringwrapper_ge(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap >= "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'operator0': '>=', 'param0': 'foo'}


def test_stringwrapper_ne(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap != "bar", {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'operator0': '!=', 'param0': 'foo'}


def test_stringwrapper_like(foostrwrap):
    res = hek.attrs.walker.create(foostrwrap.like("bar"), {})
    assert len(res) == 1
    assert res[0] == {'value0': 'bar', 'operator0': 'like', 'param0': 'foo'}


def test_err_dummyattr_create():
    with pytest.raises(TypeError):
        hek.attrs.walker.create(attr.DummyAttr(), {})


def test_err_dummyattr_apply():
    with pytest.raises(TypeError):
        hek.attrs.walker.apply(attr.DummyAttr(), {})


@pytest.mark.remote_data
def test_hek_client(hek_result):
    assert isinstance(hek_result, hek.hek.HEKTable)


@pytest.mark.remote_data
def test_hek_empty_search_result():
    startTime = '1985-05-04 00:00:00'
    endTime = '1985-05-04 00:00:00'
    eventType = 'FL'

    hekTime = attrs.Time(startTime, endTime)
    hekEvent = attrs.hek.EventType(eventType)

    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    assert isinstance(hek_query, hek.hek.HEKTable)
    assert len(hek_query) == 0


@pytest.mark.remote_data
def test_get_voevent(hek_result):
    ve = hek_result[0].get_voevent()
    assert len(ve['voe:VOEvent']) == 7


@pytest.mark.remote_data
def test_hek_time_col(hek_result):
    assert isinstance(hek_result[0]['event_starttime'], Time)
    assert isinstance(hek_result[0]['event_endtime'], Time)


@pytest.mark.remote_data
def test_vso_time(hek_result):
    ve = hek_result[0].vso_time
    assert isinstance(ve, attrs.Time)


@pytest.mark.remote_data
def test_vso_instrument(hek_result):
    vc = hek_result[1].vso_instrument
    assert isinstance(vc, attrs.Instrument)


@pytest.mark.remote_data
def test_HEKRow_get(hek_result):
    assert hek_result[0]['event_peaktime'] == hek_result[0].get('event_peaktime')
    assert hek_result[0].get('') is None


@pytest.mark.remote_data
def test_mixed_results_get():
    # To check that the following bug is fixed:
    # https://github.com/sunpy/sunpy/issues/3238
    client = hek.HEKClient()
    result = client.search(attrs.Time('2013/02/01 00:00:00', '2013/02/01 23:30:00'),
                           attrs.hek.FRM.Name == 'SPoCA')
    assert isinstance(result, hek.hek.HEKTable)
    assert len(result) == 89
    # We do not check the full timestamp as the last 8 digits change as data is reprocessed.
    assert result[0]["SOL_standard"].startswith("SOL2013-01-31T20:13:31")


@pytest.mark.remote_data
def test_mixed_results_get_2():
    # To check that the following bug is fixed:
    # # https://github.com/sunpy/sunpy/issues/3898
    client = hek.HEKClient()
    result = client.search(attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'),
                           attrs.hek.EventType("FL"))
    assert isinstance(result, hek.hek.HEKTable)
    assert len(result) == 19
    # We do not check the full timestamp as the last 8 digits change as data is reprocessed.
    assert result[0]["SOL_standard"].startswith("SOL2011-08-08T01:30:04")


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


@pytest.mark.remote_data
def test_query_multiple_operators():
    event_type = "FL"
    tstart = "2013/10/28"
    tend = "2013/10/29"
    client = hek.HEKClient()
    results = client.search(attrs.Time(tstart, tend),
                            attrs.hek.EventType(event_type),
                            attrs.hek.FL.GOESCls > "M1.0",
                            attrs.hek.OBS.Observatory == "GOES")
    assert len(results) == 7


@pytest.mark.remote_data
def test_astropy_unit_parsing(coronal_hole_search_result):
    hek_attributes = _get_unit_attributes() + _get_coord_attributes()
    attributes_with_unit = [attribute for attribute in hek_attributes if attribute.get("unit_prop", False)]
    for attribute in attributes_with_unit:
        if attribute["name"] in coronal_hole_search_result.colnames and not attribute.get("is_chaincode", False):
            assert isinstance(coronal_hole_search_result[attribute['name']], u.Quantity)


@pytest.mark.remote_data
def test_chaincode_parsing(coronal_hole_search_result):
    chaincode_attributes = [attribute for attribute in _get_coord_attributes() if attribute.get("is_chaincode", False)]
    for attribute in chaincode_attributes:
        if attribute["name"] in coronal_hole_search_result.colnames:
            if isinstance(coronal_hole_search_result[attribute['name']], MaskedColumn):
                assert all([isinstance(r, SkyCoord) for r in coronal_hole_search_result[attribute['name']]])
            else:
                assert isinstance(coronal_hole_search_result[attribute['name']], SkyCoord)


@pytest.mark.remote_data
def test_missing_times():
    # Check for https://github.com/sunpy/sunpy/pull/7627#issuecomment-2113451964
    client = hek.HEKClient()
    results = client.search(attrs.Time('2024-05-10', '2024-05-12'), attrs.hek.AR.NOAANum == 13664)
    assert isinstance(results["event_peaktime"][0], np.ma.core.MaskedConstant)
    assert results["event_peaktime"][3].isot == "2024-05-10T00:13:00.000"


@pytest.fixture(scope="session")
def coronal_hole_search_result():
    client = hek.HEKClient()
    result = client.search(attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'),
                           attrs.hek.EventType('CH'))
    return result


@pytest.mark.remote_data
def test_compare_event_coords(coronal_hole_search_result):
    event_coord = SkyCoord(-2.91584*u.arcsec,
                           940.667*u.arcsec,
                           frame=Helioprojective(observer=get_earth('2011-08-09 06:00:08.000')))
    assert coronal_hole_search_result['event_coord'][0] == event_coord


@pytest.mark.remote_data
def test_obs_meanwavel(coronal_hole_search_result):
    assert coronal_hole_search_result['obs_meanwavel'][0] == 193.0*u.angstrom


@pytest.fixture
def flare_search():
    return attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'), attrs.hek.EventType('FL')


@pytest.mark.remote_data
def test_ssw_latest_events_flares(flare_search):
    result = Fido.search(*flare_search, attrs.hek.FRM.Name == 'SSW Latest Events')
    assert len(result[0]) == 2


@pytest.mark.remote_data
def test_not_ssw_latest_events_flares(flare_search):
    result = Fido.search(*flare_search, attrs.hek.FRM.Name != 'SSW Latest Events')
    assert len(result[0]) == 19


@pytest.mark.remote_data
def test_flares_peak_flux(flare_search):
    result = Fido.search(*flare_search, attrs.hek.FL.PeakFlux > 4000.0)
    assert len(result[0]) == 1
    assert result[0]['fl_peakflux'] > 4000.0*u.DN/(u.pix*u.s)


@pytest.mark.remote_data
def test_flares_peak_flux_and_position(flare_search):
    result = Fido.search(*flare_search,
                         attrs.hek.Event.Coord1 > 800,
                         attrs.hek.FL.PeakFlux > 1000)
    assert len(result[0]) == 7
    assert all([c.Tx > 800*u.arcsec for c in result[0]['event_coord']])
    assert (result[0]['fl_peakflux'] > 1000.0*u.DN/(u.pix*u.s)).all()


@pytest.mark.remote_data
def test_flares_python_logical_ops(flare_search):
    result = Fido.search(*flare_search,
                         (attrs.hek.Event.Coord1 > 50) and (attrs.hek.FL.PeakFlux > 1000))
    assert len(result[0]) == 7
    assert all([c.Tx > 50*u.arcsec for c in result[0]['event_coord']])
    assert (result[0]['fl_peakflux'] > 1000.0*u.DN/(u.pix*u.s)).all()


@pytest.mark.remote_data
@pytest.mark.parametrize("event_type", [
    'AR',
    'CE',
    'CD',
    'CW',
    'FI',
    'FE',
    'FA',
    'LP',
    'OS',
    'SS',
    'EF',
    'CJ',
    'PG',
    'OT',
    'NR',
    'SG',
    'SP',
    'CR',
    'CC',
    'ER',
    'TO',
    'HY',
    'BU',
    'EE',
    'PB',
    'PT',
])
def test_event_types(event_type):
    # Smoke test for all event types
    _ = Fido.search(attrs.Time('2017/09/06 11:59:04', '2017/09/06 17:05:04'), attrs.hek.EventType(event_type))


@pytest.mark.remote_data
def test_raw_hek_result_preserved(hek_result):
    assert hasattr(hek_result, 'raw')
    assert isinstance(hek_result.raw, Table)
    # Check that separate event_coord columns are still present
    removed_columns = ['event_coord1', 'event_coord2', 'event_coord3', 'event_coordsys']
    for col in removed_columns:
        assert col in hek_result.raw.colnames
    # Check that times are still strings
    assert np.issubdtype(hek_result.raw['event_starttime'].dtype, np.str_)
    assert np.issubdtype(hek_result.raw['event_endtime'].dtype, np.str_)
    # Check that chaincodes are still strings
    for coord_attr in _get_coord_attributes():
        if not coord_attr.get('is_chaincode', False):
            continue
        assert np.issubdtype(hek_result.raw[coord_attr['name']].dtype, np.dtype('U'))
    # Check that quantities are still floats and that unit columns still exist
    for unit_attr in _get_unit_attributes():
        if (name := unit_attr['name']) not in hek_result.raw.colnames:
            continue
        column_dtype = hek_result.raw[name].dtype
        if unit_attr.get('unit_prop') is not None:
            # NOTE: Columns have object dtype if there are Nones present
            assert np.issubdtype(column_dtype, np.float64) | np.issubdtype(column_dtype, np.object_)
        elif unit_attr.get('is_unit_prop', False):
            assert np.issubdtype(column_dtype, np.str_)
