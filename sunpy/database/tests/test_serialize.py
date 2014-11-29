import json

from astropy import units as u

from sunpy.net import vso
from sunpy.database import attrs as db_attrs
from sunpy.database.serialize import QueryEncoder, query_decode


def test_vso_wave():
    attr = vso.attrs.Wave(100 * u.AA, 200 * u.AA)
    expected = '{"Wave": [100.0, 200.0, "Angstrom"]}'
    assert json.dumps(attr, cls=QueryEncoder) == expected


def test_vso_time():
    attr = vso.attrs.Time((1999, 12, 31), (2000, 1, 1))
    expected = '{"Time": ["1999-12-31 00:00:00", "2000-01-01 00:00:00", null]}'
    assert json.dumps(attr, cls=QueryEncoder) == expected


def test_vso_simple_attr():
    attr = vso.attrs.Instrument('EIT')
    expected = '{"Instrument": "EIT"}'
    assert json.dumps(attr, cls=QueryEncoder) == expected


def test_starred():
    expected = '{"Starred": true}'
    assert json.dumps(db_attrs.Starred(), cls=QueryEncoder) == expected


def test_starred_inverted():
    expected = '{"Starred": false}'
    assert json.dumps(~db_attrs.Starred(), cls=QueryEncoder) == expected


def test_tag():
    expected = '{"Tag": ["foo", false]}'
    assert json.dumps(db_attrs.Tag('foo'), cls=QueryEncoder) == expected


def test_tag_inverted():
    expected = '{"Tag": ["foo", true]}'
    assert json.dumps(~db_attrs.Tag('foo'), cls=QueryEncoder) == expected


def test_path():
    expected = '{"Path": ["bar", false]}'
    assert json.dumps(db_attrs.Path('bar'), cls=QueryEncoder) == expected


def test_path_inverted():
    expected = '{"Path": ["bar", true]}'
    assert json.dumps(~db_attrs.Path('bar'), cls=QueryEncoder) == expected


def test_download_time():
    attr = db_attrs.DownloadTime((1991, 8, 25, 3, 15, 40), (2001, 3, 5))
    expected = (
        '{"DownloadTime": '
        '["1991-08-25 03:15:40", "2001-03-05 00:00:00", false]}')
    assert json.dumps(attr, cls=QueryEncoder) == expected


def test_download_time_inverted():
    attr = ~db_attrs.DownloadTime((1991, 8, 25, 3, 15, 40), (2001, 3, 5))
    expected = (
        '{"DownloadTime": '
        '["1991-08-25 03:15:40", "2001-03-05 00:00:00", true]}')
    assert json.dumps(attr, cls=QueryEncoder) == expected


def test_fits_header_entry():
    attr = db_attrs.FitsHeaderEntry('key', 'value')
    expected = '{"FitsHeaderEntry": ["key", "value", false]}'
    assert json.dumps(attr, cls=QueryEncoder) == expected


def test_fits_header_entry_inverted():
    attr = ~db_attrs.FitsHeaderEntry('key', 'value')
    expected = '{"FitsHeaderEntry": ["key", "value", true]}'
    assert json.dumps(attr, cls=QueryEncoder) == expected


def test_attr_or():
    attr = vso.attrs.Source('SOHO') | vso.attrs.Provider('SDAC')
    expected = '{"AttrOr": [{"Provider": "SDAC"}, {"Source": "SOHO"}]}'
    assert json.dumps(attr, cls=QueryEncoder) == expected


def test_attr_and():
    attr = vso.attrs.Source('SOHO') & vso.attrs.Provider('SDAC')
    expected = '{"AttrAnd": [{"Provider": "SDAC"}, {"Source": "SOHO"}]}'
    assert json.dumps(attr, cls=QueryEncoder) == expected


def test_decode_wave():
    dump = '{"Wave": [10.0, 20.0, "Angstrom"]}'
    assert json.loads(dump, object_hook=query_decode) == vso.attrs.Wave(10 * u.AA, 20 * u.AA)


def test_decode_time():
    dump = '{"Time": ["1999-12-31 00:00:00", "2000-01-01 00:00:00", null]}'
    expected = vso.attrs.Time((1999, 12, 31), (2000, 1, 1))
    assert json.loads(dump, object_hook=query_decode) == expected


def test_decode_simple_attr():
    dump = '{"Instrument": "EIT"}'
    expected = vso.attrs.Instrument('EIT')
    assert json.loads(dump, object_hook=query_decode) == expected


def test_decode_starred():
    dump = '{"Starred": false}'
    assert json.loads(dump, object_hook=query_decode) == db_attrs.Starred()


def test_decode_starred_inverted():
    dump = '{"Starred": true}'
    assert json.loads(dump, object_hook=query_decode) == ~db_attrs.Starred()


def test_decode_tag():
    dump = '{"Tag": ["foo", false]}'
    assert json.loads(dump, object_hook=query_decode) == db_attrs.Tag('foo')


def test_decode_path():
    dump = '{"Path": ["bar", false]}'
    assert json.loads(dump, object_hook=query_decode) == db_attrs.Path('bar')


def test_decode_download_time():
    dump = (
        '{"DownloadTime": '
        '["1991-08-25 03:15:40", "2001-03-05 00:00:00", true]}')
    expected = ~db_attrs.DownloadTime((1991, 8, 25, 3, 15, 40), (2001, 3, 5))
    assert json.loads(dump, object_hook=query_decode) == expected


def test_decode_fits_header_entry():
    dump = '{"FitsHeaderEntry": ["key", "value", false]}'
    expected = db_attrs.FitsHeaderEntry('key', 'value')
    assert json.loads(dump, object_hook=query_decode) == expected


def test_decode_or():
    dump = '{"AttrOr": [{"Source": "SOHO"}, {"Provider": "SDAC"}]}'
    expected = vso.attrs.Source('SOHO') | vso.attrs.Provider('SDAC')
    assert json.loads(dump, object_hook=query_decode) == expected


def test_decode_and():
    dump = '{"AttrAnd": [{"Source": "SOHO"}, {"Provider": "SDAC"}]}'
    expected = vso.attrs.Source('SOHO') & vso.attrs.Provider('SDAC')
    assert json.loads(dump, object_hook=query_decode) == expected
