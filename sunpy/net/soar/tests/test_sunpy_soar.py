from pathlib import Path

import astropy.units as u
import pytest
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.util.exceptions import SunpyDeprecationWarning

from sunpy.net.soar.client import SOARClient


def test_search():
    instrument = a.Instrument("EUI")
    time = a.Time("2021-02-11", "2021-02-12")
    level = a.Level(1)
    product = a.soar.Product("eui-fsi174-image")

    res = Fido.search(instrument, time, level, product)
    assert len(res) == 1
    assert len(res[0]) == 37
    assert u.allclose(res[0, 0]["Filesize"], 3.47 * u.Mbyte)

    # check passing upper case descriptor
    product = a.soar.Product("EUI-FSI174-IMAGE")
    res = Fido.search(time, level, product)
    assert res.file_num == 37

    files = Fido.fetch(res[0, 0])
    assert len(files) == 1
    fname = files[0]
    assert u.allclose(Path(fname).stat().st_size * u.byte, res[0, 0]["Filesize"], atol=1e-3 * u.Mbyte)
    # Smoke test that we can read this into a map
    sunpy.map.Map(fname)


def test_search_low_latency():
    time = a.Time("2020-11-13", "2020-11-14")
    level = a.Level("LL02")
    product = a.soar.Product("mag")

    res = Fido.search(time, level, product)
    assert len(res) == 1
    assert len(res[0]) == 1

    files = Fido.fetch(res[0, 0])
    assert len(files) == 1


def test_deprecated_identifier():
    instrument = a.Instrument("EUI")
    time = a.Time("2021-02-11", "2021-02-12")
    level = a.Level(1)
    with pytest.warns(SunpyDeprecationWarning):
        identifier = a.soar.Identifier("eui-fsi174-image")
    product = a.soar.Product("eui-fsi174-image")
    res1 = Fido.search(instrument, time, level, identifier)
    res2 = Fido.search(instrument, time, level, product)

    assert res1.__str__() == res2.__str__()


def test_insitu_search():
    instrument = a.Instrument("MAG")
    time = a.Time("2020-04-16", "2020-04-17")
    identifier = a.soar.Product("mag-rtn-normal-1-minute")

    res = Fido.search(instrument, time, identifier)
    assert len(res) == 1
    assert len(res[0]) == 2

    files = Fido.fetch(res[0, 0])
    assert len(files) == 1


def test_no_results():
    instrument = a.Instrument("EUI")
    time = a.Time("2019-02-01", "2019-02-02")
    query = instrument & time

    res = SOARClient().search(query)
    assert len(res) == 0


def test_no_instrument():
    # Check that a time only search returns results
    time = a.Time("2020-04-16", "2020-04-17")
    res = SOARClient().search(time)
    assert len(res) == 63


def test_download_path(tmp_path):
    # Check that we can download things to a custom path using
    # the search parameters
    instrument = a.Instrument("EUI")
    time = a.Time("2021-02-01", "2021-02-02")
    level = a.Level(1)
    res = Fido.search(instrument & time & level)
    files = Fido.fetch(res[0, 0], path=tmp_path / "{instrument}")
    assert len(files) == 1
    for f in files:
        assert "EUI" in f


def test_registered_attrs():
    attr_str = str(a.soar.Product)
    # Check that at least one attr value is registered
    assert "epd_ept_asun_burst_ele_close" in attr_str


def test_registered_instr_attrs():
    # Check if the Solo instruments are registered in a.Instrument
    instr_attr = a.Instrument
    assert "SOAR" in instr_attr._attr_registry[instr_attr].client
    assert "stix" in instr_attr._attr_registry[instr_attr].name


def test_registered_soop_names():
    # Check if the soop names are registered in a.soar.SOOP
    soop_attr = str(a.soar.SOOP)
    assert "\nr_small_mres_mcad_ar_long_term" in soop_attr


def test_search_soop():
    instrument = a.Instrument("EUI")
    time = a.Time("2022-04-01 01:00", "2022-04-01 02:00")
    soop_attr = a.soar.SOOP.r_small_mres_mcad_ar_long_term
    res = Fido.search(time, instrument, soop_attr)
    assert "SOOP Name" in res[0].columns
    assert res.file_num == 16

    # test non valid soop name passed
    res = Fido.search(time, instrument, a.soar.SOOP("hello"))
    assert res.file_num == 0


def test_when_soar_provider_passed():
    # Tests when a.Provider.soar is passed that only SOARClient results are returned
    instrument = a.Instrument("EUI")
    time = a.Time("2022-04-01 00:00", "2022-04-01 01:00")
    provider = a.Provider.soar
    res = Fido.search(time & instrument & provider)
    assert len(res) == 1
    assert "soar" in res.keys()


def test_when_sdac_provider_passed():
    # tests that only VSO EUI results are returned when explicitly setting the provider to SDAC
    instrument = a.Instrument("EUI")
    time = a.Time("2022-04-01 00:00", "2022-04-01 01:00")
    provider = a.Provider.sdac
    res = Fido.search(time & instrument & provider)
    assert len(res) == 1
    assert "vso" in res.keys()


def test_when_wrong_provider_passed():
    # Tests that no results are returned when a provider is passed which does not provide EUI data.
    # This is different from the above test because the SDAC and the SOAR both provide EUI data while
    # NOAA has no overlap with the data provided by the SOAR.
    instrument = a.Instrument("EUI")
    time = a.Time("2022-04-01 00:00", "2022-04-01 01:00")
    provider = a.Provider.noaa
    res = Fido.search(time & instrument & provider)
    assert len(res) == 0
