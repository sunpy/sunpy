from pathlib import Path

import astropy.units as u
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

from sunpy.net.soar.client import SOARClient


def test_search():
    instrument = a.Instrument("EUI")
    time = a.Time("2022-02-11", "2022-02-12")
    level = a.Level(1)
    product = a.soar.Product("eui-fsi174-image")

    res = Fido.search(instrument, time, level, product)
    assert len(res) == 1
    assert len(res[0]) == 660
    assert u.allclose(res[0, 0]["Filesize"], 2.439 * u.Mbyte)

    # check passing upper case descriptor
    product = a.soar.Product("EUI-FSI174-IMAGE")
    res = Fido.search(time, level, product)
    assert res.file_num == 660

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


def test_insitu_search():
    instrument = a.Instrument("MAG")
    time = a.Time("2020-04-16", "2020-04-17")
    product = a.soar.Product("mag-rtn-normal-1-minute")

    res = Fido.search(instrument, time, product)
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
    assert "SOAR" in instr_attr._attr_registry[instr_attr].client  # NOQA: SLF001
    assert "stix" in instr_attr._attr_registry[instr_attr].name  # NOQA: SLF001


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
    assert res["soar"]


def test_when_sdac_provider_passed():
    # tests that only VSO EUI results are returned when explicitly setting the provider to SDAC
    instrument = a.Instrument("EUI")
    time = a.Time("2022-04-01 00:00", "2022-04-01 01:00")
    provider = a.Provider.sdac
    res = Fido.search(time & instrument & provider)
    assert len(res) == 1
    assert res["vso"]


def test_when_wrong_provider_passed():
    # Tests that no results are returned when a provider is passed which does not provide EUI data.
    # This is different from the above test because the SDAC and the SOAR both provide EUI data while
    # NOAA has no overlap with the data provided by the SOAR.
    instrument = a.Instrument("EUI")
    time = a.Time("2022-04-01 00:00", "2022-04-01 01:00")
    provider = a.Provider.noaa
    res = Fido.search(time & instrument & provider)
    assert len(res) == 0


def test_search_wavelength_detector_column():
    instrument = a.Instrument("EUI")
    time = a.Time("2021-02-01", "2021-02-02")
    level = a.Level(1)
    product = a.soar.Product("EUI-FSI174-IMAGE")
    res = Fido.search(instrument & time & level & product)
    assert "Wavelength" in res[0].columns
    assert "Detector" in res[0].columns


def test_search_detector_instrument_dimension_2():
    # Instruments "EUI", "METIS", "PHI" and "SOLOHI" have two dimensions in the SOAR data.
    # Selecting no dimension index in the query results in two identical output rows.
    # To avoid repeating data, we have methods to take dimension index=1, which avoids any repetition.
    instrument = a.Instrument("EUI")
    time = a.Time("2020-03-03", "2020-03-04")
    level = a.Level(1)
    detector = a.Detector("HRI_EUV")
    res = Fido.search(instrument & time & level & detector)
    assert "Detector" in res[0].columns
    assert res.file_num == 266


def test_search_detector_instrument_dimension_4():
    # The "SPICE" instrument has four dimensions in the SOAR data. As a result,
    # selecting no dimension index in the query results in four identical output rows.
    # To avoid repeating data, we have methods to take dimension index=1, which avoids any repetition.
    instrument = a.Instrument("SPICE")
    time = a.Time("2023-03-03 15:00", "2023-03-03 16:00")
    level = a.Level(1)
    detector = a.Detector("SW")
    res = Fido.search(instrument & time & level & detector)
    assert "Detector" in res[0].columns
    assert res.file_num == 11


def test_invalid_detector():
    instrument = a.Instrument("SPICE")
    time = a.Time("2023-03-03 15:00", "2023-03-03 16:00")
    level = a.Level(1)
    detector = a.Detector("hello")
    res = Fido.search(instrument & time & level & detector)
    assert "Detector" in res[0].columns
    assert res.file_num == 0


def test_wavelength_column_wavelength_exists():
    # For instruments EUI, METIS and SOLOHI "wavelength" column is available.
    # Test to check if the "Wavelength" column exists in the search results.
    instrument = a.Instrument("EUI")
    time = a.Time("2023-04-03 15:00", "2023-04-03 16:00")
    level = a.Level(1)
    wavelength = a.Wavelength(304 * u.AA)
    res = Fido.search(instrument & time & level & wavelength)
    assert "Wavelength" in res[0].columns
    assert res.file_num == 12


def test_wavelength_single():
    # Test to check if the wavelength value is filtered for a single value provided.
    instrument = a.Instrument("EUI")
    time = a.Time("2023-04-03 15:00", "2023-04-03 16:00")
    level = a.Level(1)
    wavelength = a.Wavelength(304 * u.AA)
    res = Fido.search(instrument & time & level & wavelength)
    for table in res:
        assert all(table["Wavelength"] == 304)


def test_wavelength_range():
    # Test to check if the wavelength value is filtered for wavemin and wavemax provided.
    instrument = a.Instrument("EUI")
    time = a.Time("2023-04-03 15:00", "2023-04-03 16:00")
    level = a.Level(1)
    wavelength = a.Wavelength(171 * u.AA, 185 * u.AA)
    res = Fido.search(instrument & time & level & wavelength)
    for table in res:
        assert all(table["Wavelength"] == 174)


def test_join_science_query():
    result = SOARClient._construct_payload(  # NOQA: SLF001
        [
            "instrument='EUI'",
            "begin_time>='2021-02-01+00:00:00'+AND+begin_time<='2021-02-02+00:00:00'",
            "level='L1'",
            "descriptor='eui-fsi174-image'",
        ]
    )

    assert result["QUERY"] == (
        "SELECT+h1.instrument, h1.descriptor, h1.level, h1.begin_time, h1.end_time, "
        "h1.data_item_id, h1.filesize, h1.filename, h1.soop_name, h2.detector, h2.wavelength, "
        "h2.dimension_index+FROM+v_sc_data_item AS h1 JOIN v_eui_sc_fits AS h2 USING (data_item_oid)"
        "+WHERE+h1.instrument='EUI'+AND+h1.begin_time>='2021-02-01+00:00:00'+AND+h1.begin_time<='2021-02-02+00:00:00'"
        "+AND+h2.dimension_index='1'+AND+h1.level='L1'+AND+h1.descriptor='eui-fsi174-image'"
    )


def test_join_low_latency_query():
    result = SOARClient._construct_payload(  # NOQA: SLF001
        [
            "instrument='EUI'",
            "begin_time>='2021-02-01+00:00:00'+AND+begin_time<='2021-02-02+00:00:00'",
            "level='LL01'",
            "descriptor='eui-fsi174-image'",
        ]
    )

    assert result["QUERY"] == (
        "SELECT+h1.instrument, h1.descriptor, h1.level, h1.begin_time, h1.end_time, "
        "h1.data_item_id, h1.filesize, h1.filename, h1.soop_name, h2.detector, h2.wavelength, "
        "h2.dimension_index+FROM+v_ll_data_item AS h1 JOIN v_eui_ll_fits AS h2 USING (data_item_oid)"
        "+WHERE+h1.instrument='EUI'+AND+h1.begin_time>='2021-02-01+00:00:00'+AND+h1.begin_time<='2021-02-02+00:00:00'"
        "+AND+h2.dimension_index='1'+AND+h1.level='LL01'+AND+h1.descriptor='eui-fsi174-image'"
    )
