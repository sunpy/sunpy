import re

import numpy as np
import pytest

from astropy.time import Time

from sunpy.net.base_client import QueryResponseTable
from sunpy.net.SPICE import attrs as a
from sunpy.net.SPICE.SPICEClient import SPICEClient


@pytest.fixture
def client():
    return SPICEClient()

@pytest.mark.remote_data
def test_search(client):
    query = [
        a.Mission('Solo'),
        a.Kernel_type('lsk'),
    ]
    results = client.search(*query)

    expected_response = QueryResponseTable([
        {"Mission": "Solo", "Kernel": "lsk", "Link": "aareadme.txt", "Index": np.array(0)},
        {"Mission": "Solo", "Kernel": "lsk", "Link": "naif0012.tls", "Index": np.array(1)}
    ])

    assert len(query) == len(expected_response)
    for q, expected in zip(results, expected_response):
        assert q["Mission"] == expected["Mission"]
        assert q["Kernel"] == expected["Kernel"]
        assert q["Link"] == expected["Link"]
        assert np.array_equal(q["Index"], expected["Index"])

@pytest.mark.parametrize(("kernel_type", "expected"), [
    (a.Kernel_type("lsk"), True),
    (a.Kernel_type("spk"), True),
    (a.Time(Time("2020-01-01"), Time("2020-12-31")), False)
])
def test_can_handle_query(client, kernel_type, expected):
    can_handle = client._can_handle_query(kernel_type)
    assert can_handle == expected

@pytest.mark.remote_data
def test_no_url_found(client):
    data = client.search(a.Kernel_type("lsk"), a.Time("1979-01-01"))
    assert len(data) == 0

def test_no_kernel(client):
    with pytest.raises(ValueError, match="Kernel type must be specified in the query."):
        client.search(a.Time("2024-01-01","2025-01-01"))

@pytest.mark.remote_data
def test_response_type(client):
    query = [
        a.Kernel_type("ik"),
        a.Instrument("eui"),
        a.Version("01"),
        a.Readme(True),
    ]
    response = client.search(*query)
    assert isinstance(response,QueryResponseTable)

@pytest.mark.remote_data
def test_direct_link(client):
    url = client.search(a.Kernel_type("ck"), a.Link("solo_ANC_soc-default-att-stp_20200210-20301120_272_V1_00276_V01.bc"),a.Version("01"))
    assert len(url) == 1
    assert "solo_ANC_soc-default-att-stp_20200210-20301120_272_V1_00276_V01.bc" in str(url[0]["Link"])

def test_wrong_mission(client):
    query = [
        a.Kernel_type("ik"),
        a.Mission("OHNO"),
    ]
    with pytest.raises(ValueError , match =re.escape("Unsupported mission: OHNO. Supported missions: ['PSP', 'Solo']")):
        client.search(*query)

@pytest.mark.remote_data
def test_fetch(client,tmp_path):
    query = client.search(a.Kernel_type("pck"))

    client.fetch(query,path = tmp_path)
    downloaded_files = list(tmp_path.glob("*"))
    assert len(downloaded_files) > 0
    assert any("pck" in str(file) for file in downloaded_files)

@pytest.mark.remote_data
def test_index(client):
    query = [a.Kernel_type("ik"),a.Index(1)]
    response = client.search(*query)
    expected_response = QueryResponseTable([
        {"Mission": "PSP", "Kernel": "ik", "Link": "spp_epilo_v100.ti", "Index": np.array(1)},
        {"Mission": "Solo","Kernel": "ik","Link": "solo_ANC_soc-epd-ik_V00.ti","Index": np.array(1)}
    ])
    for q, expected in zip(response, expected_response):
        assert q["Mission"] == expected["Mission"]
        assert q["Kernel"] == expected["Kernel"]
        assert q["Link"] == expected["Link"]
        assert np.array_equal(q["Index"], expected["Index"])
    assert isinstance(response,QueryResponseTable)

@pytest.mark.remote_data
def test_get_readme(client):
    query = [a.Kernel_type("ik"),a.Readme(True),a.Mission("Solo")]
    response = client.search(*query)

    expected_response = QueryResponseTable([
        {"Mission": "Solo","Kernel": "ik","Link": "aareadme.txt","Index": np.array(0)}
    ])
    for q, expected in zip(response, expected_response):
        assert q["Mission"] == expected["Mission"]
        assert q["Kernel"] == expected["Kernel"]
        assert q["Link"] == expected["Link"]
        assert np.array_equal(q["Index"], expected["Index"])
    assert isinstance(response,QueryResponseTable)

@pytest.mark.remote_data
def test_get_Analysis_fk(client):
    query = [
        a.Kernel_type("fk"),
        a.Analysis_fk(True),
        a.Version("200")
    ]
    response = client.search(*query)
    expected_response = QueryResponseTable([
        {"Mission": "PSP","Kernel": "fk","Link": "spp_dyn_v200.tf","Index": np.array(4)}
    ])
    for q, expected in zip(response, expected_response):
        assert q["Mission"] == expected["Mission"]
        assert q["Kernel"] == expected["Kernel"]
        assert q["Link"] == expected["Link"]
        assert np.array_equal(q["Index"], expected["Index"])
    assert isinstance(response,QueryResponseTable)

@pytest.mark.remote_data
def test_get_pe(client):
    query = [
        a.Kernel_type("pek"),
    ]

    response = client.search(*query)
    expected_response = QueryResponseTable([
        {"Mission": "PSP","Kernel": "pek","Link": "de430.bsp","Index": np.array(0)}
    ])
    for q, expected in zip(response, expected_response):
        assert q["Mission"] == expected["Mission"]
        assert q["Kernel"] == expected["Kernel"]
        assert q["Link"] == expected["Link"]
        assert np.array_equal(q["Index"], expected["Index"])
    assert isinstance(response,QueryResponseTable)
