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
        client.search(a.Time("2024-01-01"))
