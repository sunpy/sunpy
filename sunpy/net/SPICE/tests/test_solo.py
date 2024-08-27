import numpy as np
import pytest

from astropy.time import Time

from sunpy.net.SPICE.Solo import attrs as sa
from sunpy.net.SPICE.Solo.solar_orbiter import SoloClient, SoloKernel, SoloResponseTable


@pytest.fixture
def client():
    return SoloClient()

@pytest.fixture
def solo_kernel():
    return SoloKernel("lsk")

@pytest.mark.remote_data
def test_get_all_links(solo_kernel):
    links = solo_kernel.get_all_links()
    assert len(links) > 0
    assert isinstance(links, list)
    assert links == ["aareadme.txt", "naif0012.tls"]

@pytest.mark.remote_data
def test_get_readme(solo_kernel):
    readme = solo_kernel.get_readme()
    assert "aareadme.txt" in readme

@pytest.mark.remote_data
def test_get_link_by_index(solo_kernel):
    links_by_index = solo_kernel.get_link_by_index()
    assert isinstance(links_by_index, dict)
    assert len(links_by_index) > 0

@pytest.mark.remote_data
def test_download_by_index(solo_kernel, tmp_path):
    links_by_index = solo_kernel.get_link_by_index()
    first_index = next(iter(links_by_index.keys()))

    results = solo_kernel.download_by_index(first_index, path=tmp_path)
    assert len(results) > 0

@pytest.mark.remote_data
def test_filter_kernels(solo_kernel):
    readme = solo_kernel.filter_kernels(get_readme=True)
    assert {0: readme[0]} == {0: "aareadme.txt"}

@pytest.mark.remote_data
def test_search(client):
    query = client.search(sa.Kernel_type("lsk"))

    expected_response = SoloResponseTable([
        {"Mission": "solo", "Kernel": "lsk", "Link": "aareadme.txt", "Index": np.array(0)},
        {"Mission": "solo", "Kernel": "lsk", "Link": "naif0012.tls", "Index": np.array(1)}
    ])

    assert len(query) == len(expected_response)
    assert isinstance(query, SoloResponseTable)

    for q, expected in zip(query, expected_response):
        assert q["Mission"] == expected["Mission"]
        assert q["Kernel"] == expected["Kernel"]
        assert q["Link"] == expected["Link"]
        assert np.array_equal(q["Index"], expected["Index"])

@pytest.mark.remote_data
def test_fetch(client, tmp_path):
    query = client.search(sa.Kernel_type("lsk"))
    client.fetch(query[:1], path=tmp_path)

    downloaded_files = list(tmp_path.glob("*"))
    assert len(downloaded_files) > 0
    assert any("lsk" in str(file) for file in downloaded_files)

@pytest.mark.parametrize(("kernel_type", "expected"), [
    (sa.Kernel_type("lsk"), True),
    (sa.Kernel_type("spk"), True),
    (sa.Time(Time("2020-01-01"), Time("2020-12-31")), False)
])
def test_can_handle_query(client, kernel_type, expected):
    can_handle = client._can_handle_query(kernel_type)
    assert can_handle == expected

@pytest.mark.remote_data
def test_no_url_found(client):
    data = client.search(sa.Kernel_type("lsk"), sa.Time("1979-01-01"))
    assert len(data) == 0

def test_no_kernel(client):
    with pytest.raises(ValueError, match="Kernel type must be specified in the query."):
        client.search(sa.Time("2024-01-01"))

def test_None_kernel(client):
    with pytest.raises(ValueError, match="kernel type is required"):
        client.search(sa.Kernel_type(None))

def test_invalid_kernel(client):
    with pytest.raises(ValueError, match="Kernel type not recognized 'ok'"):
        client.search(sa.Kernel_type("ok"))
