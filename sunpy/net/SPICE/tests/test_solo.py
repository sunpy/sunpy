import numpy as np
import pytest

from sunpy.net.SPICE.Solo.solar_orbiter import SoloClient, SoloKernel, SoloResponseTable

TEST_DICT = {0: "aareadme.txt", 1: "naif0012.tls"}

@pytest.fixture
def client():
    return SoloClient()

@pytest.fixture
def source():
    return SoloKernel("lsk")

@pytest.fixture
def solo_response_table():
    kernel = SoloResponseTable([
        {"Mission": "solo", "Kernel": "lsk", "Link": "aareadme.txt", "Index": np.array(0)},
        {"Mission": "solo", "Kernel": "lsk", "Link": "naif0012.tls", "Index": np.array(1)}
    ])
    return kernel

@pytest.mark.remote_data
def test_get_all_links(source):
    links = source.get_all_links()
    assert isinstance(links, list)
    assert links == ["aareadme.txt", "naif0012.tls"]

@pytest.mark.remote_data
def test_get_readme(source):
    readme = source.get_readme()
    assert readme == "aareadme.txt"

@pytest.mark.remote_data
def test_get_link_by_index(source):
    assert source.get_link_by_index() == TEST_DICT

@pytest.mark.remote_data
def test_filter_kernels(source):
    readme = source.filter_kernels(get_readme=True)
    assert {0: readme[0]} == {0: "aareadme.txt"}
