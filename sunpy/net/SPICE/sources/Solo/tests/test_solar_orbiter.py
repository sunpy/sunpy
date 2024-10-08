import pytest

from sunpy.net.SPICE.Solo.solar_orbiter import SoloKernel


@pytest.fixture
def Solo():
    return SoloKernel("lsk")

@pytest.mark.remote_data
def test_all_links(Solo):
    links = Solo.get_all_links()
    assert isinstance(links,list)
    assert "naif0012.tls" in links

@pytest.mark.remote_data
def test_get_readme(Solo):
    readme = Solo.get_readme()
    assert readme == "aareadme.txt"

@pytest.mark.remote_data
def test_filter_kernels_readme(Solo):
    readme = Solo.filter_kernels(get_readme = True)
    assert isinstance(readme,dict)
    assert {0:"aareadme.txt"} == readme

@pytest.mark.remote_data
def test_filter_kernels_index(Solo):
    query1 = {
        "index":0
    }
    assert {0:"aareadme.txt"} == Solo.filter_kernels(query1)
