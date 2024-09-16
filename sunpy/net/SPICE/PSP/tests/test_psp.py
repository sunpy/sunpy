import pytest

from sunpy.net.SPICE.PSP.psp import PSPKernel


@pytest.fixture
def PSP():
    return PSPKernel("lsk")

@pytest.mark.remote_data
def test_all_links(PSP):
    link = PSP.get_all_links()
    assert isinstance(link,list)
    assert "naif0012.tls" in link

@pytest.mark.remote_data
def test_links_by_index(PSP):
    index = PSP.get_link_by_index()
    assert isinstance(index,dict)
    assert {0:"naif0012.tls"} == index

@pytest.mark.remote_data
def test_filter_kernels(PSP):
    query = {
        "version": "100",
    }
    filtered = PSP.filter_kernels(query)
    assert isinstance(filtered,dict)
    assert len(filtered) == 0
