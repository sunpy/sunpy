import pytest

from sunpy.net import vso


@pytest.fixture
def client():
    return vso.VSOClient()
