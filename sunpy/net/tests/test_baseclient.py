import re

import pytest

from sunpy.net import base_client, dataretriever, jsoc, vso
from sunpy.net.base_client import QueryResponseTable, convert_row_to_table
from sunpy.net.dataretriever.sources.norh import NoRHClient

_REGEX = re.compile(r"Client")

CLIENT_LIST = []

for a_import in [vso, jsoc, dataretriever]:
    for item in dir(a_import):
        if _REGEX.search(item):
            CLIENT_LIST.append(getattr(a_import, item))

CLIENT_LIST.remove(dataretriever.client.GenericClient)

# We can access the registry directly
CLIENT_NAMES = base_client.BaseClient._registry.keys()
CLIENTS_REG = base_client.BaseClient._registry.items()


@pytest.mark.parametrize("client", CLIENT_LIST)
def test_registry(client):
    """
    Check if each client has been registered.
    """
    assert client in CLIENT_NAMES
    assert (client, client._can_handle_query) in CLIENTS_REG


@pytest.fixture
def dummy_response():
    return QueryResponseTable([{'hello': 1}], client=NoRHClient())


def test_slice(dummy_response):
    assert len(dummy_response) == 1

    row = dummy_response[0]
    table = row.as_table()

    assert len(table) == 1
    assert isinstance(table.client, NoRHClient)

    col = dummy_response['hello']
    table = col.as_table()

    assert len(table) == 1
    assert isinstance(table.client, NoRHClient)


def test_path_format_keys(dummy_response):
    assert dummy_response.path_format_keys() == {'hello'}


def test_convert_row_to_table(dummy_response):

    @convert_row_to_table
    def example(self, query_results, **kwargs):
        return query_results

    assert example(None, dummy_response) is dummy_response
    # This is a single row table anyway
    assert example(None, dummy_response[0]) == dummy_response
