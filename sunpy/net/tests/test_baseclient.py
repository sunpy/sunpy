# -*- coding: utf-8 -*-
import re

import pytest

from sunpy.net import base_client, vso, jsoc, dataretriever

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
