# -*- coding: utf-8 -*-
from sunpy.net import base_client, vso, jsoc, dataretriever


def test_registry():
    """
    Just to see that some of the clients are in the registry.
    """

    # TODO: Better this
    client_list = [vso.VSOClient, jsoc.JSOCClient, dataretriever.EVEClient,
                   dataretriever.LYRAClient, dataretriever.NOAAIndicesClient,
                   dataretriever.NOAAPredictClient, dataretriever.NoRHClient,
                   dataretriever.RHESSIClient, dataretriever.SRSClient,
                   dataretriever.XRSClient]

    keys = base_client.BaseClient._registry.keys()
    items = base_client.BaseClient._registry.items()

    assert len(keys) == 10
    assert len(items) == 10

    for client in client_list:
        assert client in keys
        assert (client, client._can_handle_query) in items
