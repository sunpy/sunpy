# -*- coding: utf-8 -*-

from sunpy.net import base_client


def test_registry():
    """
    Just to see that some of the clients are in the registry.
    """

    keys = base_client.BaseClient._registry.keys()
    items = base_client.BaseClient._registry.items()

    # Todo: Make this better
    assert len(keys) == 10
    assert len(items) == 10
