from abc import ABC, abstractmethod


class BaseClient(ABC):
    """
    This defines the Abstract Base Class for each download client.

    The BaseClient has several abstract methods that ensure that any subclass enforces the bare minimum API.
    These are `search`, `fetch` and `_can_handle_query`.
    The last one ensures that each download client can be registered with Fido.

    All download clients should subclass `~sunpy.net.dataretriever.GenericClient`.
    If the structure of `~sunpy.net.dataretriever.GenericClient`
    is not useful you should use `~sunpy.net.BaseClient`.
    `~sunpy.net.vso.VSOClient` and `~sunpy.net.jsoc.JSOCClient`
    are examples of download clients that subclass `BaseClient`.
    """

    _registry = dict()

    def __init_subclass__(cls, *args, **kwargs):
        """
        An __init_subclass__ hook initializes all of the subclasses of a given class.
        So for each subclass, it will call this block of code on import.
        This replicates some metaclass magic without the need to be aware of metaclasses.
        Here we use this to register each subclass in a dict that has the `_can_handle_query` attribute.
        This is then passed into the UnifiedDownloaderFactory so we can register them.
        This means that Fido can use the clients internally.
        """
        super().__init_subclass__(**kwargs)
        # We do not want to register GenericClient since its a dummy client.
        if cls.__name__ is not 'GenericClient':
            cls._registry[cls] = cls._can_handle_query

    @abstractmethod
    def search(self, *args, **kwargs):
        """
        This enables the user to search for data using the client.
        """

    @abstractmethod
    def fetch(self, *args, **kwargs):
        """
        This enables the user to fetch the data using the client, after a search.
        """

    @classmethod
    @abstractmethod
    def _can_handle_query(cls, *query):
        """
        This enables the client to register what kind of searches it can handle, to prevent Fido using the incorrect client.
        """
