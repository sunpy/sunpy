import os
import json

from sunpy.net import attrs as a
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.SPICE.attrs import walker
# import sunpy.net.SPICE.attrs as a
from sunpy.net.SPICE.sources.PSP.psp import PSPKernel
from sunpy.net.SPICE.sources.solar_orbiter import SoloKernel

__all__ = ['SPICEClient']

class SPICEClient(BaseClient):
    """
    Provides access to  SPICE data of both PSP and Solar Orbiter (Solo).

    The `SPICEClient` provides methods to search for and download SPICE kernels from both
    the Parker Solar Probe (PSP) and Solar Orbiter (Solo) missions.

    Attributes
    ----------
    kernel_classes : dict
        A dictionary mapping mission names ('PSP', 'Solo') to their respective kernel classes.

    Methods
    -------
    search(*query, missions=None):
        Search for SPICE kernels based on mission and other query criteria.
    fetch(query_results, path=None, **kwargs):
        Download selected SPICE kernels based on search results.
    _can_handle_query(*query):
        Check if this client can handle the given query.

    Examples
    --------
    >>> from sunpy.net.SPICE import attrs as a
    >>> from sunpy.net.SPICE.SPICEClient import SPICEClient
    >>> from astropy.time import Time
    >>> client = SPICEClient()  # doctest: +REMOTE_DATA
    >>> query = [a.Mission('Solo'), a.Kernel_type('lsk')]   # doctest: +REMOTE_DATA
    >>> results = client.search(*query) # doctest: +REMOTE_DATA
    >>> print(results) # doctest: +REMOTE_DATA
    Mission Kernel     Link     Index
    ------- ------ ------------ -----
       Solo    lsk aareadme.txt     0
       Solo    lsk naif0012.tls     1

    """
    kernel_classes = {
        'psp': PSPKernel,
        'solo': SoloKernel
    }

    def search(self, *query):
        """
        Search for SPICE kernels based on mission and other criteria.

        Parameters
        ----------
        *query : `sunpy.net.attr`
            Search parameters such as kernel type, time range, instrument, etc.
        missions : tuple, optional
            A list of missions to search for. Defaults to both 'PSP' and 'Solo'.

        Returns
        -------
        `sunpy.net.base_client.QueryResponseTable`
            A table containing the search results.

        Raises
        ------
        ValueError
            If no kernel type is specified in the query, or if an unsupported mission is provided.

        Examples
        --------
        >>> from sunpy.net.SPICE import attrs as a
        >>> from sunpy.net.SPICE.SPICEClient import SPICEClient
        >>> from astropy.time import Time
        >>> client = SPICEClient()  # doctest: +REMOTE_DATA
        >>> query = [a.Mission('PSP'), a.Kernel_type('lsk')]    # doctest: +REMOTE_DATA
        >>> results = client.search(*query) # doctest: +REMOTE_DATA
        >>> print(results)  # doctest: +REMOTE_DATA
        Mission Kernel     Link     Index
        ------- ------ ------------ -----
            PSP    lsk naif0012.tls     0
        """

        missions=None
        results = []
        kernel_type = None
        query = and_(*query)
        block  = walker.create(query)[0]

        if "kernel_type" in block:
            kernel_type = block["kernel_type"]
            block.pop("kernel_type")
        else:
            raise ValueError("Kernel type must be specified in the query.")

        if "obsevotory" in block:
            missions = block["obsevotory"]
            block.pop("obsevotory")
        else:
            missions = ('psp', 'solo')

        psp_recognized = ["fk","lsk","sclk","pck","ahk","pek","ltapk","ik"]
        solo_recognized = ["ck", "fk", "ik", "lsk", "pck", "sclk", "spk","mk"]

        if "psp" in missions and kernel_type not in psp_recognized:
            missions = ("solo",)
        elif "solo" in missions and kernel_type not in solo_recognized:
            missions = ("psp",)

        if any(param in block for param in ["get_readme", "voem", "sensor"]):
            missions = ("solo",)
        elif any(param in block for param in ["Analysis_fk", "numupdates"]):
            missions = ("psp",)

        for mission in missions:
            if mission not in self.kernel_classes:
                raise ValueError(f"Unsupported mission: {mission}. Supported missions: {list(self.kernel_classes.keys())}")
            kernel_class = self.kernel_classes[mission](kernel_type)
            filtered_kernels = kernel_class.filter_kernels(**block)

            for index, link in filtered_kernels.items():
                results.append({
                    'Mission': mission,
                    'Kernel': kernel_type,
                    'Link': link,
                    'Index': index
                })

        return QueryResponseTable(results, client=self)

    @classmethod
    def _attrs_module(cls):
        return 'SPICE', 'sunpy.net.SPICE.attrs'

    def fetch(self, query_results, path=None, **kwargs):
        """
        Fetch the selected kernels from both PSP and Solo missions based on the search results.

        Parameters
        ----------
        query_results : list
            A list of query result entries returned by the `search` method.
        path : str, optional
            The directory path where the kernels will be downloaded. Defaults to the current directory.
        **kwargs : dict
            Additional download options.
        """
        for result in query_results:
            mission = result['Mission']
            kernel_type = result['Kernel']
            index = result['Index']
            kernel_class = self.kernel_classes[mission](kernel_type)
            kernel_class.download_by_index(index, overwrite=False, progress=True, wait=True, path=path)

    @staticmethod
    def load_spice_values():
        from sunpy.net import attrs as a
        here = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(here, 'data', 'attrs.json')) as attrs_file:
            keyword_info = json.load(attrs_file)
        keyword_info = list(keyword_info.items())

        attrs = {a.SPICE.Obsevotory:keyword_info}
        return attrs

    @classmethod
    def register_values(cls):
        """
        Register the SOAR specific attributes with Fido.

        Returns
        -------
        dict
            The dictionary containing the values formed into attributes.
        """
        return cls.load_spice_values()

    @staticmethod
    def _can_handle_query(*query):
        """
        Check if this client can handle the given query.

        Parameters
        ----------
        *query : `sunpy.net.attr`
            The query parameters provided to the client.

        Returns
        -------
        bool
        """
        return any(isinstance(q, a.SPICE.Kernel_type) for q in query)
