import json
import pathlib

from sunpy.net import attrs as a
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.SPICE.attrs import walker
from sunpy.net.SPICE.sources.psp import PSPKernel
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
    >>> from sunpy.net import attrs as a
    >>> from sunpy.net.SPICE.SPICEClient import SPICEClient
    >>> from astropy.time import Time
    >>> client = SPICEClient()  # doctest: +REMOTE_DATA
    >>> query = [a.SPICE.Observatory.psp,a.Instrument.sweap]   # doctest: +REMOTE_DATA
    >>> results = client.search(*query) # doctest: +REMOTE_DATA
    >>> print(results) # doctest: +REMOTE_DATA
    Mission Kernel        Link       Index
    ------- ------ ----------------- -----
        psp     ik spp_sweap_v100.ti     0
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
    >>> from sunpy.net import attrs as a
    >>> from sunpy.net.SPICE.SPICEClient import SPICEClient
    >>> from astropy.time import Time
    >>> client = SPICEClient()  # doctest: +REMOTE_DATA
    >>> query = [a.SPICE.Observatory.psp,a.Instrument.sweap]   # doctest: +REMOTE_DATA
    >>> results = client.search(*query) # doctest: +REMOTE_DATA
    >>> print(results) # doctest: +REMOTE_DATA
    Mission Kernel        Link       Index
    ------- ------ ----------------- -----
        psp     ik spp_sweap_v100.ti     0
        """
        missions=None
        results = []
        kernel_type = None
        query = and_(*query)
        block  = walker.create(query)[0]

        if "instrument" in block:
            kernel_type = "ik"

        if "observatory" in block:
            missions = block["observatory"]
            block.pop("observatory")

        if missions not in self.kernel_classes:
            raise ValueError(f"Unsupported mission: {missions}. Supported missions: {list(self.kernel_classes.keys())}")
        kernel_class = self.kernel_classes[missions](kernel_type)
        filtered_kernels = kernel_class.filter_kernels(**block)

        for index, link in enumerate(filtered_kernels.values()):
                results.append({
                    'Mission': missions,
                    'Kernel': kernel_type,
                    'Link': link,
                    'Index': index
                })

        return QueryResponseTable(results, client=self)

    @classmethod
    def _attrs_module(cls):
        return 'SPICE', 'sunpy.net.SPICE.attrs'

    def fetch(self, query_results, overwrite = False,path=None, **kwargs):
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
            link = result['Link']
            kernel_class = self.kernel_classes[mission](kernel_type)
            kernel_class.download_by_link(link, overwrite=overwrite, progress=True, wait=True, path=path)

    @staticmethod
    def load_spice_values():
        instr_path = pathlib.Path(__file__).parent / "data" / "instrument_attrs.json"
        with instr_path.open() as instr_attrs_file:
            all_instr = json.load(instr_attrs_file)
        all_instr = list(all_instr.items())
        attrs = {
            a.Instrument:all_instr,
            a.SPICE.Observatory:[('solo', 'solo'), ('psp', 'psp')]
                }
        return attrs

    @classmethod
    def register_values(cls):
        """
        Register the SPICE specific attributes with Fido.

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

        return any((isinstance(q, a.SPICE.Observatory) or isinstance(q,a.SPICE.Link))for q in query)
