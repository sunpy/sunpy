from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.SPICE import attrs as a
from sunpy.net.SPICE.PSP import attrs as ps
from sunpy.net.SPICE.PSP.psp import PSPKernel
from sunpy.net.SPICE.Solo import attrs as sa
from sunpy.net.SPICE.Solo.solar_orbiter import SoloKernel

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
        'PSP': PSPKernel,
        'Solo': SoloKernel
    }

    def search(self, *query, missions=None):
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
        if missions is None:
            missions = ('PSP', 'Solo')

        results = []
        query_params = {}
        kernel_type = None

        for q in query:
            if isinstance(q,a.Mission):
                missions = q.value
            if isinstance(q, a.Kernel_type):
                kernel_type = q.value
            if isinstance(q, a.Time):
                query_params['start'] = q.start
                query_params['end'] = q.end
            if isinstance(q, a.Instrument):
                query_params['instrument'] = q.value
            if isinstance(q, a.Link):
                query_params['link'] = q.value
            if isinstance(q, a.Version):
                query_params['version'] = q.value
            if isinstance(q,sa.Readme) and q.value and missions == ("Solo",):
                query_params["get_readme"] = True
            if isinstance(q,a.Index):
                query_params["index"] = q.value
            if isinstance(q,ps.Analysis_fk) and missions == ("PSP",):
                if q.value:
                    query_params["Analysis_fk"] = True

        if not kernel_type:
            raise ValueError("Kernel type must be specified in the query.")

        # Search for kernels in each specified mission
        psp_recognized = ["fk","lsk","sclk","pck","ahk","pek","ltapk","ik"]
        solo_recognized = ["ck", "fk", "ik", "lsk", "pck", "sclk", "spk","mk"]
        if "PSP" in missions and kernel_type not in psp_recognized:
            missions = ("Solo",)
        elif "Solo" in missions and kernel_type not in solo_recognized:
            missions = ("PSP",)

        for mission in missions:
            if mission not in self.kernel_classes:
                raise ValueError(f"Unsupported mission: {mission}. Supported missions: {list(self.kernel_classes.keys())}")
            kernel_class = self.kernel_classes[mission](kernel_type)
            filtered_kernels = kernel_class.filter_kernels(**query_params)

            for index, link in filtered_kernels.items():
                results.append({
                    'Mission': mission,
                    'Kernel': kernel_type,
                    'Link': link,
                    'Index': index
                })

        return QueryResponseTable(results, client=self)

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

        Raises
        ------
        ValueError
            If the mission is not supported.
        """
        for result in query_results:
            mission = result['Mission']
            kernel_type = result['Kernel']
            index = result['Index']

            if mission not in self.kernel_classes:
                raise ValueError(f"Unsupported mission: {mission}. Supported missions: {list(self.kernel_classes.keys())}")

            kernel_class = self.kernel_classes[mission](kernel_type)
            kernel_class.download_by_index(index, overwrite=False, progress=True, wait=True, path=path)

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
        return any(isinstance(q, a.Kernel_type) for q in query)
