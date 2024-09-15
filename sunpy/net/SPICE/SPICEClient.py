from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.SPICE import attrs as a
from sunpy.net.SPICE.PSP import attrs as ps
from sunpy.net.SPICE.PSP.psp import PSPKernel
from sunpy.net.SPICE.Solo import attrs as sa
from sunpy.net.SPICE.Solo.solar_orbiter import SoloKernel

__all__ = ['SPICEClient']

class SPICEClient(BaseClient):
    """
    A centralized SPICE client to handle both PSP and Solo kernel data.
    """

    kernel_classes = {
        'PSP': PSPKernel,
        'Solo': SoloKernel
    }

    def search(self, *query, missions=None):
        """
        Search for SPICE kernels based on mission and other criteria.

        Parameters:
        - query: search parameters for kernel type, time range, instrument, etc.
        - missions: list of missions to search for, defaults to ['PSP', 'Solo'] if not specified.
        """
        if missions is None:
            missions = ('PSP', 'Solo')

        results = []
        query_params = {}
        kernel_type = None

        # Extract query parameters

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

        # Combine results into a single QueryResponseTable
        return QueryResponseTable(results, client=self)

    def fetch(self, query_results, path=None, **kwargs):
        """
        Fetch the selected kernels from both PSP and Solo if required.
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
        """Check if this client can handle the given query."""
        return any(isinstance(q, a.Kernel_type) for q in query)
