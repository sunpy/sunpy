from solar_orbiter import SoloKernel

from sunpy.net.base_client import BaseClient, QueryResponseTable


class SPICEClient(BaseClient):

    def __init__(self, mission,kernel_type):
        self.mission = mission.upper()
        self.kernel_type = kernel_type

    def search(self,**kwargs):
        """
        Search for SPICE kernels based on mission and other criteria.

        Parameters:
        - kernel: Type of kernel to search for (e.g., "ck", "ik", "sclk").
        - instrument: Specific instrument for instrument kernels.
        - version: Version of the kernel.
        - start: Start date (for sclk kernels).
        - end: End date (for sclk kernels).

        Returns:
        - QueryResponseTable with matching kernels.
        """
        results = []

        if self.mission == "SOLO":
            solo_kernel = SoloKernel(self.kernel_type)
            filtered_kernels = solo_kernel.filter_kernels(**kwargs)

            for index, link in filtered_kernels.items():
                results.append({
                    'Mission': self.mission,
                    'Kernel': self.kernel_type,
                    'Link': link,
                    'Index': index
                })

        return QueryResponseTable(results)

    def fetch(self, query_results, path=None, **kwargs):
        """
        Fetch the selected kernels.

        Parameters:
        - query_results: Results returned from the search method.
        - path: Destination path for downloaded files.

        Returns:
        - List of paths to downloaded files.
        """
        if path is None:
            path = './{file}'

        downloaded_files = []

        for result in query_results:
            kernel_type = result['Kernel']

            solo_kernel = SoloKernel(kernel_type)
            index = result['Index']
            downloaded_files.append(solo_kernel.download_by_index(index, destination=path))

        return downloaded_files

    @staticmethod
    def _can_handle_query(*query):
        """
        Check if this client can handle the given query.
        """
        return True
