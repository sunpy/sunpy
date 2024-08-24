import os
import urllib.request
from pathlib import Path

from bs4 import BeautifulSoup

from astropy.time import Time

from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.SPICE.Solo import attrs as sa
from sunpy.util.parfive_helpers import Downloader, Results

BASE_URL = "https://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/{}"

class SoloKernel:
    def __init__(self, kernel_type):
        self.kernel_urls = BASE_URL.format(kernel_type)
        self.kernel_type = kernel_type



    def get_all_links(self):
        """
        see all the available kernels for specific kernel type
        """
        try:
            start = False
            links =[]
            with urllib.request.urlopen(self.kernel_urls) as response:
                soup = BeautifulSoup(response, "html.parser")
                for link in soup.find_all("a",href = True):
                  if link.get_text() =="aareadme.txt":
                      start = True

                  if link.get_text().endswith("/"):
                      continue

                  if start:
                      links.append(link.get_text())

            return links

        except Exception as e:
            print(e)




    def get_readme(self):
        """
        get aareadme.txt for specific kernel
        """
        readme = self.get_link_by_index()[0]
        return readme

    def get_link_by_index(self):
        """
        returns a dictionary with each link mapped with its index
        """
        mapped = {}
        for index , link in enumerate(self.get_all_links()):
            mapped[index] = link
        return mapped

    def download_by_index(self,ind,overwrite = False,progress = True,wait = True,downloader = None, path = None,):
        """
        Allows downloading links with their corresponding index.
        Args being index.
        """
        links_mapped = self.filter_kernels()


        # Create a downloader instance
        downloader = Downloader(progress = progress,overwrite= overwrite)
        index = ind.tolist()
        file_url = self.kernel_urls + "/" + links_mapped[index]

        if path is None:
            default_dir = "Downloads" + f"_{self.kernel_type}"
            path = os.path.join(default_dir,'{file}')

        elif isinstance(path,Path):
            path = str(path)
        if isinstance(path,str) and '{file}' not in path:
            path = os.path.join(path,'{file}')

        if not index:
            file_name = path.format(file = f"for {self.kernel_type} {links_mapped[index]}")
        else:
            file_name = path.format(file = links_mapped[index])
        downloader.enqueue_file(file_url, path=file_name, filename=file_name)
        print(f"Queued for download: {index} -- {file_name}")
        print(file_name)



        # Start downloading files
        if not wait:
            return Results()

        results =  downloader.download()
        return results


    def filter_kernels(self,get_readme = False,**kwargs):
        """
        Filter kernels based on search terms using specified boolean logic.
        Returns:
        - A dictionary with the filtered kernels retaining the original indices.
        """
        filtered_kernel = {}
        original_links = self.get_link_by_index()

        if get_readme:
            filtered_kernel[0] = self.get_readme()
            return filtered_kernel
        if 'start' in kwargs:
            kwargs['start'] = Time(kwargs['start']).tdb.strftime('%Y%m%d')
            print(kwargs['start'])
        if 'end' in kwargs and kwargs["end"] is not None:
            kwargs['end'] = Time(kwargs['end']).tdb.strftime('%Y%m%d')
            print(kwargs["end"])
        if "version" in kwargs:
            kwargs["version"] = "V" + str(kwargs["version"])
        for index, link in original_links.items():
            match = None


            match = all(value in link for value in kwargs.values())

            if match:
                filtered_kernel[index] = link
        if not len(filtered_kernel):
            print("no match found for the search terms!")

        return filtered_kernel

class SoloResponseTable(QueryResponseTable):
    """
    A table for storing spice kerenels
    """

class SoloClient(BaseClient):
    """
    this is for solo client (will add more in future )
    """
    @property
    def info_url(self):
        return "https://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/"

    def search(self, *query):
        """
        Search for SPICE kernels based on mission and other criteria.
        """
        results = []
        query_params = {}
        kernel_type = None

        for q in query:
            if isinstance(q, sa.Kernel_type):
                kernel_type = q.value
            if isinstance(q, sa.Time):
                query_params['start'] = q.start
                query_params['end'] = q.end
            if isinstance(q, sa.Instrument):
                query_params['instrument'] = q.value
            if isinstance(q, sa.link):
                query_params['link'] = q.value
            if isinstance(q,sa.sensor):
                query_params["sensor"] = q.value
            if isinstance(q,sa.Version):
                query_params["version"] = q.value
            if isinstance(q,sa.Readme):
                if q.value:
                    query_params["get_readme"] = True
                else:
                    continue


        if not kernel_type:
            raise ValueError("Kernel type must be specified in the query.")

        solo_kernel = SoloKernel(kernel_type)
        filtered_kernels = solo_kernel.filter_kernels(**query_params)

        for index, link in filtered_kernels.items():
            results.append({
                'Mission': "solo",
                'Kernel': kernel_type,
                'Link': link,
                'Index': index
            })

        return SoloResponseTable(results,client=self)

    def fetch(self, query_results,path=None, **kwargs):
        """
        Fetch the selected kernels.
        """
        for result in query_results:
            kernel_type = result['Kernel']
            index = result['Index']


            solo_kernel = SoloKernel(kernel_type)
            solo_kernel.download_by_index(index,overwrite = False,progress = True,wait = True,path=path)

    @staticmethod
    def _can_handle_query(*query):
        """Check if this client can handle the given query."""
        return any(isinstance(q, sa.Kernel_type) for q in query)
