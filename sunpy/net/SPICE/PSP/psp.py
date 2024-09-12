import os
import urllib.request
from pathlib import Path

from bs4 import BeautifulSoup

from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.util.parfive_helpers import Downloader, Results
from sunpy.net.SPICE.PSP import attrs as ps
from astropy.time import Time

cheat_sheet = {
    "fk":"PSP_Frame_Kernels/",
    "lsk":"Leap_Second_Kernel/",
    "sclk":"SCLK_files/",
    "reck":"Reconstructed_Ephemerides/",
    "pck":"Planetary_Constant_Kernel/",
    "ahk":"Attitude_History_Kernels/",
    "stapk":"Short_Term_Attitude_Predict_Kernels/",
    "pek":"Planetary_Ephemerides/",
    "ltapk":"Long_Term_Attitude_Predict_Kernels/",
    "ltpek":"Long_Term_Predicted_Ephemeris/",
    "ik":"PSP_Frame_Kernels/"
}

BASE_URL = "https://spdf.gsfc.nasa.gov/pub/data/psp/ephemeris/spice/{}"

class PSPKernel:
    def __init__(self, kernel_type):
        self.kernel_type = cheat_sheet[kernel_type]
        self.kernel_urls = BASE_URL.format(self.kernel_type)
        print(self.kernel_urls)

    def get_all_links(self):
        try:
            links = []
            with urllib.request.urlopen(self.kernel_urls) as response:
                soup = BeautifulSoup(response, "html.parser")
                for link in soup.find_all("a", href=True):
                    href = link.get("href")
                    link_text = link.get_text()
                    if href and not href.endswith("/") and all(x not in link_text for x in ["Name", "Last modified", "Size"]):
                        if self.kernel_type == "ik" and link_text.endswith("ti"):
                            links.append(link_text)
                        else:
                            links.append(link_text)

            return links
        except Exception as e:
            print(f"An error occurred while retrieving links: {e}")


    def get_link_by_index(self):
        mapped = {}

        for index , link in enumerate(self.get_all_links()):
            mapped[index] = link
        return mapped

    def download_by_index(self,index,overwrite = False,progress = True,wait = True,downloader = None, path = None,):
        """
        Allows downloading links with their corresponding index.
        ind being index.
        """
        links_mapped = self.filter_kernels()


        downloader = Downloader(progress = progress,overwrite= overwrite)
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



        if not wait:
            return Results()

        results =  downloader.download()
        return results

    def filter_kernels(self,**kwargs):


        filtered_kernel = {}
        original_links = self.get_link_by_index()

        if kwargs["Analysis"]:
            kwargs["Analysis"] = "spp_dyn"
        if "index" in kwargs:
            for i,j in enumerate(kwargs["index"]):
                filtered_kernel[j] = original_links[j]
            return filtered_kernel
        if 'start' in kwargs:
            kwargs['start'] = Time(kwargs['start']).tdb.strftime('%Y%m%d')
        if 'end' in kwargs and kwargs["end"] is not None:
            kwargs['end'] = Time(kwargs['end']).tdb.strftime('%Y%m%d')
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


class PSPResponseTable(QueryResponseTable):
    """
    A table for storing psp spice kernels
    """
 
class PSPClient(BaseClient):
    """
    PASS
    """
    def search(self, *query):
        """
        Search for SPICE kernels based on mission and other criteria.
        """
        results = []
        query_params = {}
        kernel_type = None

        for q in query:
            if isinstance(q, ps.Kernel_type):
                kernel_type = q.value
            if isinstance(q, ps.Time):
                query_params['start'] = q.start
                query_params['end'] = q.end
            if isinstance(q, ps.Instrument):
                query_params['instrument'] = q.value
            if isinstance(q, ps.Link):
                query_params['link'] = q.value
            if isinstance(q,ps.Version):
                query_params["version"] = q.value
            if isinstance(q,ps.Numupdates):
                query_params["Numupdates"] = q.value
            if isinstance(q,ps.Index):
                query_params["index"] = q.value
            if isinstance(q,ps.Analysis_fk):
                if q.value:
                    query_params["Analysis"] = True
                else:
                    continue


        if not kernel_type:
            raise ValueError("Kernel type must be specified in the query.")

        solo_kernel = PSPKernel(kernel_type)
        filtered_kernels = solo_kernel.filter_kernels(**query_params)

        for index, link in filtered_kernels.items():
            results.append({
                'Mission': "PSP",
                'Kernel': kernel_type,
                'Link': link,
                'Index': index
            })

        return PSPResponseTable(results,client=self)

    def fetch(self, query_results,path=None, **kwargs):
        """
        Fetch the selected kernels.
        """
        for result in query_results:
            kernel_type = result['Kernel']
            index = result['Index']


            solo_kernel = PSPKernel(kernel_type)
            solo_kernel.download_by_index(index,overwrite = False,progress = True,wait = True,path=path)

    @staticmethod
    def _can_handle_query(*query):
        """Check if this client can handle the given query."""
        return any(isinstance(q, ps.Kernel_type) for q in query)
