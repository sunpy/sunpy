import os
import urllib.request
from pathlib import Path

from bs4 import BeautifulSoup

from astropy.time import Time

from sunpy.util.parfive_helpers import Downloader, Results

BASE_URL = "https://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/{}"

class SoloKernel:
    def __init__(self, kernel_type):
        self.kernel_urls = BASE_URL.format(kernel_type)
        self.kernel_type = kernel_type



    def get_all_links(self):
        """
        See all the available kernels for a specific kernel type.
        """
        try:
            start = False
            links = []
            with urllib.request.urlopen(self.kernel_urls) as response:
                soup = BeautifulSoup(response, "html.parser")
                for link in soup.find_all("a", href=True):
                    if link.get_text() == "aareadme.txt":
                        start = True

                    if link.get_text().endswith("/"):
                        continue

                    if start:
                        links.append(link.get_text())

            return links


        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    def get_link_by_index(self):
        """
        returns a dictionary with each link mapped with its index
        """
        mapped = {}
        for index , link in enumerate(self.get_all_links()):
            mapped[index] = link
        return mapped

    def download_by_link(self,link,overwrite = False,progress = True,wait = True,downloader = None, path = None,):
            """
            Allows downloading links with their corresponding index.
            ind being index.
            """
            # Create a downloader instance
            downloader = Downloader(progress = progress,overwrite= overwrite,max_splits=1)
            file_url = self.kernel_urls + "/" + link

            if path is None:
                default_dir = "Downloads" + f"_{self.kernel_type}"
                path = os.path.join(default_dir,'{file}')

            elif isinstance(path,Path):
                path = str(path)
            if isinstance(path,str) and '{file}' not in path:
                path = os.path.join(path,'{file}')

            if "aareadme.txt" == link:
                file_name = path.format(file = f"for {self.kernel_type} {link}")
            else:
                file_name = path.format(file = link)


            downloader.enqueue_file(file_url, path=file_name, filename=file_name,max_splits=1)




            if not wait:
                return Results()

            results =  downloader.download()
            return results


    def filter_kernels(self,**kwargs):
        """
        Filter kernels based on search terms using specified boolean logic.
        Returns:
        - A dictionary with the filtered kernels retaining the original indices.
        """
        filtered_kernel = {}
        original_links = self.get_link_by_index()

        if "index" in kwargs:
            for _,j in enumerate(kwargs["index"]):
                filtered_kernel[j] = original_links[j]
            return filtered_kernel
        if 'start' in kwargs:
            kwargs['start'] = Time(kwargs['start']).tdb.strftime('%Y%m%d')
        if 'end' in kwargs and kwargs["end"] is not None:
            kwargs['end'] = Time(kwargs['end']).tdb.strftime('%Y%m%d')


        for index, link in original_links.items():
            match = None


            match = all(value in link for value in kwargs.values())

            if match:
                filtered_kernel[index] = link
        if not len(filtered_kernel):
            print("no match found for the search terms!")

        return filtered_kernel