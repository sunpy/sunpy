import os
import urllib.request
from pathlib import Path

from bs4 import BeautifulSoup

from astropy.time import Time

from sunpy.util.parfive_helpers import Downloader, Results

cheat_sheet = {
    "fk":"PSP_Frame_Kernels",
    "lsk":"Leap_Second_Kernel",
    "sclk":"SCLK_files",
    "reck":"Reconstructed_Ephemerides",
    "pck":"Planetary_Constant_Kernel",
    "ahk":"Attitude_History_Kernels",
    "stapk":"Short_Term_Attitude_Predict_Kernels",
    "pek":"Planetary_Ephemerides",
    "ltapk":"Long_Term_Attitude_Predict_Kernels",
    "ltpek":"Long_Term_Predicted_Ephemeris",
    "ik":"PSP_Frame_Kernels"
}

BASE_URL = "https://spdf.gsfc.nasa.gov/pub/data/psp/ephemeris/spice/{}"

class PSPKernel:
    def __init__(self, kernel_type):
        self.kernel_type = kernel_type
        self.kernel_urls = BASE_URL.format(cheat_sheet[self.kernel_type])

    def get_all_links(self):
        try:
            links = []
            with urllib.request.urlopen(self.kernel_urls) as response:
                soup = BeautifulSoup(response, "html.parser")
                for link in soup.find_all("a", href=True):
                    href = link.get("href")
                    link_text = link.get_text()

                    # Check for valid href and ensure it's not a directory or unwanted text
                    if href and not href.endswith("/") and all(x not in link_text for x in ["Name", "Last modified", "Size"]):
                        # Special condition for kernel_type "ik"
                        if self.kernel_type == "ik":
                            if link_text.endswith("ti"):
                                links.append(link_text)
                        else:
                            links.append(link_text)  # Append link in all other cases

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


        file_name = path.format(file = links_mapped[index])
        downloader.enqueue_file(file_url, path=file_name, filename=file_name)
        print(f"Queued for download: {index} -- {file_name}")
        print(file_name)



        if not wait:
            return Results()

        results =  downloader.download()
        return results

    def filter_kernels(self,Analysis_fk = False,**kwargs):


        filtered_kernel = {}
        original_links = self.get_link_by_index()

        if Analysis_fk:
            kwargs['Analysis'] = "spp_dyn"
        if "index" in kwargs:
            for i,j in enumerate(kwargs["index"]):
                filtered_kernel[j] = original_links[j]
            return filtered_kernel
        if 'start' in kwargs:
            kwargs['start'] = Time(kwargs['start']).tdb.strftime('%Y%m%d')
        if 'end' in kwargs and kwargs["end"] is not None:
            kwargs['end'] = Time(kwargs['end']).tdb.strftime('%Y%m%d')
        if "version" in kwargs:
            kwargs["version"] = "v" + str(kwargs["version"])


        for index, link in original_links.items():
            match = None


            match = all(value in link for value in kwargs.values())

            if match:
                filtered_kernel[index] = link
        if not len(filtered_kernel):
            print("no match found for the search terms!")

        return filtered_kernel
