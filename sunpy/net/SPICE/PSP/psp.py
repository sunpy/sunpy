import os
import urllib.request
from pathlib import Path

from bs4 import BeautifulSoup

from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.util.parfive_helpers import Downloader, Results
from sunpy.net.SPICE.PSP import attrs as ps

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

FILE_CONVENTIONS = {
"sclk" : """ File Naming Convention	spp_sclk_NNNN.tsc (where NNNN is the number of times the file has been updated since launch)
Description	The SCLK is a daily produced kernel that supports time conversions between spacecraft clock and Barycentric Dynamical Time (TDB). The kernel is produced by the PSP mission.
""",

"ltpek" :""" File Naming Convention	spp_nom_yyyymmdd_yyyymmdd_vNNN_###.bsp (where the yyyymmdd pair defines the time span of the file, the NNN is the version number and ### further describes the content of the file)
Description	This kernel contains a nominal PSP trajectory long term predicted ephemeris for the PSP spacecraft. The kernel is produced by the PSP mission.
""",

"reck" :""" File Naming Convention	spp_recon_yyyymmdd_yyyymmdd_vNNN.bsp (where the yyyymmdd pair defines the time span of the file, and the NNN is the version number)
Description	This kernel contains the reconstructed ephemeris for the PSP spacecraft. The kernel is produced by the PSP mission.
""",
"lsk" :""" File Naming Convention	naifNNNN.tls (where NNNN is the NAIF assigned version number)
Description	The kernel is used for Coordinated Universal Time (UTC) to TDB time conversions. It contains a tabulation of all leap seconds that have occurred. It is a generic SPICE kernel, independent of flight project. The kernel is produced by NAIF and is only updated as needed.
""",

"fk" :""" File Naming Convention	spp_vNNN.tf (where NNNN is the PSP assigned version number)
Description	The kernel contains definitions of and specification of relationships between reference frames. This frame kernel contains the current set of coordinate frame definitions for the Parker Solar Probe spacecraft, structures, and science instruments. The kernel is produced by the PSP mission.
""",

"yak":""" File Naming Convention	spp_nom_yyyymmdd_yyyymmdd_vNNN_###_ yyyymmdd_yyyymmdd.bc (where the yyyymmdd pair defines the time span of the file, and the NNN is the version number)
Description	This kernel contains the long-term attitude for the PSP spacecraft. The kernel is produced by the PSP mission.
""",
"ah" :""" File Naming Convention	spp_nom_yyyymmdd_yyyymmdd_vNNN_###_ yyyymmdd_yyyymmdd.bc (where the yyyymmdd pair defines the time span of the file, and the NNN is the version number)
Description	This kernel contains the long-term attitude for the PSP spacecraft. The kernel is produced by the PSP mission.
""",
"pck":""" File Naming Convention	pckNNNNN.tpc (where NNNN is the NAIF assigned version number)
Description	The kernel is used to obtain celestial body orientation, size, shape and other constants. It is a generic SPICE kernel, independent of flight project. The kernel is produced by NAIF.
""",

"pek":""" File Naming Convention	spp_yyyy_doy_NN.ah.bc (where yyyy_doy defines the time span of the data in the file and NN is the version number)
Description	This kernel contains the daily attitude history for the PSP spacecraft. The kernel is produced by the PSP mission.
"""


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
        if kwargs["Analysis"]:
            kwargs["Analysis"] = "spp_dyn"
             


class PSPResponseTable(QueryResponseTable):
    """
    A table for storing psp spice kernels
    """

class PSPClient(BaseClient):
    """
    PASS
    """
