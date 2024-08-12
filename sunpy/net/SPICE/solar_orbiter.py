import os
import warnings
import urllib.request

from bs4 import BeautifulSoup

from astropy.time import Time

BASE_URL = "https://spiftp.esac.esa.int/data/SPICE/SOLAR-ORBITER/kernels/{}"

KERNEL_TYPES = ["ck", "fk", "ik", "lsk", "pck", "sclk", "spk","mk"]

class SoloKernel:
    def __init__(self, kernel_type,instrument = None,start = None,end = None,**kwargs):
        self.kernel_type = kernel_type
        self.instrument = instrument
        self.kernel_type = kernel_type
        self.start = start
        self.end = end

        if kernel_type not in KERNEL_TYPES:
            raise ValueError(f"Kernel type not recognized '{kernel_type}'. Recognized ones: {KERNEL_TYPES}")

        self.kernel_urls = BASE_URL.format(kernel_type)



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

    def download_by_index(self,*args,destination = "Downloads"):
        """
        allows downloading links with its corresponding index.
        args being index
        """
        if self.kernel_type == "mk":
            warnings.warn("Note you are currently downloading meta kernels")
        else:
            print(f"current kerenel {self.kernel_type}")
        print("Note  aareadme.txt provide description for kernel")
        links_mapped = self.get_link_by_index()
        num_links = len(links_mapped.keys())
        print(f"number of kernels {num_links}")

        for index in args:
            if 0 <= index < num_links:
                file_url = self.kernel_urls + "/" + links_mapped[index]
                if not index:
                    file_name = os.path.join(destination, f"for_{self.kernel_type} " +links_mapped[index])
                else:
                    file_name = os.path.join(destination, links_mapped[index])
                try:
                    urllib.request.urlretrieve(file_url,file_name)
                    print(f"Downloaded:{index} -- {file_name}")

                except Exception as e:
                    print(f"error downloading {file_name} :{e}")


            else :
                raise ValueError(f"index '{index}' must valid between 0 and {num_links-1}")
    def search(self):
        pass


    def filter_kernels(self, logic="and", get_all = False,get_readme = False,**kwargs):
        """
        Filter kernels based on search terms using specified boolean logic.
        Parameters:
        - logic: "and" (default), "or", or "not" to determine filtering behavior.
        - kwargs: key-value pairs to filter by.

        Returns:
        - A dictionary with the filtered kernels retaining the original indices.
        """
        filtered_kernel = {}
        original_links = self.get_link_by_index()

        if get_readme:
            filtered_kernel[0] = self.get_readme()
            return filtered_kernel
        if get_all:
            return original_links
        if 'start' in kwargs:
            kwargs['start'] = Time(kwargs['start']).tdb.strftime('%Y%m%d')
        if 'end' in kwargs:
            kwargs['end'] = Time(kwargs['end']).tdb.strftime('%Y%m%d')
            print(kwargs["end"])
        for index, link in original_links.items():
            match = None

            if logic == "and":
                match = all(value in link for value in kwargs.values())
            elif logic == "or":
                match = any(value in link for value in kwargs.values())
            elif logic == "not":
                match = all(value not in link for value in kwargs.values())
            else:
                raise ValueError("Invalid logic type. Choose 'and', 'or', or 'not'.")

            if match:
                filtered_kernel[index] = link
        if not len(filtered_kernel):
            print("no match found for the search terms!")

        return filtered_kernel
