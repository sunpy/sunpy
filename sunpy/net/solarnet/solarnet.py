import os
import json
import pathlib
from pathlib import Path

import requests

from sunpy.net import attrs as a
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.solarnet.attrs import walker
from sunpy.util.parfive_helpers import Downloader, Results

base_url = "https://solarnet2.oma.be/service/api/svo/{}"


class SolarnetClient(BaseClient):

    def search(self,*query):
        results = []
        query = and_(*query)
        block = walker.create(query)[0]
        if "datasets" in block:
            url = base_url.format(block["datasets"])
            source = block.pop("datasets")
        self.links = self._generate_links(block,url)

        for i in range(len(self.links)):
            results.append({
                "index":i,
                "datasets":source,
                "name":os.path.splitext(os.path.basename(self.links[i]))[0]
            })
        return QueryResponseTable(results, client=self)

    def _generate_links(self,block,url):
        links = []
        req = requests.get(url,params=block)
        data = req.json()["objects"]
        for i in data:
            link = i["data_location"]["file_url"]
            links.append(link)
        return links

    def downloader(self,link,overwrite = False,progress = True,wait = True,downloader = None, path = None,):
        downloader = Downloader(progress = progress,overwrite= overwrite,max_splits=1)
        link_name = link.split('/')[-1]
        if path is None:
            default_dir = "Downloads"
            path = os.path.join(default_dir,'{file}')

        elif isinstance(path,Path):
            path = str(path)
        if isinstance(path,str) and '{file}' not in path:
            path = os.path.join(path,'{file}')
        file_name = file_name = path.format(file = link_name)
        downloader.enqueue_file(link,filename=file_name,max_splits=1)

        if not wait:
            return Results()
        results =  downloader.download()
        return results

    def fetch(self, query_results, overwrite = False,path=None, **kwargs):
        if len(query_results) == 3:
            index = query_results[0]
            self.downloader(self.links[index],overwrite=overwrite, progress=True, wait=True, path=path)
        else:
            for i in query_results:
                index = i["index"]
                self.downloader(self.links[index],overwrite=overwrite, progress=True, wait=True, path=path)

    @staticmethod
    def load_solarnet_values():

        data_sets = pathlib.Path(__file__).parent / "data" / "datasets.json"
        with data_sets.open() as data_values:
            data = json.load(data_values)
        data = list(data.items())
        attrs = {
            a.solarnet.Datasets:data,
        }

        return attrs

    @classmethod
    def register_values(cls):
        """
        Register the solarnet specific attributes with Fido.

        Returns
        -------
        dict
            The dictionary containing the values formed into attributes.
        """
        return cls.load_solarnet_values()

    @classmethod
    def _attrs_module(cls):
        return 'solarnet', 'sunpy.net.solarnet.attrs'

    @staticmethod
    def _can_handle_query(*query):
        from sunpy.net import attrs as a
        return any(isinstance(q, a.solarnet.Datasets) for q in query)
