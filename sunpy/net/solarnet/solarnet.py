import os
from pathlib import Path

import requests

from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.solarnet.attrs import walker
from sunpy.util.parfive_helpers import Downloader, Results

base_url = "https://solarnet2.oma.be/service/api/svo/{}"


class SolarnetClient(BaseClient):

    def search(self,*query,**kwargs):
        results = []
        query = and_(*query)
        block = walker.create(query)[0]
        if "datasets" in block:
            url = base_url.format(block["datasets"])
            source = block.pop("datasets")
        links = self._generate_links(block,url)

        for i in range(len(links)):
            results.append({
                "datasets":source,
            })
        return QueryResponseTable(results, client=self)

    def _generate_links(self,block,url):
        links = []
        req = requests.get(url,params=block)
        data = req.json()["objects"]
        for i in data:
            link = i["data_location"]["file_url"]
            links.append(link)
            print(link)
        return links

    def downloader(self,link,overwrite = False,progress = True,wait = True,downloader = None, path = None,):
        downloader = Downloader(progress = progress,overwrite= overwrite,max_splits=1)
        if path is None:
            default_dir = "Downloads" + f"_{link}"
            path = os.path.join(default_dir,'{file}')

        elif isinstance(path,Path):
            path = str(path)
        if isinstance(path,str) and '{file}' not in path:
            path = os.path.join(path,'{file}')
        file_name = file_name = path.format(file = link)
        downloader.enqueue_file(link, path=file_name, filename=file_name,max_splits=1)

        if not wait:
            return Results()
        results =  downloader.download()
        return results

    def fetch(self):
        pass

    @staticmethod
    def _can_handle_query():
        return True
