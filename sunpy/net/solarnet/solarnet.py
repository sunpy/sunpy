
from sunpy.net.base_client import BaseClient
from sunpy.net.solarnet import attrs as a
import requests

base_url = "https://solarnet2.oma.be/service/api/svo/{}"


class SolarnetClient(BaseClient):

    def search(self,*query,**kwargs):
        query = and_(*query)
        block = a.walker.create(query)[0]
        if "datasets" in block:
            url = base_url.format(block["datasets"])
            block.pop("datasets")
        req = requests.get(url,params=block)
        data = req.json()["objects"]
        for i in data:
            link = i["data_location]["file_url"]
            
