import os
import json
import pathlib

import requests

from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.solarnet.attrs import Dataset, walker

_BASE_URL = "https://solarnet2.oma.be/service/api/svo/{}"

__all__ = ["SOLARNETClient"]

class SOLARNETClient(BaseClient):
    """
    Provides access to query and download from the SOLARNET Virtual Observatory (SVO).

    The SVO was first supported by the `SOLARNET <http://solarnet-east.eu/>`__ project,
    funded by the European Commission’s FP7 Capacities Programme under the Grant Agreement 312495.
    It was turned operational due to the `SOLARNET2 <https://solarnet-project.eu/>`__ project, funded by the European Union's Horizon 2020
    Research and Innovation Programme under Grant Agreement 824135.
    It's purpose is to collect metadata from as many solar observations as possible,
    especially those involved in the SOLARNET projects, in a common catalog and make them available to the scientific community.

    There is various ways one can access SVO `Data <https://solarnet.oma.be/#introduction>`__, this client using the
    RESTful API (`available documentation <https://solarnet.oma.be/svo_restful_api_user_manual.html>`__).

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> query = [a.solarnet.Dataset.eui_level_2 , a.solarnet.Limit(2) , a.Detector("HRI_EUV")]
    >>> search_results = Fido.search(*query)   # doctest: +REMOTE_DATA
    >>> search_results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the SOLARNETClient:
    Source: https://solarnet2.oma.be
    <BLANKLINE>
          datasets                              name                        detector
    -------------------- -------------------------------------------------- --------
    metadata_eui_level_2 solo_L2_eui-hrieuvopn-image_20200512T122556952_V06  HRI_EUV
    metadata_eui_level_2 solo_L2_eui-hrieuvopn-image_20200512T122606952_V06  HRI_EUV
    <BLANKLINE>
    <BLANKLINE>
    """

    @property
    def info_url(self):
        return 'https://solarnet2.oma.be'

    def search(self, *query):
        """
        Searches the SVO for datasets matching the query.

        Parameters
        ----------
        *query : `~sunpy.net.attrs`
            Attributes to search for datasets.

        Returns
        -------
        `~sunpy.net.base_client.QueryResponseTable`
            A table containing the search results.

        Examples
        --------
        >>> from sunpy.net import Fido, attrs as a
        >>> query = [a.solarnet.Dataset.lyra_level_2 , a.solarnet.Limit(3) , a.Detector("HRI_EUV")]
        >>> search_results = Fido.search(*query)   # doctest: +REMOTE_DATA
        >>> search_results  # doctest: +REMOTE_DATA
        <sunpy.net.fido_factory.UnifiedResponse object at ...>
        Results from 1 Provider:
        <BLANKLINE>
        3 Results from the SOLARNETClient:
        Source: https://solarnet2.oma.be
        <BLANKLINE>
               datasets                    name             detector
        --------------------- ----------------------------- --------
        metadata_lyra_level_2 lyra_20100106-000000_lev2_std  HRI_EUV
        metadata_lyra_level_2 lyra_20100107-000000_lev2_std  HRI_EUV
        metadata_lyra_level_2 lyra_20100108-000000_lev2_std  HRI_EUV
        <BLANKLINE>
        <BLANKLINE>
        """
        results = []
        query = and_(*query)
        block = walker.create(query)[0]
        if "datasets" in block:
            url = _BASE_URL.format(block["datasets"])
            source = block.pop("datasets")
        self.links = self._generate_links(block, url)

        for i in range(len(self.links)):
            results.append({
                "datasets": source,
                "name": os.path.splitext(os.path.basename(self.links[i]))[0],
                "detector": block["detector__iexact"] if "detector__iexact" in block else "Not specified",
                "url" : self.links[i],
            })

        qrt = QueryResponseTable(results, client=self)
        qrt.hide_keys = ["url"]
        return qrt

    def _generate_links(self, block, url):
        """
        Helper function to generate links based on query blocks.

        Parameters
        ----------
        block : `dict`
            Query block parameters.
        url : `str`
            Base URL for the query.

        Returns
        -------
        links : `list`
            List of dataset file URLs.
        """
        links = []
        req = requests.get(url,params = block)
        data = req.json()["objects"]
        for i in data:
            link = i["data_location"]["file_url"]
            links.append(link)
        return links

    def fetch(self, query_results, downloader, path=None, **kwargs):
        """
        Download a set of results.

        Parameters
        ----------
        query_results : `~sunpy.net.base_client.QueryResponseTable`
            Results from a Fido search.
        path : `str` or `pathlib.Path`, optional
            Path to save data to, defaults to SunPy download dir
        downloader : `parfive.Downloader`
            Downloader instance used to download data.
        """
        for row in query_results:
            fname = row['url'].split('/')[-1]
            filepath = str(path).format(file=fname)
            downloader.enqueue_file(row['url'], filename=filepath)

    @staticmethod
    def load_solarnet_values():
        """
        Load Solarnet dataset values from a local JSON file.

        Returns
        -------
        attrs : `dict`
            Dictionary of Solarnet dataset values.
        """
        data_sets = pathlib.Path(__file__).parent / "data" / "datasets.json"
        with data_sets.open() as data_values:
            data = json.load(data_values)
        return {Dataset: list(data.items())}

    @staticmethod
    def create_parse_solarnet_values():
        """
        Gets the names of available datasets using https://solarnet.oma.be/service/api/svo/dataset.
        """
        dir = os.path.dirname(os.path.realpath(__file__))
        url = _BASE_URL.format("dataset")
        response = requests.get(url, params={"limit": 100})
        data = response.json()
        values = {}
        for obj in data.get("objects", []):
            name = obj["name"].replace(" ", "_").lower()
            description =  obj.get("description") or  obj.get("telescope", {}).get("description")
            values[name] = description.split(". ")[0].split(", ")[0]
        with open(os.path.join(dir, 'data', 'datasets.json'), 'w') as attrs_file:
            json.dump(dict(sorted(values.items())), attrs_file, indent=2)

    @classmethod
    def register_values(cls):
        return cls.load_solarnet_values()

    @classmethod
    def _attrs_module(cls):
        return 'solarnet', 'sunpy.net.solarnet.attrs'

    @staticmethod
    def _can_handle_query(*query):
        from sunpy.net import attrs as a

        # Checks for dataset instance in query
        return any(isinstance(q, a.solarnet.Dataset) for q in query)
