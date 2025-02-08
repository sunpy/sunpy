import os
import json
import pathlib
from pathlib import Path

import requests

from sunpy import config
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.solarnet.attrs import Dataset, walker
from sunpy.util.parfive_helpers import Downloader, Results

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
    Source: https://solarnet2.oma.be/
    <BLANKLINE>
    index       datasets                              name
    ----- -------------------- --------------------------------------------------
        0 metadata_eui_level_2 solo_L2_eui-hrieuvopn-image_20200512T122556952_V06
        1 metadata_eui_level_2 solo_L2_eui-hrieuvopn-image_20200512T122606952_V06
    <BLANKLINE>
    <BLANKLINE>
    """

    @property
    def info_url(self):
        return 'https://solarnet2.oma.be/'

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
        Source: https://solarnet2.oma.be/
        <BLANKLINE>
        index        datasets                    name
        ----- --------------------- -----------------------------
            0 metadata_lyra_level_2 lyra_20100106-000000_lev2_std
            1 metadata_lyra_level_2 lyra_20100107-000000_lev2_std
            2 metadata_lyra_level_2 lyra_20100108-000000_lev2_std
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
                "index": i,
                "datasets": source,
                "name": os.path.splitext(os.path.basename(self.links[i]))[0]
            })
        return QueryResponseTable(results, client=self)

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

    def _downloader(self, link, overwrite=False, progress=True, wait=True, path=None):
        """
        Download a file from a given link.

        Parameters
        ----------
        link : `str`
            URL of the file to download.
        overwrite : `bool`, optional
            Whether to overwrite existing files (default is False).
        progress : `bool`, optional
            Show download progress (default is True).
        wait : `bool`, optional
            Wait for the download to finish (default is True).
        path : `str` or `Pathlib.Path`
            Path to save data to, defaults to SunPy download dir

        Returns
        -------
        Results : `parfive.Results`
            A `parfive.Results` instance or `None` if no URLs to download
        """
        downloader = Downloader(progress=progress, overwrite=overwrite, max_splits=1)
        link_name = link.split('/')[-1]
        if path is None:
            default_dir = config.get("downloads", "download_dir")
            path = os.path.join(default_dir, '{file}')
        elif isinstance(path, Path):
            path = str(path)
        if isinstance(path, str) and '{file}' not in path:
            path = os.path.join(path, '{file}')

        file_name = path.format(file=link_name)
        downloader.enqueue_file(link, filename=file_name, max_splits=1)

        if not wait:
            return Results()
        results = downloader.download()
        return results

    def fetch(self, query_results, overwrite=False, path=None, **kwargs):
        """
        Fetch one or more datasets based on query results.

        Parameters
        ----------
        query_results : `~sunpy.net.base_client.QueryResponseTable`
            Query results to fetch.
        overwrite : `bool`, optional
            Whether to overwrite existing files (default is False).
        path : `str` or `pathlib.Path`, optional
            Path to save data to, defaults to SunPy download dir

        Returns
        -------
        Results : `parfive.Results`
            A `parfive.Results` instance or `None` if no URLs to download
        """
        indices = (
            [query_results[0]] if len(query_results) == 3
            else [i["index"] for i in query_results]
        )
        for index in indices:
            self._downloader(self.links[index], overwrite=overwrite, path=path)

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
            description =  obj.get("instrument", {}).get("description") or  obj.get("telescope", {}).get("description")
            values[name] = description.split(". ")[0]
        with open(os.path.join(dir, 'data', 'datasets.json'), 'w') as attrs_file:
            json.dump(dict(sorted(values.items())), attrs_file, indent=2)

    @classmethod
    def register_values(cls):
        #loads the solarnet values
        return cls.load_solarnet_values()

    @classmethod
    def _attrs_module(cls):
        return 'solarnet', 'sunpy.net.solarnet.attrs'

    @staticmethod
    def _can_handle_query(*query):
        from sunpy.net import attrs as a

        # Checks for dataset instance in query
        return any(isinstance(q, a.solarnet.Dataset) for q in query)
