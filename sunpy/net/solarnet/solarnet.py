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
    """
    Provides access to query and download from Solarnet API.

    Methods
    -------
    search(*query)
        Search the Solarnet database for datasets matching the query.
    downloader(link, overwrite=False, progress=True, wait=True, path=None)
        Download a file from a given link.
    fetch(query_results, overwrite=False, path=None, **kwargs)
        Fetch one or more datasets based on query results.
    load_solarnet_values()
        Load Solarnet dataset values from a local JSON file.
    register_values()
        Register Solarnet-specific attributes with Fido.
    _can_handle_query(*query)
        Determine if the query can be handled by this client.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> query = [a.solarnet.Dataset.eui_level_2 , a.solarnet.Limit(2) , a.Detector("HRI_EUV")]  # doctest: +SKIP
    >>> url = Fido.search(*query)   # doctest: +SKIP
    >>> print(url)  # doctest: +SKIP
    Index    datasets                              name
    ----- -------------------- --------------------------------------------------
        0 metadata_eui_level_2 solo_L2_eui-hrieuvopn-image_20200512T122556952_V06
        1 metadata_eui_level_2 solo_L2_eui-hrieuvopn-image_20200512T122606952_V06
    """

    def search(self, *query):
        """
        Search the Solarnet database for datasets matching the query.

        Parameters
        ----------
        *query : sunpy.net.attrs
            Attributes to search for datasets.

        Returns
        -------
        QueryResponseTable
            Table of search results.

        Examples
        --------
        >>> from sunpy.net import Fido, attrs as a
        >>> query = [a.solarnet.Dataset.eui_level_2 , a.solarnet.Limit(2) , a.Detector("HRI_EUV")]  # doctest: +SKIP
        >>> url = Fido.search(*query)   # doctest: +SKIP
        >>> print(url)  # doctest: +SKIP
        Index    datasets                              name
        ----- -------------------- --------------------------------------------------
            0 metadata_eui_level_2 solo_L2_eui-hrieuvopn-image_20200512T122556952_V06
            1 metadata_eui_level_2 solo_L2_eui-hrieuvopn-image_20200512T122606952_V06
        """
        results = []
        query = and_(*query)
        block = walker.create(query)[0]
        if "datasets" in block:
            url = base_url.format(block["datasets"])
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
        block : dict
            Query block parameters.
        url : str
            Base URL for the query.

        Returns
        -------
        list
            List of dataset file URLs.
        """
        links = []
        req = requests.get(url, params=block)
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
        link : str
            URL of the file to download.
        overwrite : bool, optional
            Whether to overwrite existing files (default is False).
        progress : bool, optional
            Show download progress (default is True).
        wait : bool, optional
            Wait for the download to finish (default is True).
        path : str or Path, optional
            Path template for saving the downloaded file.

        Returns
        -------
        Results
            Result of the download operation.
        """
        downloader = Downloader(progress=progress, overwrite=overwrite, max_splits=1)
        link_name = link.split('/')[-1]
        if path is None:
            default_dir = "Downloads"
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
        query_results : QueryResponseTable
            Query results to fetch.
        overwrite : bool, optional
            Whether to overwrite existing files (default is False).
        path : str or Path, optional
            Path template for saving the downloaded files.
        """
        if len(query_results) == 3:
            index = query_results[0]
            self._downloader(self.links[index], overwrite=overwrite, path=path)
        else:
            for i in query_results:
                index = i["index"]
                self._downloader(self.links[index], overwrite=overwrite, path=path)

    @staticmethod
    def load_solarnet_values():
        """
        Load Solarnet dataset values from a local JSON file.

        Returns
        -------
        dict
            Dictionary of Solarnet dataset values.
        """
        data_sets = pathlib.Path(__file__).parent / "data" / "datasets.json"
        with data_sets.open() as data_values:
            data = json.load(data_values)
        data = list(data.items())
        attrs = {a.solarnet.Dataset: data}
        return attrs

    @staticmethod
    def create_parse_solarnet_values():
        """
        To get the names of available datasets.
        Make a GET request to https://solarnet.oma.be/service/api/svo/dataset to do so.
        """

        dir = os.path.dirname(os.path.realpath(__file__))
        url = "https://solarnet.oma.be/service/api/svo/dataset"
        response = requests.get(url, params={"limit": 100})
        data = response.json()
        names = [obj["name"].replace(" ", "_").lower() for obj in data.get("objects", [])]
        values = {name: name for name in names}
        with open(os.path.join(dir, 'data', 'datasets.json'), 'w') as attrs_file:
            json.dump(dict(sorted(values.items())), attrs_file, indent=2)

    @classmethod
    def register_values(cls):
        """
        Register the Solarnet-specific attributes with Fido.

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
        """
        Determine if the query can be handled by this client.

        Returns
        -------
        bool
            True if the query can be handled, False otherwise.
        """
        from sunpy.net import attrs as a
        return any(isinstance(q, a.solarnet.Dataset) for q in query)
