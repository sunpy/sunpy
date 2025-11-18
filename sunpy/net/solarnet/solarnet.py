import os
import json
import urllib.request
from pathlib import Path

import requests

from sunpy.net import attrs as a
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.solarnet.attrs import Dataset, Tags, Target, walker

_BASE_URL = "https://solarnet2.oma.be/service/api/svo/{}"

__all__ = ["SOLARNETClient"]

class SOLARNETClient(BaseClient):
    """
    Provides access to query and download from the SOLARNET Virtual Observatory (SVO).

    The SVO was first supported by the `SOLARNET <http://solarnet-east.eu/>`__ project,
    funded by the European Commission's FP7 Capacities Programme under the Grant Agreement 312495.
    It was turned operational due to the `SOLARNET2 <https://solarnet-project.eu/>`__ project, funded by
    the European Union's Horizon 2020 Research and Innovation Programme under Grant Agreement 824135.
    It's purpose is to collect metadata from as many solar observations as possible,
    especially those involved in the SOLARNET projects, in a common catalog and make them available to
    the scientific community.

    There is various ways one can access SVO `Data <https://solarnet.oma.be/#introduction>`__, this
    client using the RESTful API (`available documentation <https://solarnet.oma.be/svo_restful_api_user_manual.html>`__).

    Notes
    -----
    This client by defaults limits the number of results returned from a search to 20.
    To change this limit, use the ``a.solarnet.Limit`` attribute in your query.

    In addition, the query returns a large number of columns, many of which are not useful to the user.
    These are columns found in the FITS headers of the files, and are hidden by default.
    If you want to see all the columns, you can do ``query_results.show()`` on the results
    returned from a search.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a

    >>> search_results = Fido.search(a.solarnet.Dataset.gris_level_1,a.solarnet.Limit(2), a.solarnet.Target.ar)   # doctest: +REMOTE_DATA
    >>> search_results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the SOLARNETClient:
    Source: https://solarnet2.oma.be
    <BLANKLINE>
         OID               DATE_BEG                 DATE_END          WAVEMIN   WAVEMAX          DATE_OBS         EXPTIME MEASURE STEPS WAVELENG
    -------------- ------------------------ ------------------------ --------- --------- ------------------------ ------- ------- ----- --------
    20140426094659 2014-04-26T09:46:59.000Z 2014-04-26T09:53:41.300Z 1564.0199 1568.0461 2014-04-26T00:00:00.000Z    None          None     1566
    20140426095643 2014-04-26T09:56:43.000Z 2014-04-26T10:09:50.300Z 1564.0191 1568.0461 2014-04-26T00:00:00.000Z    None          None     1566
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
        """
        results = []
        query = and_(*query)
        block = walker.create(query)[0]
        if "datasets" in block:
            url = _BASE_URL.format(block["datasets"])
            block.pop("datasets")
        try:
            req = urllib.request.urlopen(url + '?' + urllib.parse.urlencode(block))
            data = json.loads(req.read().decode())["objects"]
        except Exception as e:
            raise OSError(f"Failed to fetch data from {url}") from e
        # It returns so many keys that are not useful for the user
        allowed_keys = [
            "oid", "date_beg", "date_end", "wavemin", "wavemax", "date_obs",
            "exptime", "measure", "steps", "telescope", "waveleng",
        ]
        hide_keys = [
            key.upper() for key in data[0].keys() if key.lower() not in allowed_keys
        ]
        # Just to avoid empty tags
        hide_keys.append("Tags") if "tag" not in block else hide_keys
        for i in data:
            # Filter out keys that should be hidden or not displayed
            filtered_data = {key.upper(): value for key, value in i.items()}
            # Add specific fields from data_location
            filtered_data.update({
                "F_SIZE": i["data_location"]["file_size"],
                "URL": i["data_location"]["file_url"],
            })
            results.append(filtered_data)
        qrt = QueryResponseTable(results, client=self)
        qrt.hide_keys = hide_keys + ["URL", "F_SIZE"]
        return qrt

    def fetch(self, query_results, path, downloader, **kwargs):
        for row in query_results:
            fname = row['URL'].split('/')[-1]
            filepath = str(path).format(file=fname)
            max_splits = kwargs.get('max_splits', 5)
            downloader.enqueue_file(row['URL'], filename=filepath, max_splits=max_splits)

    @staticmethod
    def load_solarnet_values():
        """
        Load Solarnet dataset values from a local JSON file.

        Returns
        -------
        attrs : `dict`
            Dictionary of Solarnet dataset values.
        """
        data_sets = Path(__file__).parent / "data" / "datasets.json"
        tags = Path(__file__).parent / "data" / "tags.json"
        with data_sets.open() as data_values:
            data = json.load(data_values)
        with tags.open() as tags_values:
            tags = json.load(tags_values)
        data_set = {Dataset: list(data.items())}
        target = {Target: [("AR", "Active Region"),("CH","Coronal Hole"),("FS","Flare"),("QR","Quiet Region")]}
        tags = {Tags:list(tags.items())}
        attrs = data_set | target | tags
        return attrs

    @staticmethod
    def _update_solarnet_data():
        """
        Gets the names of available datasets, tags, and targets from SolarNet and saves them into
        separate JSON files.
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))
        data_dir = os.path.join(dir_path, 'data')
        os.makedirs(data_dir, exist_ok=True)
        # Datasets
        url_dataset = _BASE_URL.format("dataset")
        response = requests.get(url_dataset, params={"limit": 100})
        data = response.json()
        values = {}
        for obj in data.get("objects", []):
            name = obj["name"].replace(" ", "_").lower()
            description = obj.get("description") or obj.get("telescope", {}).get("description", "")
            values[name] = description.split(". ")[0].split(", ")[0] if description else ""
        with open(os.path.join(data_dir, 'datasets.json'), 'w') as f:
            json.dump(dict(sorted(values.items())), f, indent=2)
        # Tags
        url_tags = _BASE_URL.format("tag")
        response = requests.get(url_tags)
        data_tags = response.json()
        tag_values = {}
        for obj in data_tags.get("objects", []):
            tag_name = obj["name"].replace(" ", "_").lower()
            tag_values[tag_name] = tag_name
        with open(os.path.join(data_dir, 'tags.json'), 'w') as f:
            json.dump(dict(sorted(tag_values.items())), f, indent=2)

    @classmethod
    def register_values(cls):
        return cls.load_solarnet_values()

    @classmethod
    def _attrs_module(cls):
        return 'solarnet', 'sunpy.net.solarnet.attrs'

    @classmethod
    def _can_handle_query(cls, *query):
        return any(isinstance(q, a.solarnet.Dataset) for q in query)
