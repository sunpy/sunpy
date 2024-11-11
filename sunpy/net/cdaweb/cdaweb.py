import json
import pathlib

import requests

import astropy.table
from astropy.time import Time

import sunpy.net.attrs as a
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.net.cdaweb.attrs import Dataset
from .walker import walker

__all__ = ['CDAWEBClient']

_CDAS_BASEURL = 'https://cdaweb.gsfc.nasa.gov/WS/cdasr/1'
_CDAS_HEADERS = {'Accept': 'application/json'}
_CDAS_TIME_FMT = '%Y%m%dT%H%M%SZ'
_DATAVIEW = 'sp_phys'


class CDAWEBClient(BaseClient):
    """
    Provides access to query and download from the Coordinated Data Analysis Web (CDAWeb).

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> print(a.cdaweb.Dataset)
    sunpy.net.cdaweb.attrs.Dataset
    <BLANKLINE>
    Dataset ID.
    <BLANKLINE>
    <BLANKLINE>
                         Attribute Name                     Client ...                                   Description
    ------------------------------------------------------- ------ ... --------------------------------------------------------------------------------
    a1_k0_mpa                                               CDAWEB ... LANL 2001 Magnetospheric Plasma Analyzer Key Parameters - Mike Henderson (LANL)
    a2_k0_mpa                                               CDAWEB ... LANL 2002 Magnetospheric Plasma Analyzer Key Parameters - Mike Henderson (LANL)
    ...
    >>> res = Fido.search(a.Time('2021/07/01', '2021/07/08'),
    ...                   a.cdaweb.Dataset('SOLO_L2_MAG-RTN-NORMAL-1-MINUTE')) # doctest: +REMOTE_DATA
    >>> res # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    7 Results from the CDAWEBClient:
    Source: https://cdaweb.gsfc.nasa.gov/index.html
    <BLANKLINE>
               Dataset                    Start time               End time
    ------------------------------- ----------------------- -----------------------
    SOLO_L2_MAG-RTN-NORMAL-1-MINUTE 2021-07-01 00:00:29.000 2021-07-01 23:59:30.000
    SOLO_L2_MAG-RTN-NORMAL-1-MINUTE 2021-07-02 00:00:29.000 2021-07-02 23:59:30.000
    SOLO_L2_MAG-RTN-NORMAL-1-MINUTE 2021-07-03 00:00:29.000 2021-07-03 23:59:30.000
    SOLO_L2_MAG-RTN-NORMAL-1-MINUTE 2021-07-04 00:00:29.000 2021-07-04 23:59:30.000
    SOLO_L2_MAG-RTN-NORMAL-1-MINUTE 2021-07-05 00:00:29.000 2021-07-05 23:59:30.000
    SOLO_L2_MAG-RTN-NORMAL-1-MINUTE 2021-07-06 00:00:29.000 2021-07-06 23:59:30.000
    SOLO_L2_MAG-RTN-NORMAL-1-MINUTE 2021-07-07 00:00:29.000 2021-07-07 23:59:30.000
    <BLANKLINE>
    <BLANKLINE>

    See Also
    --------
    sunpy.net.cdaweb.get_datasets : Find dataset IDs for a given observatory
    sunpy.net.cdaweb.get_observatory_groups : Get all observatories available from CDAWeb
    """
    @property
    def info_url(self):
        return 'https://cdaweb.gsfc.nasa.gov/index.html'

    def search(self, *query, **kwargs):
        """
        Search for datasets provided by the Space Physics Data Facility.
        """
        query = and_(*query)
        queries = walker.create(query)

        results = []
        for query_parameters in queries:
            results.append(self._do_search(query_parameters))
        table = astropy.table.vstack(results)
        qrt = QueryResponseTable(table, client=self)
        qrt.hide_keys = ['URL']
        return qrt

    def _do_search(self, query):
        response = (self._get_remote_files(query['dataset'],
                                           query['begin_time'],
                                           query['end_time']))

        if 'FileDescription' not in response:
            # No results
            return astropy.table.QTable(
                {'Dataset': [],
                 'Start time': [],
                 'End time': [],
                 'URL': []})
        else:
            stimes = [f['StartTime'] for f in response['FileDescription']]
            etimes = [f['EndTime'] for f in response['FileDescription']]
            urls = [f['Name'] for f in response['FileDescription']]
            return astropy.table.QTable(
                {'Dataset': [query['dataset']] * len(stimes),
                 'Start time': Time.strptime(stimes, '%Y-%m-%dT%H:%M:%S.%fZ').iso,
                 'End time': Time.strptime(etimes, '%Y-%m-%dT%H:%M:%S.%fZ').iso,
                 'URL': urls})

    @staticmethod
    def _get_remote_files(dataset, start, end):
        # Get a list of files for a given dataset between start and end times
        start = start.strftime(_CDAS_TIME_FMT)
        end = end.strftime(_CDAS_TIME_FMT)
        url = '/'.join([
            _CDAS_BASEURL,
            'dataviews', _DATAVIEW,
            'datasets', dataset,
            'orig_data', f'{start},{end}'
        ])
        response = requests.get(url, headers=_CDAS_HEADERS)
        return response.json()

    def fetch(self, query_results, *, path, downloader, **kwargs):
        for row in query_results:
            fname = row['URL'].split('/')[-1]
            filepath = str(path).format(file=fname)
            # Manually cap max_splits at 3 to make CDAWeb happy
            max_splits = kwargs.get('max_splits', 3)
            max_splits = min(max_splits, 3)
            downloader.enqueue_file(row['URL'], filename=filepath, max_splits=max_splits)

    @classmethod
    def _can_handle_query(cls, *query):
        required = {Dataset, a.Time}
        query_attrs = {type(x) for x in query}
        return required == query_attrs

    @classmethod
    def _attrs_module(cls):
        return 'cdaweb', 'sunpy.net.cdaweb.attrs'

    @classmethod
    def register_values(cls):
        return cls.load_dataset_values()

    @staticmethod
    def load_dataset_values():
        from sunpy.net import attrs as a

        attrs_path = pathlib.Path(__file__).parent / 'data' / 'attrs.json'
        with open(attrs_path) as attrs_file:
            all_datasets = json.load(attrs_file)

        # Convert from dict to list of tuples
        all_datasets = [(id, desc) for id, desc in all_datasets.items()]
        return {a.cdaweb.Dataset: all_datasets}
