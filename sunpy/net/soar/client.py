import requests

import astropy.table

from sunpy import log
import sunpy.net.attrs as a
from sunpy.net.attr import and_
from sunpy.net.base_client import BaseClient, QueryResponseTable
from sunpy.time import parse_time

from sunpy.net.soar.attrs import Identifier, walker


__all__ = ['SOARClient']


class SOARClient(BaseClient):
    """
    Client to access the Solar Orbiter Archive (SOAR).
    """
    def search(self, *query, **kwargs):
        query = and_(*query)
        queries = walker.create(query)

        results = []
        for query_parameters in queries:
            results.append(self._do_search(query_parameters))
        table = astropy.table.vstack(results)
        qrt = QueryResponseTable(table, client=self)
        filesize = qrt['Filesize'] / (1024 * 1024)
        qrt['Filesize'] = [f'{f:.1f} MB' for f in filesize]
        qrt.hide_keys = ['Data item ID', 'Filename']
        return qrt

    @staticmethod
    def _do_search(query):
        """
        Query the SOAR server with a single query.

        Parameters
        ----------
        query : list[str]
            List of query items.

        Returns
        -------
        astropy.table.QTable
            Query results.
        """
        base_url = ('http://soar.esac.esa.int/soar-sl-tap/tap/'
                    'sync?REQUEST=doQuery&')
        # Need to manually set the intervals based on a query
        request_dict = {}
        request_dict['LANG'] = 'ADQL'
        request_dict['FORMAT'] = 'json'

        url_query = {}
        url_query['SELECT'] = '*'
        url_query['FROM'] = 'v_sc_data_item'
        url_query['WHERE'] = '+AND+'.join(query)
        request_dict['QUERY'] = '+'.join([f'{item}+{url_query[item]}' for
                                          item in url_query])

        request_str = ''
        request_str = [f'{item}={request_dict[item]}' for item in request_dict]
        request_str = '&'.join(request_str)

        url = base_url + request_str
        log.debug(f'Getting request from URL: {url}')
        # Get request info
        r = requests.get(url)
        r.raise_for_status()

        # Do some list/dict wrangling
        names = [m['name'] for m in r.json()['metadata']]
        info = {name: [] for name in names}
        for entry in r.json()['data']:
            for i, name in enumerate(names):
                info[name].append(entry[i])

        if len(info['begin_time']):
            info['begin_time'] = parse_time(info['begin_time']).iso
            info['end_time'] = parse_time(info['end_time']).iso

        return astropy.table.QTable({'Instrument': info['instrument'],
                                     'Data product': info['descriptor'],
                                     'Level': info['level'],
                                     'Start time': info['begin_time'],
                                     'End time': info['end_time'],
                                     'Data item ID': info['data_item_id'],
                                     'Filename': info['filename'],
                                     'Filesize': info['filesize']
                                     })

    def fetch(self, query_results, *, path, downloader, **kwargs):
        """
        Queue a set of results to be downloaded. `BaseClient` does the actual
        downloading, so we just have to queue up the ``downloader``.

        Parameters
        ----------
        query_results : sunpy.net.fido_factory.UnifiedResponse
            Results from a Fido search.
        path : str
            Path to download files to. Must be a format string with a ``file``
            field for the filename.
        downloader : parfive.Downloader
            Downloader instance used to download data.
        kwargs :
            Keyword arguments aren't used by this client.
        """
        base_url = ('http://soar.esac.esa.int/soar-sl-tap/data?'
                    f'retrieval_type=LAST_PRODUCT&product_type=SCIENCE&'
                    'data_item_id=')

        for row in query_results:
            url = base_url + row['Data item ID']
            filepath = str(path).format(file=row['Filename'])
            log.debug(f'Queing URL: {url}')
            downloader.enqueue_file(url, filename=filepath)

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Check if this client can handle a given Fido query.

        Returns
        -------
        bool
            True if this client can handle the given query.
        """
        required = {a.Time}
        optional = {a.Instrument, a.Level, Identifier}
        return cls.check_attr_types_in_query(query, required, optional)
