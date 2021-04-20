import requests

import astropy.table

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
        qrt.hide_keys = ['Data item ID', 'Filename']
        return qrt

    @staticmethod
    def _do_search(query):
        """
        Carry out a search with a single query.
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
                                     })

    def fetch(self, query_results, *, path, downloader, **kwargs):
        base_url = ('http://soar.esac.esa.int/soar-sl-tap/data?'
                    f'retrieval_type=PRODUCT&product_type=SCIENCE&'
                    'data_item_id=')

        for row in query_results:
            url = base_url + row['Data item ID']
            filepath = str(path).format(file=row['Filename'])
            downloader.enqueue_file(url, filename=filepath)

    @classmethod
    def _can_handle_query(cls, *query):

        required = {a.Time}
        optional = {a.Instrument, a.Level, Identifier}
        return cls.check_attr_types_in_query(query, required, optional)
