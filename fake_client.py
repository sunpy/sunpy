from functools import partial
from sunpy.net.base_client import BaseClient, convert_row_to_table
import sunpy.util.net

import sunpy.net.attrs as a
import sunpy.net
from sunpy.net.attr import AttrWalker, AttrAnd, AttrOr, DataAttr
from sunpy.net.base_client import QueryResponseTable

walker = AttrWalker()

@walker.add_creator(AttrOr)
def create_or(wlk, tree):
    results = []
    for sub in tree.attrs:
        results.append(wlk.create(sub))

    return results


@walker.add_creator(AttrAnd, DataAttr)
def create_and(wlk, tree):
    result = dict()
    wlk.apply(tree, result)
    return [result]


@walker.add_applier(a.Time)
def _(wlk, attr, params):
    return params.update({'startTime': attr.start.isot,
                            'endTime': attr.end.isot})


@walker.add_applier(a.Level)
def _(wlk, attr, params):
    return params.update({'level': attr.value})


class ExampleClient(BaseClient):
    size_column = 'Filesize'

    def search(self, query):
        queries = walker.create(query)
        print("Seaching in fake client")
        results = []
        for query_parameters in queries:
            results.append(self._make_search(query_parameters))
            # try:
            #     resp = self._make_search(query_parameters)
            #     results.append(resp)
            # except ConnectionRefusedError as e:
            #     # e.client_name = client.__class__.__name__.lower()
            #     results.append(e)

        return QueryResponseTable(results, client=self)
    
    @classmethod
    def register_values(cls):

        from sunpy.net import attrs
        adict = {
        attrs.Instrument: [("LASCO", "Large Angle and Spectrometric Coronagraph")],
        attrs.Source: [('SOHO', 'Solar and Heliospheric Observatory')],
        attrs.Provider: [('SDAC', 'Solar Data Analysis Center')],
        attrs.Detector: [('C1', 'Coronograph 1'),
                        ('C2', 'Coronograph 2'),
                        ('C3', 'Coronograph 3')]
        }

        return adict
    
    def _make_search(self, query):
        raise ConnectionRefusedError                                        

    def _make_filename(path, row, resp, url):
        # Define a fallback filename based on the information in the search results
        name = f"row['ID'].fits"

        if resp:
            cdheader = resp.headers.get("Content-Disposition", None)
            if cdheader:
                _, params = sunpy.util.net.parse_header(cdheader)
                name = params.get('filename', "")

        return path.format(file=name, **row.response_block_map)

    @convert_row_to_table
    def fetch(self, query_results, *, path, downloader, **kwargs):
        for row in query_results:
            filepath = partial(self._make_filename, path, row)

            url = f"https://sfsi.sunpy.org/download/{row['ID']}"
            downloader.enqueue_file(url, filename=filepath)

    @classmethod
    def _can_handle_query(cls, *query):
        query_attrs = set(type(x) for x in query)
        supported_attrs = {a.Time, a.Level}
        return supported_attrs.issuperset(query_attrs)
    
client = ExampleClient()
# client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00')) 
print(client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00')))








