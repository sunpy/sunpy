from functools import partial

import sunpy.net
import sunpy.net.attrs as a
import sunpy.util.net
from sunpy.net.attr import AttrAnd, AttrOr, AttrWalker, DataAttr, and_
from sunpy.net.base_client import BaseClient, QueryResponseTable, convert_row_to_table

walker = AttrWalker()

@walker.add_creator(AttrOr)
def create_or(wlk, tree):
    results = []
    for sub in tree.attrs:
        results.extend(wlk.create(sub))
    return results

@walker.add_creator(AttrAnd)
def create_and(wlk, tree):
    params = {}
    for x in tree.attrs:
        wlk.apply(x, params)
    return [params]

@walker.add_applier(DataAttr)
def _(wlk, attr, params):
    # Generic handler for any DataAttr
    return params.update({attr.__class__.__name__.lower(): attr.value})

@walker.add_applier(a.Time)
def _(wlk, attr, params):
    return params.update({
        'start_time': attr.start.isot,
        'end_time': attr.end.isot
    })

@walker.add_applier(a.Instrument)
def _(wlk, attr, params):
    if attr.value.lower() == 'lasco':
        return params.update({'instrument': 'LASCO'})
    return params.update({'instrument': attr.value})


class ExampleClient(BaseClient):
    """
    This is a fake client for testing purposes
    """
    size_column = 'Filesize'

    def search(self, *query):
        query = and_(*query)
        queries = walker.create(query)
        results = []
        for query_parameters in queries:
            results.append(self._make_search(query_parameters))

            # try:
            #     resp = self._make_search(query_parameters)
            #     results.append(resp)
            # except Exception as e:
            #     print("prvented error")
            #     results.append(e)

        return QueryResponseTable(results, client=self)

    @property
    def info_url(self):
        return 'https://cdaweb.gsfc.nasa.gov/index.html'

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
        name = "row['ID'].fits"

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
        supported_attrs = {a.Time, a.Instrument}
        return supported_attrs.issuperset(query_attrs)

# client = ExampleClient()
# client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00'))
# print(client.search(a.Time('2014-01-01T00:00:00', '2014-01-01T00:10:00')))
