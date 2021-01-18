import socket
from datetime import datetime

from sunpy.extern.parse import parse
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange, parse_time

__all__ = ['SECCHIBEACONlient']


class SECCHIBEACONlient(GenericClient):
    """
    Provides access to the STEREO/SECCHI data.

    Uses this `archive <https://stereo-ssc.nascom.nasa.gov/pub/beacon/>`.
    """

    pattern = "{}/{satellite}/secchi/img/{instrument}/"
                 "{year:4d}{month:2d}{day:2d}/{year:4d}{month:2d}{day:2d}_{}"
    baseurl = 'https://stereo-ssc.nascom.nasa.gov/pub/beacon/'

    def get_filenames(self, timerange, telescope, satellite):
        tele_url = {
            'euvi': 'n7eu'+satellite[0].upper(),
            'cor2': 'd7c2'+satellite[0].upper(),
            'hi_1': 's7h1'+satellite[0].upper(),
            'hi_2': 's7h2'+satellite[0].upper()
        }

        pattern = ("https://stereo-ssc.nascom.nasa.gov/pub/beacon/{satellite}/"
                      {instrument}/img/{telescope}/%Y%m%d/%Y%m%d_%H%M%S_{tele_url}.fts""")
        secchi_scraper = Scraper(pattern, instrument='secchi',
                         satellite=satellite, telescope=telescope, tele_url=tele_url[telescope])
        files = secchi_scraper.filelist(timerange)
        return files

    def search(self, *args, **kwargs):

        baseurl, pattern, matchdict = self.pre_search_hook(*args, **kwargs)
        timerange = matchdict['Time']
        satellite = matchdict['Satellite']
        telescope = matchdict['Telescope']
        metalist = []

        filenames = self.get_filenames(timerange, telescope[0], satellite[0])

        for url in filenames:
            exdict = parse(pattern, url).named
            exdict['url'] = url
            rowdict = super().post_search_hook(exdict, matchdict)
            metalist.append(rowdict)

        return QueryResponse(metalist, client=self)

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('SECCHI_BEACON',
                                     """Sun Earth Connection Coronal
                                        and Heliospheric Investigation.""")],
                 attrs.Telescope: [('EUVI', 'Extreme UltraViolet Imager'),
                                   ('COR2', 'Outer Coronagraph (NRL)'),
                                   ('HI_1', 'Heliospheric Imager 1'),
                                   ('HI_2', 'Heliospheric Imager 2')],
                 attrs.Satellite: [('ahead', 'Ahead'),
                                   ('behind', 'Behind')],
                 attrs.Source: [('STEREO', 'Solar Terrestrial Relations Observatory.')],
                 attrs.Provider: [('NASA', 'The National Aeronautics and Space Administration.')]}
        return adict
