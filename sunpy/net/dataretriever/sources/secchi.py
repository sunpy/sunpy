import socket
from datetime import datetime

from sunpy.extern.parse import parse
from sunpy.net.dataretriever import GenericClient, QueryResponse
from sunpy.util.scraper import Scraper
from sunpy.time import TimeRange, parse_time


__all__ = ['SECCHIlient']


class SECCHIlient(GenericClient):
    """
    Provides access to the STEREO/SECCHI data.

    Uses this `archive <https://stereo-ssc.nascom.nasa.gov/pub/beacon/ahead/secchi/img>`.
    """

    patterns = {
        'euvi': '{}/img/euvi/{year:4d}{month:2d}{day:2d}/{year:4d}{month:2d}{day:2d}_{}',
        'cor2': '{}/img/cor2/{year:4d}{month:2d}{day:2d}/{year:4d}{month:2d}{day:2d}_{}',
        'hi_1': '{}/img/hi_1/{year:4d}{month:2d}{day:2d}/{year:4d}{month:2d}{day:2d}_{}',
        'hi_2': '{}/img/hi_2/{year:4d}{month:2d}{day:2d}/{year:4d}{month:2d}{day:2d}_{}'
    }
    baseurl = 'https://stereo-ssc.nascom.nasa.gov/pub/beacon/ahead/secchi/img'

    def get_filenames(self, timerange):
        telescopes = ['euvi', 'cor2', 'hi_1', 'hi_2']
        tele_url = {
            'euvi': 'n7euA',
            'cor2': 'd7c2A',
            'hi_1': 's7h1A',
            'hi_2': 's7h2A'
        }

        filenames = {}
        for telescope in telescopes:
            pattern = ("""https://stereo-ssc.nascom.nasa.gov/pub/beacon/ahead/
                       {instrument}/img/{telescope}/%Y%m%d/%Y%m%d_%H%M%S_{tele_url}.fts""")
            solmon = Scraper(pattern, instrument='secchi',
                             telescope=telescope, tele_url=tele_url[telescope])

            files = solmon.filelist(timerange)
            filenames[telescope] = files
        return filenames

    def search(self, *args, **kwargs):

        baseurl, patterns, matchdict = self.pre_search_hook(*args, **kwargs)
        timerange = matchdict['Time']
        metalist = []
        filenames = self.get_filenames(timerange)

        for telescope, url_list in filenames.items():
            for url in url_list:
                exdict = parse(self.patterns[telescope], url).named
                exdict['url'] = url
                rowdict = super().post_search_hook(exdict, matchdict)
                metalist.append(rowdict)

        return QueryResponse(metalist, client=self)

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('SECCHI',
                                     """Sun Earth Connection Coronal 
                                     and Heliospheric Investigation.""")],
                 attrs.Source: [('STEREO', 'Solar Terrestrial Relations Observatory.')],
                 attrs.Provider: [('NASA', 'The National Aeronautics and Space Administration.')]}
        return adict
