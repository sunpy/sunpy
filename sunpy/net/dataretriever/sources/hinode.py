from sunpy.net import attrs as a
from sunpy.net.dataretriever.attrs import hinode as hattrs
from sunpy.net.dataretriever.client import GenericClient, QueryResponse
from sunpy.net.scraper import Scraper
from sunpy.time import TimeRange

__all__ = ['HinodeSOTFGClient', 'HinodeSOTSPClient']


_ARCHIVE_URLS = {
    'DARTS': 'https://data.darts.isas.jaxa.jp/pub/hinode/sot/',
    'LMSAL': 'https://sot.lmsal.com/data/sot/',
}


def _get_archive_url(matchdict):
    """
    Defaults to DARTS unless LMSAL is the only provider selected.
    """
    providers = [p.upper() for p in matchdict.get('Provider', [])]
    if providers == ['LMSAL']:
        return _ARCHIVE_URLS['LMSAL']
    return _ARCHIVE_URLS['DARTS']


class HinodeSOTFGClient(GenericClient):
    """
    Provides access to the Hinode Solar Optical Telescope (SOT) Filtergraph (FG) data archive.

    Data is available from two archives:

    - `DARTS archive <https://darts.isas.jaxa.jp/en>`__ (default)
    - `LMSAL archive <https://www.lmsal.com/solarsoft/hinode/doc/hinode_sot_analysis_guide.html>`__

    Use ``a.Provider('LMSAL')`` to retrieve from LMSAL instead.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2013-07-13 00:13:32", "2013-07-13 00:13:33"),
    ...                       a.Instrument.sot, a.Level(0), a.hinode.SOTDetector('FG'))  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the HinodeSOTFGClient:
    Source: https://data.darts.isas.jaxa.jp/pub/hinode/sot/
    <BLANKLINE>
           Start Time               End Time        Instrument Source Provider Level SOTDetector
    ----------------------- ----------------------- ---------- ------ -------- ----- -----------
    2013-07-13 00:13:32.000 2013-07-13 00:13:32.999        SOT HINODE    DARTS     0          FG
    <BLANKLINE>
    <BLANKLINE>
    """
    required = {a.Time, a.Instrument, a.Level, hattrs.SOTDetector}
    optional = {a.Source, a.Provider}
    _pattern = (
        "{archive}level{Level}/{{year:4d}}/{{month:2d}}/{{day:2d}}/FG/"
        "H{{hour:2d}}00/"
        "FG{{year:4d}}{{month:2d}}{{day:2d}}_"
        "{{hour:2d}}{{minute:2d}}{{second:2d}}.{{}}.fits"
    )

    @classmethod
    def _attrs_module(cls):
        return 'hinode', 'sunpy.net.dataretriever.attrs.hinode'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs as a

        return {
            a.Instrument: [("SOT", "Hinode Solar Optical Telescope")],
            a.Source: [("Hinode", "Hinode / Solar-B")],
            a.Provider: [("DARTS", "ISAS/JAXA DARTS archive"),
                         ("LMSAL", "Lockheed Martin Solar and Astrophysics Laboratory")],
            a.Level: [("0", "Level 0"), ("1", "Level 1")],
            a.hinode.SOTDetector: [("FG", "Filtergraph")],
        }

    @property
    def info_url(self):
        return _ARCHIVE_URLS['DARTS']

    def search(self, *args, **kwargs):
        matchdict = self._get_match_dict(*args, **kwargs)
        level = str(matchdict['Level'][0])
        if level not in ['0', '1']:
            return QueryResponse([], client=self)
        archive_url = _get_archive_url(matchdict)
        pat = self._pattern.replace('{archive}', archive_url).replace('{Level}', level)
        tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        rows = []
        scraper = Scraper(format=pat)
        try:
            filemeta = scraper._extract_files_meta(tr, matcher=matchdict)
        except Exception:
            filemeta = []
        for i in filemeta:
            if i.get('url'):
                rows.append(self.post_search_hook(i, matchdict))
        rows.sort(key=lambda x: x['Start Time'])
        return QueryResponse(rows, client=self)


class HinodeSOTSPClient(GenericClient):
    """
    Provides access to the Hinode Solar Optical Telescope (SOT) Spectral-polarimeter (SP) data archive.

    Data is available from two archives:

    - `DARTS archive <https://darts.isas.jaxa.jp/en>`__ (default)
    - `LMSAL archive <https://www.lmsal.com/solarsoft/hinode/doc/hinode_sot_analysis_guide.html>`__

    Use ``a.Provider('LMSAL')`` to retrieve from LMSAL instead.

    .. note::

        Level 1 observations are grouped by observation start time in directories .
        Each directory can contain hundreds of individual scan files spanning a time range longer
        than the observation start time. As the scraper uses the directory timestamp for all files
        within it, all returned level 1 results will share the same observation start time.
        A search must include the observation start time in its time range to find level 1 files.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2013-07-13 10:02:00", "2013-07-13 10:03:00"),
    ...                       a.Instrument.sot, a.Level(0), a.hinode.SOTDetector('SP'))  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the HinodeSOTSPClient:
    Source: https://data.darts.isas.jaxa.jp/pub/hinode/sot/
    <BLANKLINE>
           Start Time               End Time        Instrument Source Provider Level SOTDetector
    ----------------------- ----------------------- ---------- ------ -------- ----- -----------
    2013-07-13 10:02:50.000 2013-07-13 10:02:50.000        SOT Hinode    DARTS     0          SP
    2013-07-13 10:02:54.000 2013-07-13 10:02:54.000        SOT Hinode    DARTS     0          SP
    2013-07-13 10:02:58.000 2013-07-13 10:02:58.000        SOT Hinode    DARTS     0          SP
    <BLANKLINE>
    <BLANKLINE>
    """
    required = {a.Time, a.Instrument, a.Level, hattrs.SOTDetector}
    optional = {a.Source, a.Provider}
    _patterns = {
        '0': (
            "{archive}level0/{{year:4d}}/{{month:2d}}/{{day:2d}}/SP4D/H{{hour:2d}}00/"
            "SP4D{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}.{{}}.fits"
        ),
        '1': (
            "{archive}level1hao/{{year:4d}}/{{month:2d}}/{{day:2d}}/SP3D/"
            "{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}/"
            "SP3D{{}}.fits"
        ),
        '2': (
            "{archive}level2hao/{{year:4d}}/{{month:2d}}/{{day:2d}}/SP3D/"
            "{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}/"
            "{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}.fits"
        ),
        '2.1': (
            "{archive}level2.1hao/{{year:4d}}/{{month:2d}}/{{day:2d}}/SP3D/"
            "{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}/"
            "{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}_L2.1.fits"
        ),
    }

    @classmethod
    def _attrs_module(cls):
        return 'hinode', 'sunpy.net.dataretriever.attrs.hinode'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs as a

        return {
            a.Instrument: [("SOT", "Hinode Solar Optical Telescope")],
            a.Source: [("Hinode", "Hinode / Solar-B")],
            a.Provider: [("DARTS", "ISAS/JAXA DARTS archive"),
                         ("LMSAL", "Lockheed Martin Solar and Astrophysics Laboratory")],
            a.Level: [("0", "Level 0"), ("1", "Level 1"), ("2", "Level 2"), ("2.1", "Level 2.1")],
            a.hinode.SOTDetector: [("SP", "Spectro-Polarimeter")],
        }

    @property
    def info_url(self):
        return _ARCHIVE_URLS['DARTS']

    def post_search_hook(self, exdict, matchdict):
        rowdict = super().post_search_hook(exdict, matchdict)
        rowdict['End Time'] = rowdict['Start Time']
        rowdict['Source'] = 'Hinode'
        providers = [p.upper() for p in matchdict.get('Provider', [])]
        rowdict['Provider'] = 'LMSAL' if providers == ['LMSAL'] else 'DARTS'
        return rowdict

    def search(self, *args, **kwargs):
        matchdict = self._get_match_dict(*args, **kwargs)
        level = str(matchdict['Level'][0])
        pattern_template = self._patterns.get(level)
        if pattern_template is None:
            return QueryResponse([], client=self)
        archive_url = _get_archive_url(matchdict)
        pattern = pattern_template.replace('{archive}', archive_url)
        tr = TimeRange(matchdict['Start Time'], matchdict['End Time'])
        rows = []
        scraper = Scraper(format=pattern)
        try:
            filemeta = scraper._extract_files_meta(tr, matcher=matchdict)
        except Exception:
            filemeta = []
        for i in filemeta:
            if i.get('url'):
                rows.append(self.post_search_hook(i, matchdict))
        rows.sort(key=lambda x: x['Start Time'])
        return QueryResponse(rows, client=self)
