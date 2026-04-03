from sunpy.net import attrs as a
from sunpy.net.dataretriever.client import GenericClient, QueryResponse
from sunpy.net.scraper import Scraper
from sunpy.net.scraper_utils import get_timerange_from_exdict
from sunpy.time import TimeRange

__all__ = ['HinodeSOTFGClient', 'HinodeSOTSPClient']

_BASEURL = 'https://data.darts.isas.jaxa.jp/pub/hinode/sot/'
_DAY = '{{year:4d}}/{{month:2d}}/{{day:2d}}/'
_TIME = '{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}{{second:2d}}'


class HinodeSOTFGClient(GenericClient):
    """
    Provides access to the Hinode Solar Optical Telescope (SOT) Filtergraph (FG) data archive.

    Searches level 0 and 1 data hosted by `DARTS <https://darts.isas.jaxa.jp/en>`__.
    Only the plain filtergrams in the archive's ``FG`` folders are searched.
    Add ``a.Provider.darts`` to a search to leave out VSO results.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2011-12-13 06:04:30", "2011-12-13 06:04:40"), a.Instrument.sot,
    ...                       a.Level(1), a.Detector.fg, a.Provider.darts)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    4 Results from the HinodeSOTFGClient:
    Source: https://data.darts.isas.jaxa.jp/pub/hinode/sot/
    <BLANKLINE>
           Start Time               End Time        Instrument Source Provider Level Detector
    ----------------------- ----------------------- ---------- ------ -------- ----- --------
    2011-12-13 06:04:30.000 2011-12-13 06:04:30.999        SOT HINODE    DARTS     1       FG
    2011-12-13 06:04:33.000 2011-12-13 06:04:33.999        SOT HINODE    DARTS     1       FG
    2011-12-13 06:04:36.000 2011-12-13 06:04:36.999        SOT HINODE    DARTS     1       FG
    2011-12-13 06:04:39.000 2011-12-13 06:04:39.999        SOT HINODE    DARTS     1       FG
    <BLANKLINE>
    <BLANKLINE>
    """
    required = {a.Time, a.Instrument, a.Level, a.Detector}
    pattern = _BASEURL + 'level{Level}/' + _DAY + 'FG/H{{hour:2d}}00/FG' + _TIME + '.{{}}.fits'

    @classmethod
    def pre_search_hook(cls, *args, **kwargs):
        baseurl, pattern, matchdict = super().pre_search_hook(*args, **kwargs)
        # The scraper can only fill in times, so the level goes in here.
        return baseurl, pattern.replace('{Level}', matchdict['Level'][0]), matchdict

    @classmethod
    def register_values(cls):
        return {
            a.Instrument: [("SOT", "Hinode Solar Optical Telescope")],
            a.Source: [("Hinode", "Hinode / Solar-B")],
            a.Provider: [("DARTS", "ISAS/JAXA DARTS archive")],
            a.Level: [("0", "Level 0"), ("1", "Level 1")],
            a.Detector: [("FG", "Filtergraph")],
        }

    @property
    def info_url(self):
        return _BASEURL


class HinodeSOTSPClient(GenericClient):
    """
    Provides access to the Hinode Solar Optical Telescope (SOT) Spectro-Polarimeter (SP) data archive.

    Searches level 0, 1, 2 and 2.1 data hosted by `DARTS <https://darts.isas.jaxa.jp/en>`__.
    Add ``a.Provider.darts`` to a search to leave out VSO results.

    .. note::

        Level 1, 2 and 2.1 data are stored per scan.
        A search returns every file of each scan that **started** within the time range:
        hundreds of files per scan for level 1, one for levels 2 and 2.1.
        Scans that started earlier are not returned.
        Each file is listed with the start time of its scan.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2011-12-13", "2011-12-14"), a.Instrument.sot,
    ...                       a.Level(2), a.Detector.sp, a.Provider.darts)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the HinodeSOTSPClient:
    Source: https://data.darts.isas.jaxa.jp/pub/hinode/sot/
    <BLANKLINE>
           Start Time               End Time        Instrument Source Provider Level Detector
    ----------------------- ----------------------- ---------- ------ -------- ----- --------
    2011-12-13 02:50:04.000 2011-12-13 02:50:04.999        SOT HINODE    DARTS     2       SP
    2011-12-13 06:15:04.000 2011-12-13 06:15:04.999        SOT HINODE    DARTS     2       SP
    <BLANKLINE>
    <BLANKLINE>
    """
    required = {a.Time, a.Instrument, a.Level, a.Detector}
    _patterns = {
        '0': _BASEURL + 'level0/' + _DAY + 'SP4D/H{{hour:2d}}00/SP4D' + _TIME + '.{{}}.fits',
        '1': _BASEURL + 'level1hao/' + _DAY + 'SP3D/' + _TIME + '/SP3D{{}}.fits',
        '2': _BASEURL + 'level2hao/' + _DAY + 'SP3D/' + _TIME + '/' + _TIME + '.fits',
        '2.1': _BASEURL + 'level2.1hao/' + _DAY + 'SP3D/' + _TIME + '/' + _TIME + '_L2.1.fits',
    }

    @classmethod
    def register_values(cls):
        return {
            a.Instrument: [("SOT", "Hinode Solar Optical Telescope")],
            a.Source: [("Hinode", "Hinode / Solar-B")],
            a.Provider: [("DARTS", "ISAS/JAXA DARTS archive")],
            a.Level: [("0", "Level 0"), ("1", "Level 1"), ("2", "Level 2"), ("2.1", "Level 2.1")],
            a.Detector: [("SP", "Spectro-Polarimeter")],
        }

    @property
    def info_url(self):
        return _BASEURL

    def search(self, *args, **kwargs):
        matchdict = self._get_match_dict(*args, **kwargs)
        pattern = self._patterns[matchdict['Level'][0]]
        timeranges = [TimeRange(matchdict['Start Time'], matchdict['End Time'])]
        if matchdict['Level'][0] != '0':
            # Scan directories are named to the second, so letting the scraper guess them costs one request
            # per second searched. List each day's scans instead; the {} matches the "/" ending a directory link.
            scans = Scraper(format=pattern.rsplit('/', 1)[0] + '{{}}')._extract_files_meta(timeranges[0])
            timeranges = [get_timerange_from_exdict(scan) for scan in scans]
        scraper = Scraper(format=pattern)
        rows = [self.post_search_hook(meta, matchdict) for tr in timeranges for meta in scraper._extract_files_meta(tr)]
        return QueryResponse(sorted(rows, key=lambda row: row['url']), client=self)
