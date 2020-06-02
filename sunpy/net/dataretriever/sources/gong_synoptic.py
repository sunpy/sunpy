from sunpy.net.dataretriever import GenericClient
from sunpy.util.scraper import Scraper

__all__ = ['GongSynopticClient']


class GongSynopticClient(GenericClient):
    """
    Provides access to the Magnetogram products of NSO-GONG synoptic Maps.

    Searches data hosted by the `National Solar Observatory <gong2.nso.edu/oQR/zqs/>`__

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2019/12/31 22:00", "2020/1/1 02:00"),
    ...                       a.Instrument('GONG'))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the GongSynopticClient:
         Start Time           End Time      Source Instrument Wavelength
    ------------------- ------------------- ------ ---------- ----------
    2019-12-31 22:00:00 2020-01-01 02:00:00    NSO       GONG        nan
    2019-12-31 22:00:00 2020-01-01 02:00:00    NSO       GONG        nan
    2019-12-31 22:00:00 2020-01-01 02:00:00    NSO       GONG        nan
    <BLANKLINE>
    <BLANKLINE>

    """

    def _get_url_for_timerange(self, timerange, **kwargs):
        """
        Return URL(s) for corresponding timerange.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The time range required.

        Returns
        -------
        `list`
            The URL(s) for the corresponding timerange.
        """
        prefix = r'https://gong2.nso.edu/oQR/zqs/%Y%m/mrzqs%y%m%d/'
        pattern = prefix + r'mrzqs%y%m%dt%H%Mc(\d){4}_(\d){3}\.fits.gz'
        gongsynoptic_files = Scraper(pattern, regex=True)
        urls = gongsynoptic_files.filelist(timerange)
        return urls

    def _makeimap(self):
        """
        Helper Function:used to hold information about source.
        """
        self.map_['instrument'] = 'GONG'
        self.map_['source'] = 'NSO'

    @classmethod
    def _can_handle_query(cls, *query):
        """
        Answers whether client can service the query.

        Parameters
        ----------
        query : list of query objects

        Returns
        -------
        bool : Answer as to whether client can service the query.
        """
        from sunpy.net import attrs as a

        required = {a.Time, a.Instrument}
        optional = {a.Physobs, a.gong_synoptic.ExtentType}
        if not cls.check_attr_types_in_query(query, required, optional):
            return False

        matches = True
        for x in query:
            if isinstance(x, a.Instrument) and x.value.lower() != 'gong':
                matches = False
            if isinstance(x, a.Physobs) and x.value.lower() != 'los_magnetic_field':
                matches = False
            if isinstance(x, a.gong_synoptic.ExtentType) and x.value.lower() != 'synoptic':
                matches = False

        return matches

    @classmethod
    def _attrs_module(cls):
        return 'gong_synoptic', 'sunpy.net.dataretriever.attrs.gong_synoptic'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [
            ("GONG", "Global Oscillation Network Group.")],
            attrs.gong_synoptic.ExtentType: [
            ("SYNOPTIC", "Coverage of a complete solar rotation synthesized over time")],
            attrs.Physobs: [("LOS_MAGNETIC_FIELD", "Line of sight magnetic field")]}
        return adict
