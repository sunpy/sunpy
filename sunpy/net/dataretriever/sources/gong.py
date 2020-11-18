from sunpy.net.dataretriever import GenericClient

__all__ = ['GONGClient']


class GONGClient(GenericClient):
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
    3 Results from the GONGClient:
         Start Time           End Time      Instrument ... Source ExtentType
    ------------------- ------------------- ---------- ... ------ ----------
    2019-12-31 22:14:00 2019-12-31 22:14:59       GONG ...    NSO   SYNOPTIC
    2019-12-31 23:04:00 2019-12-31 23:04:59       GONG ...    NSO   SYNOPTIC
    2019-12-31 23:54:00 2019-12-31 23:54:59       GONG ...    NSO   SYNOPTIC
    <BLANKLINE>
    <BLANKLINE>

    """
    baseurl = (r'https://gong2.nso.edu/oQR/zqs/%Y%m/mrzqs%y%m%d/mrzqs%y%m%dt%H%Mc'
               r'(\d){4}_(\d){3}\.fits.gz')
    pattern = '{}/zqs/{year:4d}{month:2d}/mrzqs{:4d}{day:2d}/mrzqs{:6d}t{hour:2d}{minute:2d}c{}'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [("GONG", "Global Oscillation Network Group.")],
                 attrs.Physobs: [("LOS_MAGNETIC_FIELD", "Line of sight magnetic field")],
                 attrs.Source: [('NSO', 'National Solar Observatory.')],
                 attrs.ExtentType: [("SYNOPTIC", "Coverage of a complete solar rotation synthesized over time")]}
        return adict
