from sunpy.net.dataretriever import GenericClient

__all__ = ['GONGClient']


class GONGClient(GenericClient):
    """
    Provides access to the Magnetogram products of NSO-GONG synoptic Maps.

    Searches data hosted by the `National Solar Observatory <https://gong2.nso.edu/oQR/zqs>`__

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2019/12/31 22:00", "2020/1/1 02:00"),
    ...                       a.Instrument('GONG'))  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    5 Results from the GONGClient:
    Source: https://gong2.nso.edu/oQR/zqs
    <BLANKLINE>
           Start Time               End Time        ... Provider ExtentType
    ----------------------- ----------------------- ... -------- ----------
    2019-12-31 22:14:00.000 2019-12-31 22:14:59.999 ...      NSO   SYNOPTIC
    2019-12-31 23:04:00.000 2019-12-31 23:04:59.999 ...      NSO   SYNOPTIC
    2019-12-31 23:54:00.000 2019-12-31 23:54:59.999 ...      NSO   SYNOPTIC
    2020-01-01 00:14:00.000 2020-01-01 00:14:59.999 ...      NSO   SYNOPTIC
    2020-01-01 01:14:00.000 2020-01-01 01:14:59.999 ...      NSO   SYNOPTIC
    <BLANKLINE>
    <BLANKLINE>
    """
    pattern = ("https://gong2.nso.edu/oQR/zqs/{{year:4d}}{{month:2d}}/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}"
               "/mrzqs{{year:2d}}{{month:2d}}{{day:2d}}t{{hour:2d}}{{minute:2d}}c{{}}.fits.gz")

    @property
    def info_url(self):
        return 'https://gong2.nso.edu/oQR/zqs'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [("GONG", "Global Oscillation Network Group.")],
                 attrs.Physobs: [("LOS_MAGNETIC_FIELD", "Line of sight magnetic field")],
                 attrs.Provider: [('NSO', 'National Solar Observatory.')],
                 attrs.ExtentType: [("SYNOPTIC", "Coverage of a complete solar rotation synthesized over time")]}
        return adict
