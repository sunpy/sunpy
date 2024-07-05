#  Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#  This Module was developed under funding provided by
#  Google Summer of Code 2014

import astropy.units as u

from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient

__all__ = ['NoRHClient']


class NoRHClient(GenericClient):
    """
    Provides access to the Nobeyama RadioHeliograph (NoRH) averaged correlation
    time series data.

    Uses this `https archive <https://solar.nro.nao.ac.jp/norh/data/tcx>`__
    hosted by the `NoRH Science Center <https://solar.nro.nao.ac.jp/norh/index.html>`__.

    Queries to NoRH should specify either 17GHz or 34GHz as a Wavelength.

    Examples
    --------
    >>> import astropy.units as u
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.norh, a.Wavelength(17*u.GHz))  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the NoRHClient:
    Source: https://solar.nro.nao.ac.jp/norh/index.html
    <BLANKLINE>
           Start Time               End Time        Instrument Source Provider Wavelength
                                                                                  GHz
    ----------------------- ----------------------- ---------- ------ -------- ----------
    2016-01-01 00:00:00.000 2016-01-01 23:59:59.999       NORH   NAOJ      NRO       17.0
    2016-01-02 00:00:00.000 2016-01-02 23:59:59.999       NORH   NAOJ      NRO       17.0
    <BLANKLINE>
    <BLANKLINE>

    References
    ----------
    * `NoRH Science Center <https://solar.nro.nao.ac.jp/norh/index.html>`__
    * `NoRH Data Archive <https://solar.nro.nao.ac.jp/norh/data/tcx>`__
    * `NoRH User Guide <https://solar.nro.nao.ac.jp/norh/doc/manuale/>`__
    """
    pattern = r'https://solar.nro.nao.ac.jp/norh/data/tcx/{{year:4d}}/{{month:2d}}/{{Wavelength:3l}}{{:4d}}{{day:2d}}'

    @property
    def info_url(self):
        return 'https://solar.nro.nao.ac.jp/norh/index.html'

    @classmethod
    def pre_search_hook(cls, *args, **kwargs):
        """
        Converts the wavelength specified in the query to its
        representation in the url which can be used by the scraper.
        """
        d = cls._get_match_dict(*args, **kwargs)
        waverange = a.Wavelength(34*u.GHz, 17*u.GHz)
        req_wave = d.get('Wavelength', waverange)
        wmin = req_wave.min.to(u.GHz, equivalencies=u.spectral())
        wmax = req_wave.max.to(u.GHz, equivalencies=u.spectral())
        req_wave = a.Wavelength(wmin, wmax)
        d['Wavelength'] = []
        if 17*u.GHz in req_wave:
            d['Wavelength'].append('tca')
        if 34*u.GHz in req_wave:
            d['Wavelength'].append('tcz')
        return cls.baseurl, cls.pattern, d

    def post_search_hook(self, exdict, matchdict):
        """
        This method converts 'tca' and 'tcz' in the url's metadata
        to a frequency of '17 GHz' and '34 GHz' respectively.
        """
        rowdict = super().post_search_hook(exdict, matchdict)
        if rowdict['Wavelength'] == 'tca':
            rowdict['Wavelength'] = 17*u.GHz
        elif rowdict['Wavelength'] == 'tcz':
            rowdict['Wavelength'] = 34*u.GHz
        return rowdict

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('NORH',
                                     ('Nobeyama Radio Heliograph is an imaging radio telescope at 17 '
                                      'or 34GHz located at the Nobeyama Solar Radio Observatory.'))],
                 attrs.Source: [('NAOJ', 'The National Astronomical Observatory of Japan')],
                 attrs.Provider: [('NRO', 'Nobeyama Radio Observatory')],
                 attrs.Wavelength: [('*')]}
        return adict
