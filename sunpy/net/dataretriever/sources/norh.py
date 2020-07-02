#  Author: Rishabh Sharma <rishabh.sharma.gunner@gmail.com>
#  This Module was developed under funding provided by
#  Google Summer of Code 2014

import astropy.units as u

from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient
from sunpy.time import TimeRange

__all__ = ['NoRHClient']

BASEURL = 'ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/%Y/%m/{freq}%y%m%d'


class NoRHClient(GenericClient):
    """
    Provides access to the Nobeyama RadioHeliograph (NoRH) averaged correlation
    time series data.

    Uses this `ftp archive <ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/>`__
    hosted by the `NoRH Science Center <https://solar.nro.nao.ac.jp/norh/doc/manuale/node1.html>`__.

    Queries to NoRH should specify either 17GHz or 34GHz as a Wavelength.

    Examples
    --------
    >>> import astropy.units as u
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/2"),
    ...                       a.Instrument.norh, a.Wavelength(17*u.GHz))  #doctest: +REMOTE_DATA
    >>> results  #doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the NoRHClient:
         Start Time           End Time      Instrument Source Provider Wavelength
    ------------------- ------------------- ---------- ------ -------- ----------
    2016-01-01 00:00:00 2016-01-01 23:59:59       NORH   NAOJ      NRO   17.0 GHz
    2016-01-02 00:00:00 2016-01-02 23:59:59       NORH   NAOJ      NRO   17.0 GHz
    <BLANKLINE>
    <BLANKLINE>

    """
    baseurl = r'ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/%Y/%m/(\w){3}%y%m%d'
    pattern = '{}/tcx/{year:4d}/{month:2d}/{Wavelength:3l}{:4d}{day:2d}'
    optional = {a.Wavelength, a.Source, a.Provider}

    @classmethod
    def pre_hook(cls, *args, **kwargs):
        d = {'Source': ['NAOJ'], 'Instrument': ['NORH'], 'Provider': ['NRO']}
        d['Wavelength'] = []
        waverange = a.Wavelength(34*u.GHz, 17*u.GHz)
        req_wave = waverange
        for elem in args:
            if isinstance(elem, a.Wavelength):
                req_wave = elem
            elif isinstance(elem, a.Time):
                timerange = TimeRange(elem.start, elem.end)
                d['timerange'] = timerange
        if hasattr(kwargs, 'Wavelength'):
            req_wave = kwargs['Wavelength']
        wmin = req_wave.min.to(u.GHz, equivalencies=u.spectral())
        wmax = req_wave.max.to(u.GHz, equivalencies=u.spectral())
        req_wave = a.Wavelength(wmin, wmax)
        if waverange.min in req_wave:
            d['Wavelength'].append('tcz')
        if waverange.max in req_wave:
            d['Wavelength'].append('tca')
        return cls.baseurl, cls.pattern, d

    def post_hook(self, exdict, matchdict):
        map_ = super().post_hook(exdict, matchdict)
        if map_['Wavelength'] == 'tcz':
            map_['Wavelength'] = 17*u.GHz
        elif map_['Wavelength'] == 'tca':
            map_['Wavelength'] = 34*u.GHz
        return map_

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('NORH',
                                     ('Nobeyama Radio Heliograph is an imaging radio telescope at 17 '
                                      'or 34GHz located at the Nobeyama Solar Radio Observatory.'))],
                 attrs.Source: [('NAOJ', 'The National Astronomical Observatory of Japan')],
                 attrs.Provider: [('NRO', 'Nobeyama Radio Observatory')]}
        return adict
