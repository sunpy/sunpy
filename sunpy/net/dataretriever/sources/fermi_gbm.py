from sunpy.net.dataretriever import GenericClient

__all__ = ['GBMClient']


class GBMClient(GenericClient):
    """
    Provides access to data from the Gamma-Ray Burst Monitor (GBM) instrument
    on board the Fermi satellite.

    Although GBMs primary objective is to detect gamma-ray bursts,
    it provides high quality high energy solar flare observations.

    The instrument consists of 12 Sodium Iodide (NaI) scintillation
    detectors, which are sensitive to an energy range of 4keV to 1MeV.
    At any one time, 6 of the NaI detectors are Sunward facing.
    The detectors are numbered 'n1' to 'n11'. This client supports the user
    to choose which detector to use through the `a.Detector <sunpy.net.attrs.Detector>` attribute.
    The default detector is 'n5'.

    The GBM data comes in daily version files in two formats:

        * CSPEC - counts accumulated every  4.096 seconds in 128 energy channels for each detector.
        * CTIME - counts accumulated every 0.256 seconds in 8 energy channels

    Both of which can be accessed through the attrs `a.Resolution <sunpy.net.attrs.Resolution>`.
    The default data type is CSPEC unless the user defines.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> res = Fido.search(a.Time('2015-06-21 00:00', '2015-06-23 23:59'),
    ...                   a.Instrument.gbm, a.Detector.n3,
    ...                   a.Resolution.ctime) # doctest: +REMOTE_DATA
    >>> res # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the GBMClient:
    Source: https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily
    <BLANKLINE>
           Start Time               End Time        ... Resolution Detector
    ----------------------- ----------------------- ... ---------- --------
    2015-06-21 00:00:00.000 2015-06-21 23:59:59.999 ...      ctime       n3
    2015-06-22 00:00:00.000 2015-06-22 23:59:59.999 ...      ctime       n3
    2015-06-23 00:00:00.000 2015-06-23 23:59:59.999 ...      ctime       n3
    <BLANKLINE>
    <BLANKLINE>

    """

    pattern = ('https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/{{year:4d}}/{{month:2d}}/{{day:2d}}/current/'
                'glg_{{Resolution:5}}_{{Detector:2}}_{{year:2d}}{{month:2d}}{{day:2d}}_v00.pha')

    @property
    def info_url(self):
        return 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily'

    @classmethod
    def register_values(cls):
        from sunpy.net import attrs
        adict = {attrs.Instrument: [('GBM', 'Gamma-Ray Burst Monitor on board the Fermi satellite.')],
                 attrs.Physobs: [('flux', 'a measure of the amount of radiation received by an object from a given source.')],
                 attrs.Source: [('FERMI', 'The Fermi Gamma-ray Space Telescope.')],
                 attrs.Provider: [('NASA', 'The National Aeronautics and Space Administration.')],
                 attrs.Resolution: [
            ("cspec", "CSPEC 128 channel spectra every 4.096 seconds."),
            ("ctime", "CTIME provides 8 channel spectra every 0.256 seconds.")],
            attrs.Detector: [(f"n{x}", f"GBM Detector short name for the detector NAI_{x:02}") for x in range(12)]}
        return adict
