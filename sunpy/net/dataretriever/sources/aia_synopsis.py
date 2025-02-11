from sunpy.net import attrs as a
from sunpy.net.dataretriever import GenericClient

__all__ = ["AIASynopsisClient"]


class AIASynopsisClient(GenericClient):
    """
    A client for retrieving AIA "synoptic" data from the JSOC.

    This dataset is not synoptic like HMI and MDI Synoptic images which are images of the solar surface reconstructed from many observations over a solar rotation but rather a synopsis of AIA data.

    The AIA synoptic data are calibrated Level 1.5 images with reduced 1k x 1k resolution at regular 2-minute cadence.

    References
    ----------
    * `Readme <https://jsoc1.stanford.edu/data/aia/synoptic/README.html>`__

    Examples
    --------
    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2016/1/1", "2016/1/1 00:01:00"),
    ...                       a.Instrument.aia, a.Level("1.5s"))  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    9 Results from the AIASynopsisClient:
    Source: https://jsoc1.stanford.edu/data/aia/synoptic/
    <BLANKLINE>
           Start Time               End Time        Instrument  Physobs  Source Provider Level wavelength
    ----------------------- ----------------------- ---------- --------- ------ -------- ----- ----------
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S         94
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        131
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        171
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        193
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        211
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        304
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S        335
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S       1600
    2016-01-01 00:00:00.000 2016-01-01 00:00:59.999        AIA intensity    SDO     JSOC  1.5S       4500
    <BLANKLINE>
    <BLANKLINE>
    """
    pattern = ('https://jsoc1.stanford.edu/data/aia/synoptic/'
               '{{year:4d}}/{{month:2d}}/{{day:2d}}/H{{hour:2d}}00/AIA{{year:4d}}{{month:2d}}{{day:2d}}_{{hour:2d}}{{minute:2d}}_{{wavelength:04d}}.fits')
    known_wavelengths = [94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500]
    # required = {a.Time, a.Instrument, a.Level}

    @property
    def info_url(self):
        return "https://jsoc1.stanford.edu/data/aia/synoptic/"

    @classmethod
    def register_values(cls):
        adict = {
            a.Instrument: [("AIA", "Data from the Atmospheric Imaging Assembly instrument.")],
            a.Physobs: [
                ("intensity", "Brightness or intensity of the solar atmosphere at different wavelengths.")
            ],
            a.Source: [("SDO", "The Solar Dynamics Observatory.")],
            a.Wavelength: [(f"{wv:04d}", f"{wv} Å") for wv in cls.known_wavelengths],
            a.Provider: [("JSOC", "Joint Science Operations Center at Stanford.")],
            a.Level: [("1.5s", "Level 1.5 data processed for quicker analysis.")],
        }
        return adict
