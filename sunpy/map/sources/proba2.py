"""PROBA2 Map subclass definitions"""

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map import GenericMap

__all__ = ['SWAPMap']


class SWAPMap(GenericMap):
    """PROBA2 SWAP Image Map.

    The Sun Watcher using Active Pixel System detector and Image Processing (SWAP)
    SWAP provides images of the solar corona at about 17.4 nm, a bandpass
    that corresponds to a temperature of roughly 1 million degrees,
    with a cadence of 1 image per 1-2 minutes, and field of view (FOV) of 54 arcmin.
    It is derived from the SOHO EIT telescope concept design.

    PROBA2 was launched on 2 November 2009.

    References
    ----------
    * `Proba2 SWAP Science Center <http://proba2.sidc.be/about/SWAP/>`__
    * `Fits headers reference <http://proba2.oma.be/data/SWAP/level0>`__
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self._nickname = self.detector
        self.plot_settings['cmap'] = 'sdoaia171'

    @property
    def observatory(self):
        return "PROBA2"

    @property
    def detector(self):
        return "SWAP"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an SWAP image"""
        return header.get('instrume') == 'SWAP'
