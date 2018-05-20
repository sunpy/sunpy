"""PROBA2 Map subclass definitions"""
from __future__ import absolute_import, print_function, division
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import matplotlib.pyplot as plt

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
    * `Proba2 SWAP Science Center <http://proba2.sidc.be/about/SWAP/>`_
    * `Fits headers reference <http://proba2.oma.be/data/SWAP/level0>`_
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        # It needs to be verified that these must actually be set and
        # are not already in the header.
        self.meta['detector'] = "SWAP"
#        self.meta['instrme'] = "SWAP"
        self.meta['obsrvtry'] = "PROBA2"

        self._nickname = self.detector
        self.plot_settings['cmap'] = plt.get_cmap(name='sdoaia171')

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an SWAP image"""
        return header.get('instrume') == 'SWAP'
