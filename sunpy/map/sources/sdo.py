"""SDO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np
from matplotlib import colors

from sunpy.map import GenericMap
from sunpy.cm import cm

__all__ = ['AIAMap', 'HMIMap']

class AIAMap(GenericMap):
    """AIA Image Map definition

    References
    ----------
    For a description of AIA headers
    http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_A_AIA-SDO_FITS_Keyword_Documents.pdf
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        # Fill in some missing info
        self.meta['detector'] = "AIA"
#        self.meta['instrme'] = "AIA"

        self._nickname = self.detector
        self.plot_settings['cmap'] = cm.get_cmap(self._get_cmap_name())

        norm = {304: colors.PowerNorm(0.35), 211: colors.PowerNorm(0.5),
                193: colors.PowerNorm(0.5), 171: colors.PowerNorm(0.5),
                131: colors.PowerNorm(0.5), 335: colors.PowerNorm(0.4),
                94: colors.PowerNorm(0.45)}
        self.plot_settings['norm'] = norm[int(self.wavelength.to('Angstrom').value)]

    @property
    def observatory(self):
        return self.meta['telescop'].split('/')[0]

    @property
    def processing_level(self):
        return self.meta['lvl_num']

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume', '').startswith('AIA')

class HMIMap(GenericMap):
    """HMI Image Map definition"""

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        self.meta['detector'] = "HMI"
#        self.meta['instrme'] = "HMI"
#        self.meta['obsrvtry'] = "SDO"
        self.meta['waveunit'] = 'Angstrom'
        self._nickname = self.detector

    @property
    def measurement(self):
        return self.meta['content'].split(" ")[0].lower()

    @property
    def observatory(self):
        return self.meta['telescop'].split('/')[0]

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an HMI image"""
        return header.get('instrume', '').startswith('HMI')
