"""STEREO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import astropy.units as u
import numpy as np

from sunpy.map import GenericMap
from sunpy.cm import cm
from matplotlib.colors import PowerNorm

__all__ = ['EUVIMap', 'CORMap']

class EUVIMap(GenericMap):
    """EUVI Image Map definition"""

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        self._nickname = "{0}-{1}".format(self.detector, self.observatory[-1])
        self.plot_settings['cmap'] = cm.get_cmap('sohoeit{wl:d}'.format(wl=int(self.wavelength.value)))
        self.meta['waveunit'] = 'Angstrom'

        # Try to identify when the FITS meta data does not have the correct
        # date FITS keyword
        if ('date_obs' in self.meta) and not('date-obs' in self.meta):
            self.meta['date-obs'] = self.meta['date_obs']

    @property
    def rsun_arcseconds(self):
        """
        Radius of the sun in arcseconds.

        References
        ----------
        http://sohowww.nascom.nasa.gov/solarsoft/stereo/secchi/doc/FITS_keywords.pdf
        """
        return self.meta.get('rsun', None)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EUVI image"""
        return header.get('detector') == 'EUVI'

    def _get_mpl_normalizer(self):
        """Returns a Normalize object to be used with XRT data"""
        return PowerNorm(0.25, self.data.min(), self.data.max())


class CORMap(GenericMap):
    """COR Image Map definition"""

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        self._nickname = "{0}-{1}".format(self.detector, self.observatory[-1])
        self.plot_settings['cmap'] = cm.get_cmap('stereocor{det!s}'.format(det=self.detector[-1]))

        # Try to identify when the FITS meta data does not have the correct
        # date FITS keyword
        if ('date_obs' in self.meta) and not('date-obs' in self.meta):
            self.meta['date-obs'] = self.meta['date_obs']

    @property
    def measurement(self):
        # TODO: This needs to do more than white-light.  Should give B, pB, etc.
        return "white-light"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an COR image"""
        return header.get('detector', '').startswith('COR')

    def _get_mpl_normalizer(self):
        """Returns a Normalize object to be used with XRT data"""
        return PowerNorm(0.25, self.data.min(), self.data.max())

class HIMap(GenericMap):
    """HI Image Map definition"""

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        self._nickname = "{0}-{1}".format(self.detector, self.observatory[-1])
        self.plot_settings['cmap'] = cm.get_cmap('stereohi{det!s}'.format(det=self.detector[-1]))

        # Try to identify when the FITS meta data does not have the correct
        # date FITS keyword
        if ('date_obs' in self.meta) and not('date-obs' in self.meta):
            self.meta['date-obs'] = self.meta['date_obs']

    @property
    def measurement(self):
        # TODO: This needs to do more than white-light.  Should give B, pB, etc.
        return "white-light"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an COR image"""
        return header.get('detector', '').startswith('HI')

    def _get_mpl_normalizer(self):
        """Returns a Normalize object to be used with XRT data"""
        return PowerNorm(0.25, self.data.min(), self.data.max())
