"""STEREO Map subclass definitions"""
from __future__ import absolute_import, print_function, division
#pylint: disable=W0221,W0222,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

import numpy as np

from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map import GenericMap
from sunpy.cm import cm


__all__ = ['EUVIMap', 'CORMap', 'HIMap']

class EUVIMap(GenericMap):
    """STEREO-SECCHI EUVI Image Map

    EUVI is an extreme ultraviolet (EUV) imager. Part of the STEREO-SECCHI
    suite it observes the Sun from 1 to 1.7 solar radii. It is capable of
    observing at 304 (He II), 171 (Fe IX), 195 (Fe XII), and 284 (Fe XV)
    Angstroms.

    References
    ----------
    * `STEREO Mission Page <http://stereo.gsfc.nasa.gov>`_
    * `STEREO SECCHI <http://secchi.nrl.navy.mil>`_
    * `Instrument Page <http://secchi.lmsal.com/EUVI/>`_
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        self._nickname = "{0}-{1}".format(self.detector, self.observatory[-1])
        self.plot_settings['cmap'] = cm.get_cmap('sohoeit{wl:d}'.format(wl=int(self.wavelength.value)))
        self.plot_settings['norm'] = ImageNormalize(stretch=PowerStretch(0.25))
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


class CORMap(GenericMap):
    """STEREO-SECCHI CORonograph Image Map.

    Part of the STEREO-SECCHI suite of remote sensing telescopes,
    COR is a set of two coronographs (COR1, COR2) onboard STEREO.
    They are both traditional Lyot coronagraphs.

    The COR1 detectors observes from 1.3 to 4 solar radii while the
    COR2 detectors observe a range from 2 to 15 solar radii.

    References
    ----------
    * `STEREO Mission Page <http://stereo.gsfc.nasa.gov>`_
    * `STEREO SECCHI <http://secchi.nrl.navy.mil>`_
    * `COR1 Instrument Page <http://cor1.gsfc.nasa.gov>`_
    * `COR2 Instrument Page <http://secchi.nrl.navy.mil/index.php?p=cor2>`_
    * `COR1 User Guide <http://cor1.gsfc.nasa.gov/guide/>`_
    """

    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)

        self._nickname = "{0}-{1}".format(self.detector, self.observatory[-1])
        self.meta['wavelnth'] = np.nan
        self.meta['waveunit'] = 'nm'
        self.plot_settings['cmap'] = cm.get_cmap('stereocor{det!s}'.format(det=self.detector[-1]))
        self.plot_settings['norm'] = ImageNormalize(stretch=PowerStretch(0.5))

        # Try to identify when the FITS meta data does not have the correct
        # date FITS keyword
        if ('date_obs' in self.meta) and not('date-obs' in self.meta):
            self.meta['date-obs'] = self.meta['date_obs']

    @property
    def measurement(self):
        """
        Returns the type of data observed.
        """
        # TODO: This needs to do more than white-light.  Should give B, pB, etc.
        return "white-light"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an COR image"""
        return header.get('detector', '').startswith('COR')


class HIMap(GenericMap):
    """STEREO-SECCHI Heliospheric Imager (HI) Map.

    The HI is a wide-angle visible-light imaging system
    for the detection of coronal mass ejection (CME) events
    in interplanetary space and, in particular, of events
    directed towards the Earth.

    The Heliospheric imager consists of two instruments, the HI-1 and HI-2.
    The HI1 observes from 15-80 solar radii while HI2 observes from 80-215
    solar radii.

    References
    ----------
    * `STEREO Mission Page <http://stereo.gsfc.nasa.gov>`_
    * `STEREO SECCHI <http://secchi.nrl.navy.mil>`_
    * `HI Instrument Page <http://www.stereo.rl.ac.uk>`_
    """
    def __init__(self, data, header, **kwargs):

        GenericMap.__init__(self, data, header, **kwargs)
        self.meta['wavelnth'] = np.nan
        self.meta['waveunit'] = 'nm'
        self._nickname = "{0}-{1}".format(self.detector, self.observatory[-1])
        self.plot_settings['cmap'] = cm.get_cmap('stereohi{det!s}'.format(det=self.detector[-1]))
        self.plot_settings['norm'] = ImageNormalize(stretch=PowerStretch(0.25))

        # Try to identify when the FITS meta data does not have the correct
        # date FITS keyword
        if ('date_obs' in self.meta) and not('date-obs' in self.meta):
            self.meta['date-obs'] = self.meta['date_obs']

    @property
    def measurement(self):
        """
        Returns the type of data observed.
        """
        # TODO: This needs to do more than white-light.  Should give B, pB, etc.
        return "white-light"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an COR image"""
        return header.get('detector', '').startswith('HI')
