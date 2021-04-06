"""SDO Map subclass definitions"""

import numpy as np

import astropy.units as u
from astropy.coordinates import CartesianRepresentation, HeliocentricMeanEcliptic
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy import log
from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ['AIAMap', 'HMIMap', 'HMISynopticMap']


class AIAMap(GenericMap):
    """AIA Image Map.

    The Atmospheric Imaging Assembly is a set of four telescopes that employ
    normal-incidence, multi-layer coated optics to provide narrow-band imaging
    of the Sun. It provides high resolution full-disk images of the corona and
    transition region up to 0.5 solar radii above the solar limb with 1.5
    arcsecond angular resolution and 12-second temporal resolution. It observes
    the Sun in the following seven extreme ultraviolet bandpasses: 94 A
    (Fe XVIII), 131 A (Fe VIII, XXI), 171 A (Fe IX), 193 A (Fe XII, XXIV),
    211 A (Fe XIV), 304 A (He II), 335 A (Fe XVI). One telescope observes
    in the visible 1600 A (C IV) and the nearby continuum (1700 A).

    Notes
    -----
    Observer location: The standard AIA FITS header provides the spacecraft location in multiple
    coordinate systems, including Heliocentric Aries Ecliptic (HAE) and Heliographic Stonyhurst
    (HGS).  SunPy uses the provided HAE coordinates due to accuracy concerns with the provided
    HGS coordinates, but other software packages may make different choices.

    References
    ----------
    * `SDO Mission Page <https://sdo.gsfc.nasa.gov/>`_
    * `Instrument Page <https://aia.lmsal.com>`_
    * `Fits Header keywords <http://jsoc.stanford.edu/doc/keywords/AIA/AIA02840_A_AIA-SDO_FITS_Keyword_Documents.pdf>`_
    * `Analysis Guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/>`_
    * `Instrument Paper <https://doi.org/10.1007/s11207-011-9776-8>`_
    * `wavelengths and temperature response reference <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/figures/aia_tel_resp.png>`_
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        # Fill in some missing info
        self.meta['detector'] = self.meta.get('detector', "AIA")
        if 'bunit' not in self.meta and 'pixlunit' in self.meta:
            # PIXLUNIT is not a FITS standard keyword
            self.meta['bunit'] = self.meta['pixlunit']
        self._nickname = self.detector
        self.plot_settings['cmap'] = self._get_cmap_name()
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False)
        # DN is not a FITS standard unit, so convert to counts
        if self.meta.get('bunit', None) == 'DN':
            self.meta['bunit'] = 'ct'
        if self.meta.get('bunit', None) == 'DN/s':
            self.meta['bunit'] = 'ct/s'

    @property
    def _supported_observer_coordinates(self):
        return [(('haex_obs', 'haey_obs', 'haez_obs'), {'x': self.meta.get('haex_obs'),
                                                        'y': self.meta.get('haey_obs'),
                                                        'z': self.meta.get('haez_obs'),
                                                        'unit': u.m,
                                                        'representation_type': CartesianRepresentation,
                                                        'frame': HeliocentricMeanEcliptic})
                ] + super()._supported_observer_coordinates

    @property
    def observatory(self):
        """
        Returns the observatory.
        """
        return self.meta.get('telescop', '').split('/')[0]

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an AIA image"""
        return str(header.get('instrume', '')).startswith('AIA')


class HMIMap(GenericMap):
    """HMI Image Map.

    HMI consists of a refracting telescope, a polarization selector,
    an image stabilization system, a narrow band tunable filter
    and two 4096 pixel CCD cameras. It observes the full solar disk in the Fe I
    absorption line at 6173 Angstrom with a resolution of 1 arc-second.
    HMI takes images in a sequence of tuning and polarizations at a 4-second
    cadence for each camera. One camera is dedicated to a 45 s Doppler and
    line-of-sight field sequence while the other to a 90 s vector field
    sequence.

    References
    ----------
    * `SDO Mission Page <https://sdo.gsfc.nasa.gov/>`_
    * `Instrument Page <http://hmi.stanford.edu>`_
    * `Analysis Guide <http://hmi.stanford.edu/doc/magnetic/guide.pdf>`_
    """

    def __init__(self, data, header, **kwargs):

        super().__init__(data, header, **kwargs)

        self.meta['detector'] = self.meta.get('detector', "HMI")
        self._nickname = self.detector

    @property
    def measurement(self):
        """
        Returns the measurement type.
        """
        return self.meta.get('content', '').split(" ")[0].lower()

    @property
    def observatory(self):
        """
        Returns the observatory.
        """
        return self.meta.get('telescop', '').split('/')[0]

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an HMI image"""
        return (str(header.get('TELESCOP', '')).endswith('HMI') and
                not HMISynopticMap.is_datasource_for(data, header))


class HMISynopticMap(HMIMap):
    """
    SDO/HMI Synoptic Map.

    Synoptic maps are constructed from HMI 720s line-of-sight magnetograms
    collected over a 27-day solar rotation.

    See `~sunpy.map.sources.sdo.HMIMap` for information on the HMI instrument.

    References
    ----------
    * `SDO Mission Page <https://sdo.gsfc.nasa.gov/>`__
    * `JSOC's HMI Synoptic Charts <http://jsoc.stanford.edu/HMI/LOS_Synoptic_charts.html>`__
    """
    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        if self.meta['cunit1'] == 'Degree':
            self.meta['cunit1'] = 'deg'

        if self.meta['cunit2'] == 'Sine Latitude':
            log.debug("Editing CUNIT2, CDELT1, CDLET2 keywords to the correct "
                      "values for a CEA projection.")
            self.meta['cunit2'] = 'deg'

            # Since, this map uses the cylindrical equal-area (CEA) projection,
            # the spacing should be modified to 180/pi times the original value
            # Reference: Section 5.5, Thompson 2006
            self.meta['cdelt2'] = 180 / np.pi * self.meta['cdelt2']
            self.meta['cdelt1'] = np.abs(self.meta['cdelt1'])

        if 'date-obs' not in self.meta and 't_obs' in self.meta:
            log.debug('Setting "DATE-OBS" keyword from "T_OBS"')
            self.meta['date-obs'] = self.meta['t_obs']

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if header corresponds to an HMI synoptic map.
        """
        return (str(header.get('TELESCOP', '')).endswith('HMI') and
                str(header.get('CONTENT', '')) ==
                'Carrington Synoptic Chart Of Br Field')
