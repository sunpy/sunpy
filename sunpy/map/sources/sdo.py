"""SDO Map subclass definitions"""

import numpy as np

import astropy.units as u
from astropy.coordinates import CartesianRepresentation, HeliocentricMeanEcliptic
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map.mapbase import GenericMap, SpatialPair
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
        self._nickname = self.detector
        self.plot_settings['cmap'] = self._get_cmap_name()
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False)

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

    @property
    def detector(self):
        return self.meta.get("detector", "AIA")

    @property
    def unit(self):
        unit_str = self.meta.get('bunit', self.meta.get('pixlunit'))
        if unit_str is None:
            return

        return self._parse_fits_unit(unit_str)

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
        if self.unit is not None and self.unit.is_equivalent(u.T):
            # Avoid JP2K images not having a norm due to UNIT8 data
            # This means they are not scaled correctly.
            if self.plot_settings.get('norm') is not None:
                # Magnetic field maps, not intensity maps
                self._set_symmetric_vmin_vmax()
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

    @property
    def detector(self):
        return self.meta.get("detector", "HMI")

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an HMI image"""
        return (str(header.get('INSTRUME', '')).startswith('HMI') and
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
        self.plot_settings['cmap'] = 'hmimag'
        self.plot_settings['norm'] = ImageNormalize(vmin=-1.5e3, vmax=1.5e3)

    @property
    def spatial_units(self):
        cunit1 = self.meta['cunit1']
        if cunit1 == 'Degree':
            cunit1 = 'deg'

        cunit2 = self.meta['cunit2']
        if cunit2 == 'Sine Latitude':
            cunit2 = 'deg'

        return SpatialPair(u.Unit(cunit1), u.Unit(cunit2))

    @property
    def scale(self):
        if self.meta['cunit2'] == 'Sine Latitude':
            # Since, this map uses the cylindrical equal-area (CEA) projection,
            # the spacing should be modified to 180/pi times the original value
            # Reference: Section 5.5, Thompson 2006
            return SpatialPair(np.abs(self.meta['cdelt1']) * self.spatial_units[0] / u.pixel,
                               180 / np.pi * self.meta['cdelt2'] * u.deg / u.pixel)

        return super().scale

    @property
    def date(self):
        """
        Image observation time.

        This is taken from the 'DATE-OBS' or 'T_OBS' keywords.
        """
        date = self._get_date('DATE-OBS')
        if date is None:
            return self._get_date('T_OBS')
        else:
            return date

    @property
    def unit(self):
        unit_str = self.meta.get('bunit', None)
        if unit_str is None:
            return
        # Maxwells aren't in the IAU unit style manual and therefore not a valid FITS unit
        # The mapbase unit property forces this validation, so we must override it to prevent it.
        return u.Unit(unit_str)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if header corresponds to an HMI synoptic map.
        """
        return str(header.get('TELESCOP', '')).endswith('HMI') and 'carrington synoptic chart' in str(header.get('CONTENT', '')).lower()
