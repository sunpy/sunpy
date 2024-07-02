"""SOHO Map subclass definitions"""

import numpy as np

import astropy.units as u
from astropy.coordinates import CartesianRepresentation, HeliocentricMeanEcliptic
from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy import log
from sunpy.map.mapbase import GenericMap, SpatialPair
from sunpy.map.sources.source_type import source_stretch
from sunpy.time import parse_time

__all__ = ['EITMap', 'LASCOMap', 'MDIMap', 'MDISynopticMap']


class EITMap(GenericMap):
    """
    SOHO EIT Image Map.

    SOHO EIT is an extreme ultraviolet (EUV) imager able to image the solar
    transition region and inner corona in four selected bandpasses,
    171 (Fe IX/X), 195 (Fe XII), 284 (Fe XV), and 304 (He II) Angstrom.

    SOHO was launched on 2 December 2 1995 into a sun-synchronous orbit and
    primary mission operations for SOHO EIT ended at the end of July 2010.

    References
    ----------
    * `SOHO Mission Page <https://sohowww.nascom.nasa.gov/>`_
    * `SOHO EIT Instrument Page <https://umbra.nascom.nasa.gov/eit/>`_
    * `SOHO EIT User Guide <https://umbra.nascom.nasa.gov/eit/eit_guide/>`_

    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self._nickname = self.detector
        self.plot_settings['cmap'] = self._get_cmap_name()
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.5)), clip=False)

    @property
    def date(self):
        # Old EIT data has date-obs in format of dd-JAN-yy so we use date_obs where available
        return self._get_date('date_obs') or super().date

    @property
    def spatial_units(self):
        """
        If not present in CUNIT{1,2} keywords, defaults to arcsec.
        """
        return SpatialPair(u.Unit(self.meta.get('cunit1', 'arcsec')),
                           u.Unit(self.meta.get('cunit2', 'arcsec')))

    @property
    def waveunit(self):
        """
        If WAVEUNIT FITS keyword isn't present, defaults to Angstrom.
        """
        unit = self.meta.get("waveunit", "Angstrom") or "Angstrom"
        return u.Unit(unit)

    @property
    def detector(self):
        return "EIT"

    @property
    def rsun_obs(self):
        return u.Quantity(self.meta['solar_r'] * self.meta['cdelt1'], 'arcsec')

    @property
    def _supported_observer_coordinates(self):
        return [(('hec_x', 'hec_y', 'hec_z'), {'x': self.meta.get('hec_x'),
                                               'y': self.meta.get('hec_y'),
                                               'z': self.meta.get('hec_z'),
                                               'unit': u.km,
                                               'representation_type': CartesianRepresentation,
                                               'frame': HeliocentricMeanEcliptic})
                ] + super()._supported_observer_coordinates

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EIT image"""
        return header.get('instrume') == 'EIT'


class LASCOMap(GenericMap):
    """
    SOHO LASCO Image Map

    The Large Angle and Spectrometric COronagraph (LASCO) is a set of three
    Lyot-type coronagraphs (C1, C2, and C3) that image the solar corona from
    1.1 to 32 solar radii.

    The C1 images rom 1.1 to 3 solar radii. The C2 telescope images the corona
    from 2 to 6 solar radii, overlapping the outer field-of-view of C1 from 2 to
    3 solar radii. The C3 telescope extends the field-of-view to 32 solar radii.

    SOHO was launched on 2 December 2 1995 into a sun-synchronous orbit.

    References
    ----------
    * `SOHO Mission Page <https://sohowww.nascom.nasa.gov/>`_
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self.plot_settings['cmap'] = f'soholasco{self.detector[1]!s}'
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.5)), clip=False)

    @property
    def spatial_units(self):
        return SpatialPair(u.Unit(self.meta.get('cunit1').lower()),
                           u.Unit(self.meta.get('cunit2').lower()))

    @property
    def rotation_matrix(self):
        # For Helioviewer images, clear rotation metadata, as these have already been rotated.
        # Also check that all CROTAn keywords exist to make sure that it's an untouched
        # Helioviewer file.
        if ('helioviewer' in self.meta and
                'crota' in self.meta and
                'crota1' in self.meta and
                'crota2' in self.meta):
            log.debug("LASCOMap: Ignoring CROTAn keywords "
                      "because the map has already been rotated by Helioviewer")
            return np.identity(2)
        else:
            return super().rotation_matrix

    @property
    def date(self):
        date = self.meta.get('date-obs', self.meta.get('date_obs'))
        # In case someone fixes the header
        if 'T' in date:
            return parse_time(date)

        time = self.meta.get('time-obs', self.meta.get('time_obs'))
        return parse_time(f"{date}T{time}")

    @property
    def nickname(self):
        filter = self.meta.get('filter', '')
        return f'{self.instrument}-{self.detector} {filter}'

    @nickname.setter
    def nickname(self, value):
        raise AttributeError("Cannot manually set nickname for LASCOMap")

    @property
    def measurement(self):
        # TODO: This needs to do more than white-light.  Should give B, pB, etc.
        return "white-light"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an LASCO image."""
        return header.get('instrume') == 'LASCO'


class MDIMap(GenericMap):
    """
    SOHO MDI Image Map

    The Michelson Doppler Imager (MDI) is a white light refracting telescope
    which feeds sunlight through a series of filters onto a CCD camera. Two
    tunable Michelson interformeters define a 94 mAngstrom bandpass that can be
    tuned across the Ni 6768 Angstrom solar absorption line.

    MDI measures line-of-sight motion (Dopplergrams), magnetic field
    (magnetograms), and brightness images of the full solar disk at several
    resolutions (4 arc-second to very low resolution) and a fixed selected
    region in higher resolution (1.2 arc-second).

    SOHO was launched on 2 December 2 1995 into a sun-synchronous orbit and
    SOHO MDI ceased normal science observations on 12 April 2011.

    References
    ----------
    * `SOHO Mission Page <https://sohowww.nascom.nasa.gov/>`_
    * `SOHO MDI Instrument Page <http://soi.stanford.edu>`_
    * `SOHO MDI Fits Header keywords <http://soi.stanford.edu/sssc/doc/keywords.html>`_
    * `SOHO MDI Instrument Paper <https://doi.org/10.1007/978-94-009-0191-9_5>`_
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        if self.unit is not None and self.unit.is_equivalent(u.T):
            # Magnetic field maps, not intensity maps
            self._set_symmetric_vmin_vmax()

    @property
    def _date_obs(self):
        if 'T' in self.meta.get('date-obs', ''):
            # Helioviewer MDI files have the full date in DATE_OBS, but we still
            # want to let normal FITS files use DATE-OBS
            return parse_time(self.meta['date-obs'])
        elif 'date_obs' in self.meta:
            return parse_time(self.meta['date_obs'])

    @property
    def unit(self):
        bunit = self.meta.get('bunit', None)
        if bunit is not None and bunit.lower() == 'arbitrary intensity units':
            return u.dimensionless_unscaled
        return super().unit

    @property
    def spatial_units(self):
        """
        If not present in CUNIT{1,2} keywords, defaults to arcsec.
        """
        return SpatialPair(u.Unit(self.meta.get('cunit1', 'arcsec')),
                           u.Unit(self.meta.get('cunit2', 'arcsec')))

    @staticmethod
    def _is_mdi_map(header):
        return header.get('instrume') == 'MDI' or header.get('camera') == 'MDI'

    @staticmethod
    def _is_synoptic_map(header):
        return 'Synoptic Chart' in header.get('CONTENT', '')

    @property
    def _supported_observer_coordinates(self):
        return [(('obs_l0', 'obs_b0', 'obs_dist'), {'lon': self.meta.get('obs_l0'),
                                                    'lat': self.meta.get('obs_b0'),
                                                    'radius': self.meta.get('obs_dist'),
                                                    'unit': (u.deg, u.deg, u.AU),
                                                    'frame': "heliographic_carrington"}),
                ] + super()._supported_observer_coordinates

    @property
    def instrument(self):
        return "MDI"

    @property
    def waveunit(self):
        """
        Always assumed to be Angstrom.
        """
        return "Angstrom"

    @property
    def measurement(self):
        """
        Returns the measurement type.
        """
        return self.meta.get('CONTENT', '')

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an MDI image"""
        return cls._is_mdi_map(header) and not cls._is_synoptic_map(header)


class MDISynopticMap(MDIMap):
    """
    SOHO MDI synoptic magnetogram Map.

    See the docstring of `MDIMap` for information on the MDI instrument.
    """
    @property
    def date(self):
        """
        Image observation time.

        This is taken from the 'DATE-OBS' or 'T_OBS' keywords.
        """
        time = self._get_date('date-obs')
        if time is None:
            return self._get_date('t_obs')

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
    def unit(self):
        bunit = self.meta.get('bunit', None)
        if bunit is None:
            return
        # Maxwells aren't in the IAU unit style manual and therefore not a valid FITS unit
        # The mapbase unit property forces this validation, so we must override it to prevent it.
        return u.Unit(bunit)

    @property
    def scale(self):
        if self.meta['cunit2'] == 'Sine Latitude':
            # Since, this map uses the cylindrical equal-area (CEA) projection,
            # the spacing should be modified to 180/pi times the original value
            # Reference: Section 5.5, Thompson 2006
            return SpatialPair(np.abs(self.meta['cdelt1']) * self.spatial_units[0] / u.pixel,
                               180 / np.pi * self.meta['cdelt2'] * u.deg / u.pixel)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an MDI image"""
        return cls._is_mdi_map(header) and cls._is_synoptic_map(header)
