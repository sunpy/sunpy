"""SOHO Map subclass definitions"""

import numpy as np
from matplotlib import colors

import astropy.units as u
from astropy.coordinates import CartesianRepresentation, HeliocentricMeanEcliptic
from astropy.visualization import PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy import log
from sunpy.map import GenericMap
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
        # Assume pixel units are arcesc if not given
        header['cunit1'] = header.get('cunit1', 'arcsec')
        header['cunit2'] = header.get('cunit2', 'arcsec')
        if 'waveunit' not in header or not header['waveunit']:
            header['waveunit'] = "Angstrom"
        super().__init__(data, header, **kwargs)

        self._nickname = self.detector
        self.plot_settings['cmap'] = self._get_cmap_name()
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.5)), clip=False)

    @property
    def detector(self):
        return "EIT"

    @property
    def rsun_obs(self):
        """
        Returns the solar radius as measured by EIT in arcseconds.
        """
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
    from 2 to 6 solar radii, overlaping the outer field-of-view of C1 from 2 to
    3 solar radii. The C3 telescope extends the field-of-view to 32 solar radii.

    SOHO was launched on 2 December 2 1995 into a sun-synchronous orbit.

    References
    ----------
    * `SOHO Mission Page <https://sohowww.nascom.nasa.gov/>`_
    """

    def __init__(self, data, header, **kwargs):

        header['cunit1'] = header['cunit1'].lower()
        header['cunit2'] = header['cunit2'].lower()

        super().__init__(data, header, **kwargs)

        # Fill in some missing or broken info
        # Test if change has already been applied
        if 'T' not in self.meta['date-obs']:
            datestr = "{date}T{time}".format(date=self.meta.get('date-obs',
                                                                self.meta.get('date_obs')
                                                                ),
                                             time=self.meta.get('time-obs',
                                                                self.meta.get('time_obs')
                                                                )
                                             )
            self.meta['date-obs'] = datestr

        # If non-standard Keyword is present, correct it too, for compatibility.
        if 'date_obs' in self.meta:
            self.meta['date_obs'] = self.meta['date-obs']
        self._nickname = self.instrument + "-" + self.detector
        self.plot_settings['cmap'] = 'soholasco{det!s}'.format(det=self.detector[1])
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, PowerStretch(0.5)), clip=False)

        # For Helioviewer images, clear rotation metadata, as these have already been rotated.
        # Also check that all CROTAn keywords exist to make sure that it's an untouched
        # Helioviewer file.
        if ('helioviewer' in self.meta
                and 'crota' in self.meta and 'crota1' in self.meta and 'crota2' in self.meta):
            log.debug("LASCOMap: Cleaning up CROTAn keywords "
                      "because the map has already been rotated by Helioviewer")
            self.meta.pop('crota')
            self.meta.pop('crota1')
            self.meta['crota2'] = 0

    @property
    def measurement(self):
        """
        Returns the type of data taken.
        """
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
        # Assume pixel units are arcesc if not given
        header['cunit1'] = header.get('cunit1', 'arcsec')
        header['cunit2'] = header.get('cunit2', 'arcsec')
        super().__init__(data, header, **kwargs)

        # Fill in some missing or broken info
        self._nickname = self.detector + " " + self.measurement
        vmin = np.nanmin(self.data)
        vmax = np.nanmax(self.data)
        threshold = max([abs(vmin), abs(vmax)])
        self.plot_settings['norm'] = colors.Normalize(-threshold, threshold)

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
    def detector(self):
        return "MDI"

    @property
    def waveunit(self):
        return "Angstrom"

    @property
    def measurement(self):
        """
        Returns the type of data in the map.
        """
        return "magnetogram" if self.meta.get('dpc_obsr', " ").find('Mag') != -1 else "continuum"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an MDI image"""
        return cls._is_mdi_map(header) and not cls._is_synoptic_map(header)


class MDISynopticMap(MDIMap):
    """
    SOHO MDI synoptic magnetogram Map.

    See the docstring of `MDIMap` for information on the MDI instrument.
    """
    def __init__(self, data, header, **kwargs):
        # FITS doesn't like "Degree" as a unit
        if header['cunit1'] == 'Degree':
            header['cunit1'] = 'deg'
        # Sine Latitude is not a valid unit - see Thompson 2006, section 5.5
        # for how to properly represent the cylindrical equal area (CEA) projection
        if header['cunit2'] == 'Sine Latitude':
            header['cdelt2'] = 180 / np.pi * header['cdelt2']
            header['cunit2'] = 'deg'
        if 'date-obs' not in header:
            header['date-obs'] = header['t_obs']
            header['date-obs'] = parse_time(header['date-obs']).isot
        for i in [1, 2]:
            if header[f'CRDER{i}'] == 'nan':
                header.pop(f'CRDER{i}')
        super().__init__(data, header, **kwargs)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an MDI image"""
        return cls._is_mdi_map(header) and cls._is_synoptic_map(header)
