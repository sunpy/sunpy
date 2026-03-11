"""
Solar Orbiter Map subclass definitions.
"""
import numpy as np
from matplotlib.colors import TwoSlopeNorm

import astropy.units as u
from astropy.coordinates import CartesianRepresentation
from astropy.visualization import AsinhStretch, ImageNormalize

from sunpy.coordinates import HeliocentricInertial
from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch
from sunpy.util.exceptions import warn_user

__all__ = ['EUIMap', 'PHIMap']


class EUIMap(GenericMap):
    """
    EUI Image Map

    The Extreme Ultraviolet Imager (EUI) is a remote sensing instrument onboard the
    Solar Orbiter (SolO) spacecraft. EUI has three telescopes that image the Sun in
    Lyman-alpha (1216 Å) and the EUV (174 Å and 304 Å). The three telescopes are the
    Full Sun Imager (FSI) and two High Resolution Imagers (HRI). The FSI images the
    whole Sun in both 174 Å and 304 Å. The EUV and Lyman-alpha HRI telescopes image a
    1000"-by-1000" patch in 174 Å and 1216 Å, respectively.

    References
    ----------
    * `Solar Orbiter Mission Page <https://sci.esa.int/web/solar-orbiter/>`__
    * `EUI Instrument Page <https://www.sidc.be/EUI/about/instrument>`__
    * Instrument Paper: :cite:t:`rochus_solar_2020`
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname = self.detector
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False)

    @property
    def _rotation_matrix_from_crota(self):
        return super()._rotation_matrix_from_crota(crota_key='CROTA')

    @property
    def processing_level(self):
        if self.meta.get('level'):
            # The level number is prepended by the letter L
            return int(self.meta.get('level')[1:])

    @property
    def waveunit(self):
        # EUI JP2000 files do not have the WAVEUNIT key in the metadata.
        # However, the FITS files do.
        # The EUI metadata spec says the WAVELNTH key is always expressed
        # in Angstroms so we assume this if the WAVEUNIT is missing.
        return super().waveunit or u.Angstrom

    @property
    def _supported_observer_coordinates(self):
        return [(('hcix_obs', 'hciy_obs', 'hciz_obs'),
                 {'x': self.meta.get('hcix_obs'),
                  'y': self.meta.get('hciy_obs'),
                  'z': self.meta.get('hciz_obs'),
                  'unit': u.m,
                  'representation_type': CartesianRepresentation,
                  'frame': HeliocentricInertial})] + super()._supported_observer_coordinates

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EUI image"""
        is_solo = 'solar orbiter' in str(header.get('obsrvtry', '')).lower()
        is_eui = str(header.get('instrume', '')).startswith('EUI')
        return is_solo and is_eui


class PHIMap(GenericMap):
    """
    PHI Image Map

    The Polarimetric and Helioseismic Imager (PHI) is a remote sensing instrument
    onboard the Solar Orbiter (SolO) spacecraft. PHI measures the Zeeman
    polarization and Doppler shifts of the Fe I 6173 Å spectral line to study the
    magnetic field and velocity field in the solar photosphere.

    It has two telescopes: the High Resolution Telescope (HRT) and the Full Disk
    Telescope (FDT). The HRT provides high-resolution observations of a small region
    of the solar surface, while the FDT always captures full-disk images of the Sun,
    no matter the distance of Solar Orbiter from the Sun. They cannot observe
    simultaneously.

    At perihelion (0.28 AU), the HRT achieves a spatial resolution (two pixels) of
    about 200 km on the solar surface.
    At perihelion (0.28 AU), the FDT achieves a spatial resolution (two pixels) of
    about 780 km on the solar surface.

    The Level 2 PHI STOKES data product (the raw Stokes polarimetric images) normally
    contain a 4D array of dimensions (y,x,4,6) (4 polarization states x 6 wavelengths).
    The STOKES currently cannot be handled by PHIMap.

    References
    ----------
    * `Solar Orbiter Mission Page <https://sci.esa.int/web/solar-orbiter/>`__
    * `PHI Instrument Page <https://www.mps.mpg.de/solar-orbiter/phi>`__
    * `PHI-HRT Data Quick Look Page <https://www.uv.es/jublanro/phidata_hrt.html>`__
    * `PHI-FDT Data Quick Look Page <https://www.uv.es/jublanro/phidata_fdt.html>`__
    * Instrument Paper: :cite:t:`solanki_polarimetric_2020`
    * HRT Instrument Paper: :cite:t:`gandorfer_high_resolution_2018`
    * HRT On-Ground Pipeline Paper: :cite:t:`sinjan_on-ground_2022`
    * HRT Magnetic Field Comparison with HMI Paper: :cite:t:`sinjan_mag_hrt_hmi_comparison_2023`
    * HRT Velocity Comparison with HMI Paper: :cite:t:`calchetti_vlos_hrt_hmi_comparison_2025`
    * FDT Magnetic Field Comparison with HMI Paper: :cite:t:`moreno_vacas_mag_fdt_hmi_comparison_2024`
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname = self.detector

        btype = self.meta.get('btype','').strip().lower()

        if btype == 'los magnetic field strength' or btype == 'blos':  # older versions may have blos
            self.plot_settings['cmap'] = 'hmimag'
            self.plot_settings['norm'] = ImageNormalize(vmin=-1.5e3, vmax=1.5e3, clip=False)
        elif btype == 'magnetic field strength' or btype == 'bmag':
            self.plot_settings['cmap'] = 'rainbow'
        elif btype == 'magnetic field inclination' or btype == 'binc':
            self.plot_settings['cmap'] = 'RdGy'
            self.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=180, clip=True)
        elif btype == 'magnetic field azimuth' or btype == 'bazi':
            self.plot_settings['cmap'] = 'hsv'
            self.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=180, clip=True)
        elif btype == 'los velocity' or btype == 'vlos':
            self.plot_settings['cmap'] = 'RdBu_r'
            v=np.nanmax(np.abs(self.data))
            self.plot_settings['norm'] = TwoSlopeNorm(vcenter=0, vmin=-v, vmax=v)
        elif btype == 'intensity' or btype == 'icnt':
            self.plot_settings['cmap'] = 'gist_heat'

        if self.detector == 'HRT':
            try:
                cal_wcs = self.meta.get('cal_wcs', None)
                if type(cal_wcs) is str:  # older data versions have CAL_WCS as a string
                    cal_wcs = str(cal_wcs).lower() in ['true', '1', 't']
                if not cal_wcs:
                    warn_user("The WCS of this SO/PHI-HRT PHIMap may not be fully calibrated. "
                              "Use caution when using the WCS for scientific analysis.")
            except Exception:
                warn_user("Could not determine if the WCS of this SO/PHI-HRT PHIMap is calibrated. "
                    "Use caution when using the WCS for scientific analysis.")

        # Check for on-board averaged dataset - there is no way to know this from the filename itself
        accumulations = int(self.meta.get('accrowit', None))
        if accumulations > 1:
            warn_user(f"This dataset was generated using on-board averaging of {accumulations} individual observations.")

    @property
    def processing_level(self):
        if self.meta.get('level'):
            # JP2 Low Latency data products have levels L3
            # Some low latency FITS files may have level LL01 (raw) or LL02 (reduced) - but not publicly available on SOAR
            # LL01 (raw) are rarely downlinked as LL is processed on board
            if self.meta.get('level').startswith('LL'):
                return int(self.meta.get('level')[3:])
            else:
                # For Regular data products, the level number is prepended by the letter L
                return int(self.meta.get('level')[1:])

    @property
    def waveunit(self):
        """
        The `~astropy.units.Unit` of the wavelength of this observation.

        PHI JP2000 files (FDT Low Latency files) may not have the WAVEUNIT key in the metadata.
        However, the FITS files do.
        The PHI metadata spec says the WAVELNTH key is always expressed
        in Angstroms so we assume this if the WAVEUNIT is missing.
        """
        return super().waveunit or u.Angstrom

    @property
    def unit(self):
        """
        Returns the unit
        """
        unit_str = self.meta.get('bunit', None)
        if unit_str is None:
            return
        elif unit_str == "Normalised Intensity" or unit_str == "Normalized":
            return u.dimensionless_unscaled
        elif unit_str == 'Degrees':
            return u.deg

        return u.Unit(unit_str)

    @property
    def measurement(self):
        """
        Returns the measurement type.
        """
        return self.meta.get('btype', 'Unknown')

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to a PHI image"""
        is_solo = 'solar orbiter' in str(header.get('obsrvtry', '')).lower()
        is_phi = str(header.get('instrume', '')).startswith('PHI')
        is_not_phi_stokes = str(header.get('btype', '')).lower() != 'stokes'
        return is_solo and is_phi and is_not_phi_stokes
