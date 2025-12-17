"""
Solar Orbiter Map subclass definitions.
"""
import astropy.units as u
from astropy.coordinates import CartesianRepresentation
from astropy.visualization import AsinhStretch, ImageNormalize
import warnings

from sunpy.coordinates import HeliocentricInertial
from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch
from sunpy.util.exceptions import SunpyMetadataWarning, warn_user

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
    * HRT On-Ground Pipeline Paper: :cite:t:`sinjan_ground_2022`
    * HRT Magnetic Field Comparison with HMI Paper: :cite:t:`sinjan_mag_hrt_hmi_comparison_2023`
    * HRT Velocity Comparison with HMI Paper: :cite:t:`calchetti_vlos_hrt_hmi_comparison_2025`
    * FDT Magnetic Field Comparison with HMI Paper: :cite:t:`moreno_vacas_mag_fdt_hmi_comparison_2024`
    """

    """
    TODO:
    - stokes files as separate map class, as it has 4D array, need NDCube
    - get all (normally 6) observed wavelengths (they are corrected for the orbital
      velocity of the spacecraft)
    - test compatibility with both HRT and FDT data
    - raise SunpyMetadataWarning if CAL_WCS is not True (and is HRT)
    """

    def __init__(self, data, header, **kwargs):
        header = header.copy()
        if header.get('BUNIT') == "Normalised Intensity":
            header['BUNIT'] = "" # dimensionless
        elif header.get('BUNIT') == 'Degrees':
            header['BUNIT'] = "deg"

        super().__init__(data, header, **kwargs)
        self._nickname = self.detector

        if self.meta.get('btype', '').lower() == 'blos':
            self.plot_settings['cmap'] = 'hmimag'
            self.plot_settings['norm'] = ImageNormalize(vmin=-1.5e3, vmax=1.5e3)
        elif self.meta.get('btype', '').lower() == 'bmag':
            self.plot_settings['cmap'] = 'rainbow'
            self.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=2.5e3)
        elif self.meta.get('btype', '').lower() == 'binc':
            self.plot_settings['cmap'] = 'RdGy'
            self.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=180)
        elif self.meta.get('btype', '').lower() == 'bazi':
            self.plot_settings['cmap'] = 'hsv'
            self.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=180)
        elif self.meta.get('btype', '').lower() == 'vlos':
            self.plot_settings['cmap'] = 'RdBu_r'
            self.plot_settings['norm'] = ImageNormalize(vmin=-2, vmax=2)
        elif self.meta.get('btype', '').lower() == 'icnt':
            self.plot_settings['cmap'] = 'gist_heat'
            self.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=1.2)

        # Warn user if WCS is not calibrated for HRT data
        def str_to_bool(s):
            return str(s).lower() in ['true', '1', 't']

        if self.detector == 'HRT':
            try:
                cal_wcs = self.meta.get('cal_wcs', None) 
                if type(cal_wcs) is str: #older data versions have CAL_WCS as a string
                    cal_wcs = str_to_bool(cal_wcs)
                if not cal_wcs:
                    warn_user("The WCS of this SO/PHI-HRT PHIMap may not be fully calibrated. "
                              "Use caution when using the WCS for scientific analysis.")
            except Exception:
                warn_user("Could not determine if the WCS of this SO/PHI-HRT PHIMap is calibrated. "
                    "Use caution when using the WCS for scientific analysis.")

            
    @property
    def _rotation_matrix_from_crota(self):
        return super()._rotation_matrix_from_crota(crota_key='CROTA')
    
    @property
    def processing_level(self):
        if self.meta.get('level'):
            # Low Latency data products have levels LL02, LL03
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
    def measurement(self):
        """
        Returns the measurement type.
        """
        return self.meta.get('btype', 'Unknown')
    
    @property
    def observatory(self):
        """
        Returns the observatory.
        """
        return self.meta.get('obsrvtry', 'Solar Orbiter')

    
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to a PHI image"""
        is_solo = 'solar orbiter' in str(header.get('obsrvtry', '')).lower()
        is_phi = str(header.get('instrume', '')).startswith('PHI')
        is_not_phi_stokes = str(header.get('btype', '')).lower() != 'stokes'
        return is_solo and is_phi and is_not_phi_stokes 
        #future higher level data products will need to be checked here