"""
Solar Orbiter Map subclass definitions.
"""
import numpy as np
from matplotlib import cm
from matplotlib.colors import CenteredNorm

import astropy.units as u
from astropy.coordinates import CartesianRepresentation
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval

from sunpy import log
from sunpy.coordinates import HeliocentricInertial
from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch
from sunpy.util.exceptions import warn_user

__all__ = ['EUIMap', 'PHIMap', 'METISMap']


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
    * `PHI Instrument Page <https://www.mps.mpg.de/solar-physics/solar-orbiter-phi>`__
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
            self.plot_settings['cmap'] = 'grey'
            self.plot_settings['norm'] = CenteredNorm()
        elif btype == 'magnetic field strength' or btype == 'bmag':
            self.plot_settings['cmap'] = 'plasma'
        elif btype == 'magnetic field inclination' or btype == 'binc':
            self.plot_settings['cmap'] = 'RdGy'
            self.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=180, clip=True)
        elif btype == 'magnetic field azimuth' or btype == 'bazi':
            self.plot_settings['cmap'] = 'twilight'
            self.plot_settings['norm'] = ImageNormalize(vmin=0, vmax=180, clip=True)
        elif btype == 'los velocity' or btype == 'vlos':
            self.plot_settings['cmap'] = 'RdBu_r'
            self.plot_settings['norm'] = CenteredNorm()
        elif btype == 'intensity' or btype == 'icnt':
            self.plot_settings['cmap'] = 'inferno'

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

class METISMap(GenericMap):
    r"""
    Metis Map.

    Metis is the multi-wavelength coronagraph on board the Solar Orbiter mission,
    dedicated to the study of the solar corona. It observes the outer atmosphere
    of the Sun simultaneously in:

    - Total and polarized visible light (580-640 nm), scattered by free electrons
      through the Thomson scattering
    - Ultraviolet band in the hydrogen Lyman-alpha line (121.6 nm), emitted by
      the few neutral hydrogen atoms surviving in the hot corona.

    By occulting the solar disk, Metis observes the faint coronal emission in an
    annular zone 1.6-2.9 deg wide, around the disk center. When Solar Orbiter
    is at its closest approach to the Sun, at the minimum perihelion of 0.28
    astronomical units, the annular zone is within 1.7 and 3.1 solar radii from
    disk center.

    Solar Orbiter was successfully launched on February 10th, 2020.

    Notes
    -----

    This map source clips out the upper and lower :math:`0.5\%` of the pixels as many Metis files have some extreme spikes.
    These clipped pixels are shown in red by default.

    References
    ----------
    * `Solar Orbiter Mission Page <https://sci.esa.int/web/solar-orbiter/>`_
    * `Metis Instrument Page <https://metis.oato.inaf.it/index.html>`_
    * Instrument Paper: :cite:t:`antonucci_metis_2020`
    """
    _BTYPE_SUFF_DICT = {
        "VL total brightness":              ("-TB",  "-TB"),
        "VL polarized brightness":          ("-PB",  "-PB"),
        "VL fixed-polarization intensity":  ("-FP",  "-Fix. Pol."),
        "VL polarization angle":            ("-PA",  "-Pol. Angle"),
        "Stokes I":                         ("-SI",  "-Stokes I"),
        "Stokes Q":                         ("-SQ",  "-Stokes Q"),
        "Stokes U":                         ("-SU",  "-Stokes U"),
        "Pixel quality":                    ("-PQ",  "-Pixel quality"),
        "Absolute error":                   ("-AE",  "-Abs. err."),
        "Relative error":                   ("-RE",  "-Rel. err."),
        "UV Lyman-alpha intensity":         ("",     ""),
    }

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        btype = self.meta["btype"]
        nickname_add = ""
        if btype in METISMap._BTYPE_SUFF_DICT:
            _, nickname_add = METISMap._BTYPE_SUFF_DICT[btype]
        self._nickname = f"{self.instrument}/{self.meta['filter']}{nickname_add}"
        cmap = cm.get_cmap(f"solo{self.instrument}{self.measurement}".lower())

        self.plot_settings["norm"] = ImageNormalize(
            self.data,
            interval=PercentileInterval(99.5),
            clip=False,
        )
        # Set the under and over value for the colormap to be red, this will highlight pixels where the value is clipped
        cmap = cmap.with_extremes(under="red", over="red")
        self.plot_settings["cmap"] = cmap

    @property
    def rsun_obs(self):
        """
        Angular radius of the observation from Sun center.
        This value is taken (in order of preference) from the 'RSUN_OBS',
        'SOLAR_R', 'RADIUS', or 'RSUN_ARC' FITS keywords. If none of these keys are present,
        the angular radius is calculated from
        `~sunpy.map.GenericMap.rsun_meters` and `~sunpy.map.GenericMap.dsun`.
        """
        rsun_arcseconds = self._rsun_obs_no_default
        if rsun_arcseconds is not None:
            return rsun_arcseconds * u.arcsec

        if "RSUN_ARC" in self.meta:
            return self.meta["RSUN_ARC"] * u.arcsec

        return super().rsun_obs

    @property
    def measurement(self):
        """
        Derive the product type string from the BTYPE header keyword.

        Returns
        -------
        str
            Product type string.
        """

        if "filter" not in self.meta:
            return super().measurement
        prodtype = self.meta["filter"]
        btype = self.meta.get("btype", "")

        suff, _ = METISMap._BTYPE_SUFF_DICT.get(btype, ("", None))
        return prodtype + suff


    @property
    def mask(self):
        """
        Mask pixels obscured by the internal and external occulters.

        Returns a mask set to ``True`` for pixels inside the inner occulter
        and outside the outer field of view.

        Warnings
        --------
        If CDELT1 and CDELT2 are not equal (non-square pixels), this method
        will warn and return ``None`` without computing a mask, as the
        circular mask calculation assumes square pixels.
        """
        if missing_keys := {"inn_fov", "out_fov",
                            "io_xcen", "io_ycen",
                            "fs_xcen", "fs_ycen"}.difference(self.meta.keys()):
            warn_user(f"Missing {', '.join(missing_keys)} keys required to calculate occulter mask")
            return super().mask

        if self._mask is None:
            if self.scale[0] != self.scale[1]:
                log.info(
                    "CDELT1 != CDELT2, can not automatically compute the occulter mask."
                )
                return None

            # Calculate occulter radii in pixels
            inn_fov = (self.meta["inn_fov"] * u.deg / self.scale[0])
            out_fov = (self.meta["out_fov"] * u.deg / self.scale[1])

            # Create coordinate grids
            x = np.arange(0, self.meta["naxis1"], 1)
            y = np.arange(0, self.meta["naxis2"], 1)
            xx, yy = np.meshgrid(x, y, sparse=True)

            # Calculate distance from internal occulter center
            in_xcen = (self.meta["io_xcen"] - 1) * u.pix
            in_ycen = (self.meta["io_ycen"] - 1) * u.pix
            dist_inncen = np.sqrt((xx * u.pix - in_xcen) ** 2 + (yy * u.pix - in_ycen) ** 2)

            # Calculate distance from external occulter/field stop center
            # NOTE: Workaround for DR1 data where fs_*cen keywords are not correctly defined.
            # In DR1, fs_xcen and fs_ycen are not available and as consequence the corresponding keywords in a FITS file are assigned to crpix1 and crpix2, respectively.
            # When this occurs, sun_xcen and sun_ycen should be used as approximate coordinates of the field stop center.
            # This problem is solved in DR2.
            if self.meta["fs_xcen"] == self.meta["crpix1"] and self.meta["fs_ycen"] == self.meta["crpix2"]:
                # DR1 workaround: use sun center instead; For the DR1 data fs_*cen keywords are not defined correctly in a FITS file header
                out_xcen = (self.meta["sun_xcen"] - 1) * u.pix
                out_ycen = (self.meta["sun_ycen"] - 1) * u.pix
            else:
                # Normal case: use field stop center
                out_xcen = (self.meta["fs_xcen"] - 1) * u.pix
                out_ycen = (self.meta["fs_ycen"] - 1) * u.pix

            dist_outcen = np.sqrt((xx * u.pix - out_xcen) ** 2 + (yy * u.pix - out_ycen) ** 2)

            mask_inner = dist_inncen <= inn_fov
            mask_outer = dist_outcen >= out_fov
            self._mask = mask_inner | mask_outer

        return self._mask


    @mask.setter
    def mask(self, value):
        self._mask = value


    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Return ``True`` if ``header`` corresponds to a METIS L2 product.

        Parameters
        ----------
        data : numpy.ndarray
        header : dict-like
        """
        instrume = str(header.get("INSTRUME", "")).strip().upper()
        obsrvtry = str(header.get("OBSRVTRY", "")).strip().upper()
        level = str(header.get("LEVEL", "")).strip().upper()

        return ("METIS" in instrume and "SOLAR ORBITER" in obsrvtry and "L2" in level)
