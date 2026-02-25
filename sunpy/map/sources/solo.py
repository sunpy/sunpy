"""
Solar Orbiter Map subclass definitions.
"""

import warnings
import numpy as np

import astropy.units as u
from astropy.coordinates import CartesianRepresentation
from astropy.visualization import AsinhStretch, ImageNormalize, AsymmetricPercentileInterval

from sunpy.coordinates import HeliocentricInertial
from sunpy.map import GenericMap
from sunpy.map.sources.source_type import source_stretch

__all__ = ["EUIMap", "METISMap"]


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
        self.plot_settings["norm"] = ImageNormalize(stretch=source_stretch(self.meta, AsinhStretch(0.01)), clip=False)

    @property
    def _rotation_matrix_from_crota(self):
        return super()._rotation_matrix_from_crota(crota_key="CROTA")

    @property
    def processing_level(self):
        if self.meta.get("level"):
            # The level number is prepended by the letter L
            return int(self.meta.get("level")[1:])

    @property
    def waveunit(self):
        # EUI JP2000 files do not have the WAVEUNIT key in the metadata.
        # However, the FITS files do.
        # The EUI metadata spec says the WAVELNTH key is always expressed
        # in Angstroms so we assume this if the WAVEUNIT is missing.
        return super().waveunit or u.Angstrom

    @property
    def _supported_observer_coordinates(self):
        return [
            (
                ("hcix_obs", "hciy_obs", "hciz_obs"),
                {
                    "x": self.meta.get("hcix_obs"),
                    "y": self.meta.get("hciy_obs"),
                    "z": self.meta.get("hciz_obs"),
                    "unit": u.m,
                    "representation_type": CartesianRepresentation,
                    "frame": HeliocentricInertial,
                },
            )
        ] + super()._supported_observer_coordinates

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EUI image"""
        is_solo = "solar orbiter" in str(header.get("obsrvtry", "")).lower()
        is_eui = str(header.get("instrume", "")).startswith("EUI")
        return is_solo and is_eui


class METISMap(GenericMap):
    """
    Metis Image Map.

    Metis is the multi-wavelength coronagraph on board the Solar Orbiter mission,
    dedicated to the study of the solar corona. It observes the outer atmosphere
    of the Sun simultaneously in:
    - Total and polarized visible light (580–640 nm), scattered by free electrons
    through the Thomson scattering
    - Ultraviolet band in the hydrogen Lyman-alpha line (121.6 nm), emitted by
    the few neutral hydrogen atoms surviving in the hot corona.

    By occulting the solar disk, Metis observes the faint coronal emission in an
    annular zone 1.6-2.9 deg wide, around the disk center. When Solar Orbiter
    is at its closest approach to the Sun, at the minimum perihelion of 0.28
    astronomical units, the annular zone is within 1.7 and 3.1 solar radii from
    disk center.

    Solar Orbiter was successfully launched on February 10th, 2020.

    References
    ----------
    * `Solar Orbiter Mission Page <https://sci.esa.int/web/solar-orbiter/>`_
    * `Metis Instrument Page <https://metis.oato.inaf.it/index.html>`_
    * Instrument Paper: :cite:t:`antonucci_metis_2020`
    """

    def __init__(self, data, header, **kwargs):
        """
        Initialize the METISMap class with the provided data and header.

        Validate that the header contains the required parameters.
        """
        if not any(k in header for k in ("RSUN_OBS", "SOLAR_R", "RADIUS")):
            header["RSUN_OBS"] = header["RSUN_ARC"]

        # Call the superclass (GenericMap) to initialize the map
        super().__init__(data, header, **kwargs)

        self._nickname = f"{self.instrument}/{self.meta['filter']}"
        self._prodtype = self.get_prodtype()
        self._contr_cut = self.get_contr_cut()
        self.update_plot_norm_settings()
        self.plot_settings["cmap"] = self._get_cmap_name()

    def get_prodtype(self):
        """
        Define the type of the Metis data product.

        Returns
        -------
        prodtype : `str`
            Name of the Metis data product.

        """
        btype_suff_dict = {
            "VL total brightness": ("-TB", "-TB"),
            "VL polarized brightness": ("-PB", "-PB"),
            "VL fixed-polarization intensity": ("-FP", "-Fix. Pol."),
            "VL polarization angle": ("-PA", "-Pol. Angle"),
            "Stokes I": ("-SI", "-Stokes I"),
            "Stokes Q": ("-SQ", "-Stokes Q"),
            "Stokes U": ("-SU", "-Stokes U"),
            "Pixel quality": ("-PQ", "-Pixel quality"),
            "Absolute error": ("-AE", "-Abs. err."),
            "Relative error": ("-RE", "-Rel. err."),
            "UV Lyman-alpha intensity": ("", ""),
        }

        btype = self.meta["btype"]
        prodtype = self.meta["filter"]

        if btype in btype_suff_dict:
            suff, nickname_add = btype_suff_dict[btype]
            prodtype += suff
            self._nickname += nickname_add
        else:
            raise ValueError(f"Error. self.meta['btype']='{btype}' is not known.")

        return prodtype

    @property
    def prodtype(self):
        """Product type identifier for this METIS data."""
        return self._prodtype

    @prodtype.setter
    def prodtype(self, value):
        raise AttributeError("Cannot manually set prodtype for METISMap")

    def get_contr_cut(self):
        """
        Define the contrast of the Metis data product.

        Returns
        -------
        contr_cut : `float` or `None`
            Contrast of the Metis data product, used for intensity scaling.

        """
        contr_cut_dict = {
            "VL-TB": 0.05,
            "VL-PB": 0.005,
            "VL-FP": 0.01,
            "VL-PA": 0.01,
            "VL-SQ": 0.01,
            "VL-SU": 0.01,
            "UV": 0.05,
            "VL-PQ": 0.0,
            "VL-AE": 0.1,
            "VL-RE": 0.02,
        }
        contr_cut_dict["VL-SI"] = contr_cut_dict["VL-TB"]
        contr_cut_dict["UV-PQ"] = contr_cut_dict["VL-PQ"]
        contr_cut_dict["UV-AE"] = contr_cut_dict["VL-AE"]
        contr_cut_dict["UV-RE"] = contr_cut_dict["VL-RE"]

        contr_cut = contr_cut_dict.get(self.prodtype, 0.0) if "L2" in self.meta["level"] else 0.0

        return contr_cut

    @property
    def contr_cut(self):
        """Contrast cutoff value for intensity scaling."""
        return self._contr_cut

    @contr_cut.setter
    def contr_cut(self, value):
        """
        Set the contrast cutoff value.

        Parameters
        ----------
        value : `float`
            New contrast cutoff value (0.0 to 1.0).

        Notes
        -----
        Changing this value will automatically update plot normalization settings.
        """
        if not isinstance(value, (int, float)):
            raise TypeError(f"contr_cut must be numeric, got {type(value)}")
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"contr_cut must be between 0.0 and 1.0, got {value}")
        self._contr_cut = value
        self.update_plot_norm_settings()

    def update_plot_norm_settings(self):
        """
        Update vmin and vmax values of plot_settings['norm'].

        Updates the plot normalization settings based on current data and
        contrast cutoff.
        """
        img_vlim = self.get_img_vlim()
        self.plot_settings["norm"] = ImageNormalize(vmin=img_vlim[0], vmax=img_vlim[1])

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determine whether the data is a L2 Metis product.

        Parameters
        ----------
        data : `numpy.ndarray`
            The image data array.
        header : `dict`
            The FITS header metadata.
        **kwargs : `dict`
            Additional keyword arguments.

        Returns
        -------
        `bool`
            ``True`` if data corresponds to a Metis product, otherwise ``False``.

        """
        instrume = header.get("INSTRUME", "").strip().upper()
        obsrvtry = header.get("OBSRVTRY", "").strip().upper()
        level = header.get("LEVEL", "").strip().upper()

        return ("METIS" in instrume) and ("SOLAR ORBITER" in obsrvtry) and ("L2" in level)

    def get_fov_rsun(self):
        """
        Return the Metis field of view in solar radii.

        Returns
        -------
        `tuple` : `(float, float, float)`
            Inner radius, outer radius, and detector board radius of the field,
            determined by the internal occulter, field stop and detector size,
            respectively. All values are in units of solar radii.

        """
        rmin_rsun = (self.meta["inn_fov"] * u.deg / self.rsun_obs).decompose()  # Inner FOV in rsun
        rmax_rsun = (self.meta["out_fov"] * u.deg / self.rsun_obs).decompose()  # Outer FOV in rsun
        board_deg = 2.9 * u.deg  # Detector board edge in deg
        board_rsun = (board_deg / self.rsun_obs).decompose()  # Board radius in rsun

        return rmin_rsun, rmax_rsun, board_rsun

    def mask_occs(self, mask_val=np.nan):
        """
        Mask the data in regions obscured by internal and external occulters.

        This method modifies the data array in-place, setting pixels outside
        the field of view to the specified mask value.

        Parameters
        ----------
        mask_val : `float`, optional
            The value to assign to masked pixels (outside the field of view).
            Default is ``np.nan``.

        Notes
        -----
        This method directly modifies ``self.data`` in-place. The plot
        normalization settings are automatically updated after masking.

        Warnings
        --------
        If CDELT1 and CDELT2 are not equal (non-square pixels), this method
        will warn and exit without performing masking, as the circular mask
        calculation assumes square pixels.

        """
        if self.scale[0] != self.scale[1]:
            warnings.warn(
                "Warning: CDELT1 != CDELT2 for {fname}. Exiting mask_occs method...".format(fname=self.meta["filename"])
            )
            return

        # Calculate occulter radii in pixels
        inn_fov = (self.meta["inn_fov"] * u.deg / self.scale[0]).decompose()  # in pix
        out_fov = (self.meta["out_fov"] * u.deg / self.scale[1]).decompose()  # in pix

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
        # Apply masks
        self.data[dist_inncen <= inn_fov] = mask_val
        self.data[dist_outcen >= out_fov] = mask_val
        # Update plot normalization after masking
        self.update_plot_norm_settings()

    def mask_bad_pix(self, qmat, mask_val=np.nan):
        """
        Mask bad-quality pixels in the Metis image.

        This method modifies the data array in-place, setting bad-quality
        pixels to the specified mask value based on the provided quality matrix.

        Parameters
        ----------
        qmat : `numpy.ndarray`
            Pixel quality matrix with the same shape as the image data.
            Expected values:
            - ``1`` : linear range (good-quality pixels)
            - ``0`` : close to 0 counts or close to saturation (bad-quality)
            - ``np.nan`` : exactly 0 count or saturated pixels (bad-quality)
        mask_val : `float`, optional
            The value to assign to masked bad pixels. Default is ``np.nan``.

        Raises
        ------
        ValueError
            If the quality matrix shape does not match the data shape.
        TypeError
            If qmat is not a numpy array.

        Notes
        -----
        This method directly modifies ``self.data`` in-place. The plot
        normalization settings are automatically updated after masking.

        """
        # Validate input type
        if not isinstance(qmat, np.ndarray):
            raise TypeError(f"qmat must be a numpy.ndarray, got {type(qmat).__name__}")
        # Validate shape compatibility
        if qmat.shape != self.data.shape:
            raise ValueError(
                f"Pixel quality matrix shape {qmat.shape} does not match "
                f"METISMap data shape {self.data.shape}. Cannot apply mask."
            )
        # Create mask: keep only pixels with value 1 (good quality)
        qmat_mask = qmat == 1
        self.data[~qmat_mask] = mask_val
        # Update plot normalization after masking
        self.update_plot_norm_settings()

    def _get_cmap_name(self):
        """
        Override the default implementation to handle Metis color maps.

        Returns
        -------
        cmap_string : `str`
            Name of the colormap for this Metis data product, formatted as
            'solo{instrument}{prodtype}' in lowercase.

        """
        cmap_string = f"solo{self.instrument}{self.prodtype}"
        cmap_string = cmap_string.lower()
        return cmap_string

    def get_img_vlim(self):
        """
        Return the intensity limits of the Metis image based on the contrast cutoff.

        Uses asymmetric percentile intervals to determine appropriate vmin and vmax
        values that exclude outliers while preserving the dynamic range of the data.

        Returns
        -------
        `tuple` : `(float, float)`
            The minimum and maximum intensity values for display, calculated
            using the contrast cutoff percentage.

        """
        vlim = AsymmetricPercentileInterval(self.contr_cut * 100, (1 - self.contr_cut) * 100).get_limits(self.data)

        return vlim
