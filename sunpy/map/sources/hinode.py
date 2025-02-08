"""Hinode XRT and SOT Map subclass definitions"""
import astropy.units as u
from astropy.visualization import ImageNormalize, LogStretch

from sunpy.map.mapbase import GenericMap, SpatialPair
from sunpy.map.sources.source_type import source_stretch

__author__ = ["Jack Ireland, Jose Ivan Campos-Rozo, David Perez-Suarez"]
__email__ = "jack.ireland@nasa.gov"
__all__ = ['XRTMap', 'SOTMap']


def _lower_list(alist):
    return [item.lower() for item in alist]


class XRTMap(GenericMap):
    """
    Hinode XRT map definition.

    The X-Ray Telescope (XRT) is a high resolution grazing incidence telescope,
    which is a successor to Yohkoh. It provides 2-arcsecond resolution images
    of the highest temperature solar coronal material,
    from 1,000,000 to 10,000,000 Kelvin.

    Hinode was launched on 22 September 2006 into a sun-synchronous orbit.

    Notes
    -----
    XRT files do not normally specify the heliographic longitude of the spacecraft,
    so sunpy silently assumes that the spacecraft is at zero Stonyhurst heliographic
    longitude (i.e., the same longitude as Earth) for L1 files. This assumption is
    safe for nearly all analyses due to Hinode's orbital altitude of only ~600 km.
    This assumption is not made if ``HGLN_OBS`` and ``HGLT_OBS`` values have been
    explicitly added to the metadata.

    References
    ----------
    * `Hinode Mission Page <https://solarb.msfc.nasa.gov/index.html>`__
    * `XRT Instrument Page <https://xrt.cfa.harvard.edu/>`__
    * `Fits header reference <https://hinode.nao.ac.jp/uploads/2016/04/22/SB_MW_Key13.pdf>`__
    * `Hinode User Guide <https://hinode.nao.ac.jp/en/for-researchers/analysis-guide/>`__
    * `XRT Analysis Guide <https://xrt.cfa.harvard.edu/science/tutorials.php>`__
    * `Coronal Temperature Diagnostic Capability of the Hinode/X-Ray Telescope Based on Self-Consistent Calibration <https://arxiv.org/abs/1011.2867>`__
    """
    filter_wheel1_measurements = ["Al_med", "Al_poly", "Be_med",
                                  "Be_thin", "C_poly", "Open"]
    filter_wheel2_measurements = ["Open", "Al_mesh", "Al_thick",
                                  "Be_thick", "Gband", "Ti_poly"]

    def __init__(self, data, header, **kwargs):
        fw1 = header.get('EC_FW1_', '')
        if fw1.lower() not in _lower_list(self.filter_wheel1_measurements):
            raise ValueError(f'Unexpected filter wheel 1 {fw1} in header.')
        fw2 = header.get('EC_FW2_', '')
        if fw2.lower() not in _lower_list(self.filter_wheel2_measurements):
            raise ValueError(f'Unexpected filter wheel 2 {fw2} in header.')
        super().__init__(data, header, **kwargs)
        self.plot_settings['cmap'] = 'hinodexrt'
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, LogStretch()), clip=False)

    @property
    def coordinate_system(self):
        """
        Override the default implementation to handle SOTMap-specific logic for CTYPE values.
        """
        ctype1, ctype2 = self.meta['ctype1'], self.meta['ctype2']
        if ctype1.lower() in ("solar-x", "solar_x"):
            ctype1 = 'HPLN-TAN'
        if ctype2.lower() in ("solar-y", "solar_y"):
            ctype2 = 'HPLT-TAN'
        return SpatialPair(ctype1, ctype2)

    @property
    def _timesys(self):
        if self.meta.get('timesys', '').upper() == 'UTC (TBR)':
            return 'UTC'
        return super()._timesys

    @property
    def _supported_observer_coordinates(self):
        # Assume observer is at zero Stonyhurst heliographic longitude if not otherwise specified
        # https://community.openastronomy.org/t/sunpymetadatawarnings-when-using-hinode-xrt-data/393/7 for more information
        return (super()._supported_observer_coordinates
                + [(('solar_b0', 'dsun_obs'), {'lon': 0*u.deg,
                                               'lat': self.meta.get('solar_b0'),
                                               'radius': self.meta.get('dsun_obs'),
                                               'unit': (u.deg, u.deg, u.m),
                                               'frame': "heliographic_stonyhurst"})])

    @property
    def detector(self):
        return "XRT"

    @property
    def observatory(self):
        return "Hinode"

    @property
    def measurement(self):
        fw1 = self.meta.get('EC_FW1_').replace("_", " ")
        fw2 = self.meta.get('EC_FW2_').replace("_", " ")
        return f"{fw1}-{fw2}"

    @property
    def processing_level(self):
        lvl = self.meta.get('DATA_LEV', None)
        if lvl is None:
            return
        return int(lvl)

    @property
    def unit(self):
        # XRT data values are in DN and are converted into DN/s if the data has been normalized.
        # A tag starting with "XRT_RENORMALIZE" is added to the HISTORY tag in that case.
        # See Table 1.1 and Section 2.11 of the XRT Analysis Guide.
        unit = super().unit
        if not unit:
            unit = u.DN
            history = self.meta.get('HISTORY', '')
            if "xrt_renormalize" in history.lower():
                unit = u.DN / u.second
        return unit

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an XRT image"""
        return header.get('instrume') == 'XRT'


class SOTMap(GenericMap):
    """Hinode SOT Image Map definition.

    The Hinode Solar Optical Telescope (SOT) consists of a 50 cm
    diffraction-limited Gregorian telescope. It is optimized for
    accurate measurement of the vector magnetic
    field in the photosphere and dynamics of both the photosphere and
    chromosphere associated with the magnetic fields.

    Hinode was launched on 22 September 2006 into a sun-synchronous orbit.

    References
    ----------
    * `Hinode Mission Page <http://solarb.msfc.nasa.gov/index.html>`__
    * `Hinode SOT Instrument Page <http://sot.lmsal.com>`__
    * `Hinode SOT Instrument Paper <https://arxiv.org/abs/0711.1715>`__
    * `Data Analsis Guide <https://sot.lmsal.com/doc/rep/sot254/fid366/SOT00042_33_SOT_Analysis_Guide_SAG.pdf>`__
    """
    # TODO: get a link for the SOT FITS headers
    # Add in some information about the the possible instrument, observation
    # type, observable ion and wavelength

    Instruments = ['SOT/WB', 'SOT/NB', 'SOT/SP', 'SOT/CT']

    Waves = ['6302A', 'BFI no move', 'CN bandhead 3883',
             'Ca II H line', 'G band 4305', 'NFI no move',
             'TF Fe I 6302', 'TF Mg I 5172', 'TF Na I 5896',
             'blue cont 4504', 'green cont 5550', 'red cont 6684']

    Observation_Type = ['FG (simple)', 'FG focus scan',
                        'FG shuttered I and V', 'FG shutterless I and V',
                        'FG shutterless I and V with 0.2s intervals',
                        'FG shutterless Stokes', 'SP IQUV 4D array']

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname = self.detector

        # TODO (add other options, Now all treated as intensity. This follows
        # Hinode SDC archive) StokesQUV -> grey, Velocity -> EIS, Width -> EIS,
        # Mag Field Azi -> IDL 5 (STD gamma II)
        # 'WB' -> red
        # 'NB'(0 = red); (>0 = gray), # nb has 1 stokes I, the rest quv
        # 'SP' (<=1 = red); (>1 = gray) #sp has 2 stokes I, the rest quv
        color = {'SOT/WB': 'intensity',
                 'SOT/NB': 'intensity',  # For the 1st dimension
                 'SOT/SP': 'intensity',  # For the 1st 2 dimensions
                 }

        self.plot_settings['cmap'] = 'hinodesot' + color[self.instrument]

    @property
    def coordinate_system(self):
        """
        Override the default implementation to handle SOTMap-specific logic for CTYPE values.
        """
        ctype1, ctype2 = self.meta['ctype1'], self.meta['ctype2']
        if ctype1.lower() in ("solar-x", "solar_x"):
            ctype1 = 'HPLN-TAN'
        if ctype2.lower() in ("solar-y", "solar_y"):
            ctype2 = 'HPLT-TAN'
        return SpatialPair(ctype1, ctype2)

    @property
    def detector(self):
        return "SOT"

    @property
    def observatory(self):
        return "Hinode"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an SOT image."""
        return header.get('instrume') in cls.Instruments
