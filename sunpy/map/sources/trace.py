"""TRACE Map subclass definitions"""

__author__ = "Jack Ireland"
__email__ = "jack.ireland@nasa.gov"

import astropy.units as u
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from sunpy.map.mapbase import GenericMap, SpatialPair
from sunpy.map.sources.source_type import source_stretch

__all__ = ['TRACEMap']


class TRACEMap(GenericMap):
    """TRACE Image Map

    The Transition Region and Coronal Explorer was a
    NASA Small Explorer (SMEX) mission to image the
    solar corona and transition region at high angular and temporal resolution.
    TRACE observed the Sun in the following passbands, 5000 A, 1700 A, 1600 A,
    1550 A (C IV), 1216 A (H1 Lyman-alpha), 173 A (Fe IX), 195 A (Fe XII),
    and 284 A (Fe XV). TRACE provides solar images with an 8.5 x 8.5 arcminute
    field of view and 0.5 arcsecond pixels. It was placed in a sun-synchronous
    orbit, enabling it to make continuous solar observations.

    The TRACE mission operated was launched on 2 April 1998 and obtained its
    last science image on 6 June 2010 23:56 UT.

    References
    ----------
    * `Mission/Instrument Page <https://sdowww.lmsal.com/TRACE>`__
    * `Fits headers <https://sdowww.lmsal.com/TRACE/Project/Instrument/cal/>`__
    * `Analysis Guide <https://sdowww.lmsal.com/TRACE/tag/>`__
    * `Passband reference <https://sdowww.lmsal.com/TRACE/Project/Instrument/inspass.htm>`__

    .. note::

        Note that this map definition is currently only being tested on JPEG2000
        files. TRACE FITS data is stored in a more complex format. Typically
        TRACE data is stored in hourly "tri" files that store all the data taken
        by TRACE in the hour indicated by the filename. Those files must first be
        understood and parsed to obtain the science data. The ability to do this
        is not yet in SunPy, but is available in SSWIDL. Please refer to the links
        above concerning how to read "tri" files in SSWIDL.
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)

        self._nickname = self.detector
        # Colour maps
        self.plot_settings['cmap'] = 'trace' + str(self.meta['WAVE_LEN'])
        self.plot_settings['norm'] = ImageNormalize(
            stretch=source_stretch(self.meta, LogStretch()), clip=False)

    @property
    def coordinate_system(self):
        """
        Override the default implementation to handle TRACEMAP-specific logic for CTYPE values.
        """
        ctype1, ctype2 = self.meta['ctype1'], self.meta['ctype2']
        if ctype1.lower() in ("solar-x", "solar_x"):
            ctype1 = 'HPLN-TAN'
        if ctype2.lower() in ("solar-y", "solar_y"):
            ctype2 = 'HPLT-TAN'
        return SpatialPair(ctype1, ctype2)

    @property
    def spatial_units(self):
        """
        If not present in CUNIT{1,2} keywords, defaults to arcsec.
        """
        return SpatialPair(u.Unit(self.meta.get('cunit1', 'arcsec')),
                           u.Unit(self.meta.get('cunit2', 'arcsec')))

    @property
    def observatory(self):
        return "TRACE"

    @property
    def detector(self):
        return "TRACE"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an TRACE image"""
        return header.get('instrume') == 'TRACE'

    @property
    def measurement(self):
        s = self.meta['WAVE_LEN']
        if s == 'WL':
            s = 'white-light'
        return s
