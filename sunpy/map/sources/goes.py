"""GOES SXI Map subclass definitions"""
from __future__ import absolute_import, division, absolute_import
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "Steven Christe"
__email__ = "steven.christe@nasa.gov"

from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from matplotlib.colors import LogNorm
from sunpy.map import GenericMap
from sunpy.cm import cm
from sunpy.map.sources.source_type import source_stretch
from sunpy.coordinates import get_sunearth_distance, get_sun_B0
from sunpy import sun
import astropy.units as u

__all__ = ['SXIMap']


class SXIMap(GenericMap):
    """GOES SXI Image Map

    The Solar X-ray Imager (SXI) provides uninterrupted, full-disk, soft X-ray
    solar images, with a 1 min cadence and a single-image (adjustable) dynamic
    range near 100. A set of metallic thin-film filters provides temperature
    discrimination in the 0.6–6.0 nm bandpass. The spatial resolution of
    approximately 10 arcsec FWHM is sampled with 5 arcsec pixels.

    SXI was first launched on 23 July 2001 on NOAA’s GOES-12 satellite and has
    obtained uninterrupted images starting on 22 January 2003.

    References
    ----------
    * `Mission/Instrument Page <https://sxi.ngdc.noaa.gov/index.html>`_
    * `Data notes and Fits headers <https://sxi.ngdc.noaa.gov/sxi_data_notes.html>`_
    * `Official SXI Data Book <https://sxi.ngdc.noaa.gov/images/section06.pdf>`_
    * `Instrument Paper<https://sxi.ngdc.noaa.gov/images/2005_SxiInstrumentOperationsAndData.pdf>`_
    * `SXI Performance Paper <https://sxi.ngdc.noaa.gov/images/2005_SxiPerformance.pdf>`__
    """

    def __init__(self, data, header, **kwargs):

        # the following are necessary fixes to the header
        header['WAVELNTH'] = 0
        header['WAVE_LEN'] = 0
        header['ctype1'] = "HPLN-TAN"
        header['ctype2'] = "HPLT-TAN"
        obs_time = header['date_obs']
        header['HGLT_OBS'] = get_sun_B0(obs_time)
        header['HGLN_OBS'] = 0
        header['RSUN_OBS'] = sun.solar_semidiameter_angular_size(obs_time).value
        header['RSUN_REF'] = sun.constants.radius.value
        header['DSUN_OBS'] = get_sunearth_distance(obs_time).value * u.astrophys.au
        GenericMap.__init__(self, data, header, **kwargs)

        self.meta['instrume'] = "SXI"
        self.meta['detector'] = "SXI"
        self._nickname = self.detector
        # Colour maps
        self.plot_settings['cmap'] = cm.get_cmap('sohoeit195')
        self.plot_settings['norm'] = LogNorm(vmin=1, vmax=1e4)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an GOES SXI image"""
        return header.get('INSTRUME').count('SXI')

    @property
    def measurement(self):
        """
        Returns the measurement type.
        """
        s = self.meta['WAVE_LEN']
        if s == 'WL':
            s = 'white-light'
        return s
