"""GOES SXI Map subclass definitions"""
from __future__ import absolute_import, division, absolute_import
#pylint: disable=W0221,W0222,E1101,E1121

__author__ = "[Iain Hannah, Steven Christe]"
__email__ = "steven.christe@nasa.gov"

from matplotlib.colors import LogNorm

import astropy.units as u

from sunpy.map import GenericMap
from sunpy.cm import cm
from sunpy.coordinates import get_sunearth_distance, get_sun_B0
from sunpy import sun

__all__ = ['SXIMap']


class SXIMap(GenericMap):
    """GOES SXI Image Map

    The GOES Solar X-ray Imager (SXI) provides uninterrupted, full-disk, soft X-ray
    solar images, with a 1 min cadence and a single-image (adjustable) dynamic
    range near 100. A set of metallic thin-film filters provides temperature
    discrimination in the 0.6–6.0 nm bandpass. The spatial resolution of
    approximately 10 arcsec Full width half maximum (FWHM) is sampled with 5
    arcsec pixels.

    SXI was first launched on 23 July 2001 on NOAA’s GOES-12 satellite and has
    obtained uninterrupted images starting on 22 January 2003.

    References
    ----------
    * `Mission Instrument Page <https://sxi.ngdc.noaa.gov/index.html>`_
    * `Data notes and Fits headers <https://sxi.ngdc.noaa.gov/sxi_data_notes.html>`_
    * `Official SXI Data Book <https://sxi.ngdc.noaa.gov/images/section06.pdf>`_
    * `Instrument Paper DOI: 10.1007/s11207-005-7416-x <https://link.springer.com/article/10.1007/s11207-005-7416-x>`_
    * `SXI Performance Paper DOI: 10.1007/s11207-005-7417-9 <https://link.springer.com/article/10.1007/s11207-005-7417-9>`_
    """

    filter_wheel_measurements = ["OPEN", "PTHNA", "PMEDA", "PTHK",
                                 "BTHNA", "BMED", "BTHK", "UV", "RDSH"]

    def __init__(self, data, header, **kwargs):

        # check if the filter setting is in the list of expected values
        fw = header.get('WAVELNTH')
        if fw.lower().strip() not in [item.lower() for item in self.filter_wheel_measurements]:
            raise ValueError('Unexpected filter setting in header.')

        # the following are necessary fixes to the header
        header['DETECTOR'] = 'SXI'
        # save the filter setting
        header['FILTER'] = header.get('WAVELNTH')
        # update the wavelength with a number since a str is not supported.
        # TODO: change the number according to the filter
        header['WAVELNTH'] = 14
        header['WAVEUNIT'] = 'Angstrom'
        # need to switch these as they are not correct - Hannah
        crota2 = header['CROTA1']
        crota1 = header['CROTA2']
        header['CROTA1'] = crota1
        header['CROTA2'] = crota2
        header['CTYPE1'] = "HPLN-TAN"
        header['CTYPE2'] = "HPLT-TAN"
        obs_time = header['date_obs']
        header['HGLT_OBS'] = get_sun_B0(obs_time).to('deg').value
        header['HGLN_OBS'] = 0
        header['RSUN_OBS'] = sun.solar_semidiameter_angular_size(obs_time).to('arcsec').value
        header['RSUN_REF'] = sun.constants.radius.to('m').value
        header['DSUN_OBS'] = get_sunearth_distance(obs_time).to('m').value

        GenericMap.__init__(self, data, header, **kwargs)

        self._nickname = "{0}/{1}".format(self.meta['telescop'], self.detector)
        # Colour maps
        self.plot_settings['cmap'] = cm.get_cmap('sohoeit195')
        self.plot_settings['norm'] = LogNorm(vmin=1, vmax=1e4)

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an GOES SXI image"""
        return header.get('instrume').count('SXI')

    @property
    def measurement(self):
        """
        Returns the measurement type.
        """
        s = u.Quantity(self.meta['WAVELNTH'], self.meta['WAVEUNIT'])
        return s
