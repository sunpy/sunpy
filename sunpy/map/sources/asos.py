"""
ASO-S Map subclass definitions.
"""

import astropy.units as u

from sunpy.coordinates import get_earth
from sunpy.map.mapbase import GenericMap, SpatialPair
from sunpy.time import parse_time

__all__ = ['HXIMap']


class HXIMap(GenericMap):
    """
    HXI Image Map Source.

    The Hard X-ray Imager (HXI, :cite:t:`Zhang_asos_mission_hxi_2019, Su_asos_hxi_simu_soft_2019`)
    is one of the three payloads of the Advanced Space-based Solar Observatory
    (ASO-S, :cite:t:`Gan_asos_mission_overview_2019, Gan_asos_issue_overview_2023`), which is designed
    to observe Hard X-Ray (HXR) spectra and images of solar flares. Having 91 subcollimators to
    modulate incident x-rays, HXI can obtain 91 modulation data and 45 visibilities to reconstruct
    images with a spatial resolution as high as ~3.1 arcsec.

    In addition, the energy ranges are ~15-300 keV for imaging and ~10-300 keV for spectra, more information
    please refer to :cite:t:`Su_asos_hxi_tests_calibrations_2024`.

    References
    ----------
    * `ASO-S Mission Page <http://aso-s.pmo.ac.cn>`__
    * `ASO-S Data Center <http://aso-s.pmo.ac.cn/sodc/dataArchive.jsp>`__
    * `ASO-S/HXI Flare List <http://aso-s.pmo.ac.cn/hxi_flare/hxi_flare_list.html>`__
    * The tests and calibrations paper: :cite:t:`Su_asos_hxi_tests_calibrations_2024`
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        self._nickname =  self.meta.get('nickname', 'ASO-S/HXI')
        self.plot_settings['cmap'] = 'rhessi'

    @property
    def _default_observer_coordinate(self):
        if not (set(['HGLN_OBS','HGLT_OBS','DSUN_OBS']) < set(self.fits_header)):
            return  get_earth(self.reference_date)

    def _get_cmap_name(self):
        return "rhessi"

    def _rotation_matrix_from_crota(self):
        """
        Rotation from CROTA.
        """
        return super()._rotation_matrix_from_crota(crota_key='CROTA')

    @property
    def spatial_units(self):
        units = [self.meta.get('cunit1', None), self.meta.get('cunit2', None)]
        if self.meta['cunit1'] == 'arcsecs':
            units[0] = 'arcsec'
        if self.meta['cunit2'] == 'arcsecs':
            units[1] = 'arcsec'
        units = [None if unit is None else u.Unit(unit.lower()) for unit in units]
        return SpatialPair(units[0], units[1])

    @property
    def coordinate_system(self):
        ctypes = [self.meta['ctype1'], self.meta['ctype2']]
        if ctypes[0] == 'solar_x':
            ctypes[0] = 'HPLN-TAN'
        if ctypes[1] == 'solar_y':
            ctypes[1] = 'HPLT-TAN'
        return SpatialPair(*ctypes)

    @property
    def date_start(self):
        """
        DATE-OBS is the start_date.
        Here do some thing to deal with date string like '01-May-23 13:08:16.962' or
        ' 1-May-2023 13:08:16.962' (the latter one is output from SSWIDL's map2fits)
        """
        date_obs = self.meta.get('DATE-OBS').strip()
        tmp = date_obs.split(" ")[0]
        if len(tmp)==9 and len(tmp.split("-")[-1])==2:
            date_obs = date_obs[:7]+'20'+date_obs[7:]
        return parse_time(date_obs)

    @property
    def date_end(self):
        return self.date_start + self.meta.get('exptime')*u.second

    @property
    def _date_obs(self):
        return self.date_start

    @property
    def reference_date(self):
        return self.date_start

    @property
    def waveunit(self):
        return u.Unit(self.meta.get("waveunit", 'keV'))

    @property
    def wavelength(self):
        waves = [self.meta.get('energy_l', None), self.meta.get('energy_h', None)]
        if None in waves:
            return None
        else:
            return waves*self.waveunit

    @property
    def observatory(self):
        return self.meta.get('observatory', 'ASO-S')

    @property
    def instrument(self):
        return self.meta.get('INSTRUME', 'HXI')

    @property
    def detector(self):
        return self.meta.get('detector', 'HXI')

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """
        Determines if header corresponds to a HXI image.
        """
        return 'HXI' in header.get('ORIGIN', '') or header.get('INSTRUME') == 'HXI'
