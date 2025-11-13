'''ASO-S Map subclass definitions, currently for the HXI payload is available.'''

import astropy.units as u
#from sunpy import log
from sunpy.time import parse_time
from sunpy.map.mapbase import GenericMap, SpatialPair

__author__="FuYu@PMO, from ASO-S/HXI Team"

__all__ = ['HXIMap']


class HXIMap(GenericMap):
    """
    HXI Image Maps (single map or multi-maps from one .fits file).
    
    The Hard X-ray Imager (HXI, Zhang et al. 2019b; Su et al. 2019) is one of the three payloads 
    of the advanced space-based solar observatory (ASO-S, Gan et al. 2019, 2023).
    
    The energy ranges are ~15-300 keV for imaging and ~10-300 keV for spectra (Su et al. 2024).

    ASO-S was launched on 9 October 2022 (UTC+8).

    References
    ----------
    * `ASO-S Mission Page <http://aso-s.pmo.ac.cn>`_
    * `ASO-S Data Center <http://aso-s.pmo.ac.cn/sodc/dataArchive.jsp>`_
    * `ASO-S/HXI Flare List <http://aso-s.pmo.ac.cn/hxi_flare/hxi_flare_list.html>`_
    * `The Tests and Calibrations of the Hard X-ray Imager Aboard ASO-S <https://link.springer.com/article/10.1007/s11207-024-02392-x>`_
    
    History:
    --------
        2025 08 27, FuYu@PMO, ASO/HXI Team. rsun_meters use the default value from sunpy (not fix).
             11 11, add 'fix_hxi_observer' keyword.
    """

    def __init__(self, data, header, **kwargs):
        super().__init__(data, header, **kwargs)
        asos_hxi, hxi = self.meta.get('ORIGIN').split(' ')[0], self.meta.get('INSTRUME','') #deal with different headers
        self._nickname =  asos_hxi if 'HXI' in asos_hxi else asos_hxi+'/'+hxi              #self.detector
        self.plot_settings['cmap'] = 'rhessi' #Using default cmap the same as RHESSI
        
        if kwargs.get('fix_hxi_observer',True):
            self._fix_hxi_observer()
            
    def _fix_hxi_observer(self):
        """
        Add the missing info. of 'dsun_obs','hgln_obs' and 'hglt_obs' for the fits_header, manually.
        ASO-S is located near Earth, so the observer is assumed as Earth to get the info. 
        according to the 'date_start' with sunpy.coordinates.get_earth
        
        Usage: set keyword 'fix_hxi_observer' when using sunpy.map.Map(), default:True.
               if not fix, the sunpy will fix it automatically in some cases.
        """
        if not (set(['HGLN_OBS','HGLT_OBS','DSUN_OBS']) < set(self.fits_header)):
            
            from sunpy.coordinates import get_earth
            observer = get_earth(self.date_start)
            
            if observer.frame.name=='heliographic_stonyhurst':
                self.meta['HGLN_OBS'] = observer.lon.to_value(u.deg)
                self.meta['HGLT_OBS'] = observer.lat.to_value(u.deg)
                self.meta['DSUN_OBS'] = observer.radius.to_value(u.m)
            else:
                raise ValueError("observer.frame.name not eq 'heliographic_stonyhurst', not support")

    def _get_cmap_name(self):
        """
        Using default cmap the same as RHESSI
        """
        return "rhessi" 

    def _rotation_matrix_from_crota(self):
        """
        Rotation in CROTA. HXI may not need
        """
        return super()._rotation_matrix_from_crota(crota_key='CROTA')

    @property
    def spatial_units(self):
        """
        If cunit{1 or 2} are equal to 'arcsecs', assumes that cunit{1 or 2}
        respectively are intended to be 'arcsec'.
        """
        units = [self.meta.get('cunit1', None), self.meta.get('cunit2', None)]
        if self.meta['cunit1'] == 'arcsecs':
            units[0] = 'arcsec'
        if self.meta['cunit2'] == 'arcsecs':
            units[1] = 'arcsec'
        units = [None if unit is None else u.Unit(unit.lower()) for unit in units]
        return SpatialPair(units[0], units[1])

    @property
    def coordinate_system(self):
        """
        If CTYPE{1 or 2} are equal to 'solar_x' or 'solar_y', assumes that CTYPE{1 or 2}
        respectively are intended to be 'HPLN-TAN' or 'HPLT-TAN'.
        """
        ctypes = [self.meta['ctype1'], self.meta['ctype2']]
        if ctypes[0] == 'solar_x':
            ctypes[0] = 'HPLN-TAN'
        if ctypes[1] == 'solar_y':
            ctypes[1] = 'HPLT-TAN'
        return SpatialPair(*ctypes)
    
    @property
    def date_start(self):
        """DATE-OBS is the start date"""
        date_obs = self.meta.get('DATE-OBS')
        if len(date_obs.split(" ")[0])==9:
            date_obs = date_obs[:7]+'20'+date_obs[7:]
        return parse_time(date_obs)

    @property
    def date_end(self):
        return self.date_start + self.meta.get('exptime')*u.second
    
    @property
    def _date_obs(self):
        """using date_start"""
        return self.date_start
        
    @property
    def reference_date(self):
        """using date_start, compatible with the old version of sunpy """
        return self.date_start
    
    @property
    def img_bkg(self):
        return self.meta.get('BKG')

    @property
    def img_counts(self):
        return self.meta.get('COUNTS')

    @property
    def img_algorith(self):
        return self.meta.get('ALGORITH')
    
    @property
    def waveunit(self):
        """
        If the WAVEUNIT FITS keyword is not present, defaults to keV.
        """
        return u.Unit(self.meta.get("waveunit", 'keV'))

    @property
    def wavelength(self):
        return [self.meta['energy_l'], self.meta['energy_h']]*self.waveunit

    @property
    def observatory(self):
        return self.meta.get('ORIGIN').split('/')[0]

    @property
    def instrument(self):
        #return self.meta['INSTRUME']
        origin = self.meta.get('ORIGIN')
        return origin.split(' ')[0].replace('ASO-S/','') if 'HXI' in origin else self.meta.get('INSTRUME')
    
    @property
    def detector(self):
        return self.instrument

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an HXI image"""
        #print(header.get('instrume') == 'HXI',header.get('ORIGIN'))
        #return header.get('instrume') == 'HXI'
        tag_ori = 'HXI' in header.get('ORIGIN') if header.get('ORIGIN') is not None else False
        return tag_ori or header.get('INSTRUME') == 'HXI'
        
