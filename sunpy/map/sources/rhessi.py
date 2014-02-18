"""RHESSI Map subclass definitions"""
#pylint: disable=W0221,W0222,E1121

__author__ = "Steven Christe"
__email__ = "steven.d.christe@nasa.gov"

from sunpy.map import GenericMap

__all__ = ['RHESSIMap']

class RHESSIMap(GenericMap):
    """RHESSI Image Map definition
    
    References
    ----------
    For a description of RHESSI image fits headers
    ???

    TODO: Currently (8/29/2011), cannot read fits files containing more than one 
    image (schriste)
    """
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        self._name = "RHESSI %d - %d keV" % (self.measurement[0], self.measurement[1])
        self._nickname = self.detector
        
        # Fix some broken/misapplied keywords
        if self.meta['ctype1'] == 'arcsec':
            self.meta['cunit1'] = 'arcsec'
            self.meta['ctype1'] = 'HPLN-TAN'
        if self.meta['ctype2'] == 'arcsec':
            self.meta['cunit2'] = 'arcsec'
            self.meta['ctype2'] = 'HPLT-TAN'
        
    @property
    def measurement(self):
        return [self.meta['energy_l'], self.meta['energy_h']]
    
    @property
    def detector(self):
        return self.meta['telescop']    
    
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an RHESSI image"""
        return header.get('instrume') == 'RHESSI'

