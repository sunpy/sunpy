"""STEREO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map import GenericMap
from sunpy.cm import cm

__all__ = ['EUVIMap', 'CORMap']

class EUVIMap(GenericMap):
    """EUVI Image Map definition"""
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        self._name = self.observatory + " " + self.detector + " " + str(self.measurement)
        self._nickname = "{0}-{1}".format(self.detector, self.observatory[-1])
        
        self.cmap = cm.get_cmap('sohoeit%d' % self.wavelength)

        # Try to identify when the FITS meta data ddes not have the correct
        # date FITS keyword
        if ('date_obs' in self.meta) and not('date-obs' in self.meta):
            self.meta['date-obs'] = self.meta['date_obs']

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an EUVI image"""
        return header.get('detector') == 'EUVI'
        
class CORMap(GenericMap):
    """COR Image Map definition"""
    
    def __init__(self, data, header, **kwargs):
        
        GenericMap.__init__(self, data, header, **kwargs)
        
        self._name = self.observatory + " " + self.detector + " " + str(self.measurement)
        self._nickname = "{0}-{1}".format(self.detector, self.observatory[-1])
        
        self.cmap = cm.get_cmap('stereocor%s' % self.detector[-1])

        # Try to identify when the FITS meta data ddes not have the correct
        # date FITS keyword
        if ('date_obs' in self.meta) and not('date-obs' in self.meta):
            self.meta['date-obs'] = self.meta['date_obs']
        
    @property
    def measurement(self):
        # TODO: This needs to do more than white-light.  Should give B, pB, etc.
        return "white-light"

    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        """Determines if header corresponds to an COR image"""
        return header.get('detector', '').startswith('COR')
        
