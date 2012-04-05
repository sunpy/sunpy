"""STEREO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1121,W0613

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.time import parse_time
from sunpy.cm import cm

class EUVIMap(BaseMap):
    """EUVI Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
    
    def __init__(self, data, header):
        BaseMap.__init__(self, header)
        self.date = parse_time(header.get('date_obs'))
        self.det = "EUVI"
        self.inst = "SECCHI"
        self.obs = header.get('obsrvtry')
        self.name = "EUVI %s" % self.meas
        self.cmap = cm.get_cmap('sohoeit%d' % header.get('wavelnth'))

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an EUVI image"""
        return header.get('detector') is 'EUVI'
        
class CORMap(BaseMap):
    """COR Image Map definition"""
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
    
    # @TODO: Deal with invalid values for exptime. E.g. STEREO-B COR2
    # on 2012/03/20 has -1 for some images.
    def __init__(self, data, header):
        BaseMap.__init__(self, header)
        self.date = parse_time(header.get('date_obs'))
        self.det = header.get('detector')
        self.inst = "SECCHI"
        self.obs = header.get('obsrvtry')
        self.meas = "white-light"
        self.name = "SECCHI %s" % header.get('detector')
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an COR image"""
        return header.get('detector', '').startswith('COR')

