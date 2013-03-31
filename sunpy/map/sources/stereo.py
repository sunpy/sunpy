"""STEREO Map subclass definitions"""
#pylint: disable=W0221,W0222,E1121

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

from sunpy.map import Map
from sunpy.time import parse_time
from sunpy.cm import cm

__all__ = ['EUVIMap', 'CORMap']

class EUVIMap(Map):
    """EUVI Image Map definition"""
    @classmethod
    def get_properties(cls, header):
        """Parses EUVI image header"""
        properties = Map.get_properties(header)
        
        properties.update({
            "date": parse_time(header.get('date-obs',header.get('date_obs'))),
            "detector": "EUVI",
            "instrument": "SECCHI",
            "observatory": header.get('obsrvtry'),
            "cmap": cm.get_cmap('sohoeit%d' % header.get('wavelnth')),
            "nickname": "EUVI-" + header.get('obsrvtry')[-1]
        })
        return properties

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an EUVI image"""
        return header.get('detector') == 'EUVI'
        
class CORMap(Map):
    """COR Image Map definition"""
    @classmethod
    def get_properties(cls, header):
        """Parses COR image header"""
        properties = Map.get_properties(header)
        
        # @TODO: Deal with invalid values for exptime. E.g. STEREO-B COR2
        # on 2012/03/20 has -1 for some images.
        properties.update({
            "date": parse_time(header.get('date-obs',header.get('date_obs'))),
            "detector": header.get('detector'),
            "instrument": "SECCHI",
            "observatory": header.get('obsrvtry'),
            "measurement": "white-light",
            "name": "SECCHI %s" % header.get('detector'),
            "nickname": "%s-%s" % (header.get('detector'), 
                                   header.get('obsrvtry')[-1]),
            "cmap": cm.get_cmap('stereocor%s' % properties['detector'][-1])
        })
        return properties

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an COR image"""
        return header.get('detector', '').startswith('COR')