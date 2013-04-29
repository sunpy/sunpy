"""RHESSI Map subclass definitions"""
#pylint: disable=W0221,W0222,E1121

__author__ = "Steven Christe"
__email__ = "steven.d.christe@nasa.gov"

from sunpy.map import Map
from sunpy.cm import cm
from sunpy.time import parse_time

__all__ = ['RHESSIMap']

class RHESSIMap(Map):
    """RHESSI Image Map definition
    
    References
    ----------
    For a description of RHESSI image fits headers
    ???

    TODO: Currently (8/29/2011), cannot read fits files containing more than one 
    image (schriste)
    """
    @classmethod
    def get_properties(cls, header):
        """Parses RHESSI image header"""
        properties = Map.get_properties(header)
        
        properties.update({
            "date": parse_time(header.get('date_obs')),
            
            "detector": header.get('telescop'),
            "instrument": header.get('telescop'),
            "measurement": [header.get('energy_l'), header.get('energy_h')],
            "observatory": "SDO",
            "name": "RHESSI %d - %d keV" % (header.get('energy_l'), 
                                            header.get('energy_h')),
            "cmap": cm.get_cmap('rhessi'),
            "exposure_time": (parse_time(header.get('date_end')) - 
                              parse_time(header.get('date_obs'))).seconds,
            "coordinate_system": {
                'x': 'HPLN-TAN',
                'y': 'HPLT-TAN'
            }
        })
        return properties

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume') == 'RHESSI'