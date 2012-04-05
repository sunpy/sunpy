"""RHESSI Map subclass definitions"""
#pylint: disable=W0221,W0222,E1121,W0613

__author__ = "Steven Christe"
__email__ = "steven.d.christe@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.cm import cm
from sunpy.time import parse_time

class RHESSIMap(BaseMap):
    """RHESSI Image Map definition
    
    Reference
    ---------
    For a description of RHESSI image fits headers
    ???

    TODO: Currently (8/29/2011), cannot read fits files containing more than one 
    image (schriste)
    """
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)
    
    def __init__(self, data, header):
        BaseMap.__init__(self, header)
        
        self.date = parse_time(header.get('date_obs'))
        self.det = header.get('telescop')
        self.inst = header.get('telescop')
        self.meas = [header.get('energy_l'), header.get('energy_h')]
        self.name = "RHESSI %d - %d keV" % (header.get('energy_l'), 
                                            header.get('energy_h'))
        self.cmap = cm.get_cmap('rhessi')
        self.exptime = (parse_time(header.get('date_end')) - 
                        parse_time(header.get('date_obs'))).seconds

    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume') is 'RHESSI'