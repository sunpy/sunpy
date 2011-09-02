"""RHESSI Map subclass definitions"""
#pylint: disable=W0221,W0222

__author__ = "Steven Christe"
__email__ = "steven.d.christe@nasa.gov"

from sunpy.map.basemap import BaseMap
from sunpy.cm import cm
from sunpy.util import util as util

class RHESSIMap(BaseMap):
    """RHESSI Image Map definition
    
    Reference
    ---------
    For a description of RHESSI image fits headers
    ???

    TODO
    ----
    Currently (8/29/2011), cannot read fits files containing more than one image (schriste)
    """
    def __new__(cls, data, header):
        return BaseMap.__new__(cls, data)

    @classmethod
    def get_properties(cls, header):
        """Returns the default and normalized values to use for the Map"""
        properties = BaseMap.get_properties()
        properties.update({
            'date': util.anytim(header.get('date_obs')),
            'det': header.get('telescop'),
            'inst': header.get('telescop'),
            'meas': [header.get('energy_l'), header.get('energy_h')],
            'obs': header.get('telescop'),
            'name': "RHESSI " + str(header.get('energy_l')) + '-' + str(header.get('energy_h')) + ' keV',
            'cmap': cm.get_cmap(name = 'rhessi'),
            # 'norm': mpl.colors.Normalize(vmin=cls.min(), vmax=cls.max())
        })
        return properties
        
    @classmethod
    def is_datasource_for(cls, header):
        """Determines if header corresponds to an AIA image"""
        return header.get('instrume') == 'RHESSI'